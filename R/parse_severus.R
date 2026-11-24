suppressPackageStartupMessages({
  library(data.table)
})

# Severus writes both mates as separate records and tags single breakends "sBND"; both are collapsed here so downstream sees one row per rearrangement
BND_SVTYPES = c("BND", "sBND")

# BND-derived svclass values (two loci, no span/size); BND_CIRCOS_CLASSES is this set minus "single breakend"
SV_JUNCTION_CLASSES = c("translocation", "intra-chr breakend", "single breakend")

# Column contract of both SV table builders (the section adds `caller`, and Padfoot columns when
# a Padfoot file exists). Downstream (overview, panel matching, circos) reads these BY NAME:
# append, never rename.
SV_TABLE_COLS = c("id", "svclass", "svtype", "chrom_a", "pos_a", "gene_a",
                  "chrom_b", "pos_b", "gene_b", "sv_len", "vaf",
                  "consequence", "impact", "support", "svtype_raw",
                  "filter", "wl_rescue", "orient_a", "orient_b")

# FILTER values that mean the record passed (VCF "." = no filters applied)
SV_PASS_FILTERS = c("PASS", ".")

# Whitelist-rescued rows (Severus INFO/WL_RESCUE): kept only because they lie in the
# --severus_whitelist regions. FALSE for every row of a table without the column
sv_is_rescued = function(t) {
  if (is.null(t) || nrow(t) == 0 || !"wl_rescue" %in% names(t)) return(rep(FALSE, NROW(t)))
  !is.na(t$wl_rescue) & nzchar(t$wl_rescue)
}

# Panel-matching windows: single source of truth for sv_panel_hits() and the client-side filter
SV_PANEL_WINDOW_BND   = 1e6   # distance from either breakend of a BND
SV_PANEL_WINDOW_OTHER = 1e5   # distance from the span of a DEL/DUP/INV/INS

# SV ring palette and y-positions. These live here, not in R/circos.R, because
# severus_circos_tracks() below is what stamps them onto every row and the tests source this
# file without circlize. There must be exactly one definition of each: draw_circos() labels the
# SV ring's y-axis from SV_YPOS, so a second copy drifting from this one silently mislabels the
# axis, and --circos-sv-* in assets/styles/report.scss mirrors SV_COLOURS for the HTML legend,
# so a second copy there renders the legend in colours the plot never used.
SV_COLOURS = c(INS = "#f97e02", DEL = "#020272", INV = "#e7cc02", DUP = "#e41a1c")
SV_YPOS    = c(INS = 1.0,       DEL = 0.66,      INV = 0.33,      DUP = 0.05)

# Read a VCF's data records; empty data.table for a header-only VCF
.severus_read_vcf = function(vcf_file, col_names) {
  con = gzfile(vcf_file, "rb")
  skip_n = 0L
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#CHROM")) break
    skip_n = skip_n + 1L
  }
  close(con)

  # fread() errors when skip lands on the last line (no records)
  dt = tryCatch(
    fread(vcf_file, skip = skip_n + 1L, sep = "\t", header = FALSE),
    error = function(e) data.table()
  )
  if (nrow(dt) == 0) return(data.table())

  # Name columns positionally: a VCF without FORMAT/SAMPLE columns must not read as zero SVs
  n = min(ncol(dt), length(col_names))
  setnames(dt, seq_len(n), col_names[seq_len(n)])
  for (nm in setdiff(col_names, names(dt))) dt[, (nm) := NA_character_]
  dt[, ..col_names]
}

# One INFO sub-field per record, NA where absent (length-preserving, unlike bare regmatches())
.info_val = function(info_vec, key) {
  m   = regexpr(paste0("(?:^|;)", key, "=([^;]+)"), info_vec, perl = TRUE)
  out = rep(NA_character_, length(info_vec))
  hit = m > 0
  if (any(hit)) out[hit] = sub(".*=", "", regmatches(info_vec, m))
  out
}

# Partner locus from ALT bracket notation (all four forms, with or without "chr"); NA for a single breakend
.bnd_partner = function(alt) {
  pat = "^.*?[\\[\\]]([^\\[\\]:]+):([0-9]+)[\\[\\]].*$"
  has = grepl(pat, alt, perl = TRUE)
  list(
    chrom = ensure_chr_prefix(ifelse(has, sub(pat, "\\1", alt, perl = TRUE), NA_character_)),
    pos   = ifelse(has, suppressWarnings(as.integer(sub(pat, "\\2", alt, perl = TRUE))),
                   NA_integer_)
  )
}

# Breakend orientation of each side, in the HMF/LINX convention: +1 when the retained
# sequence lies left of (ends at) the breakend, -1 when it lies right of (starts at) it.
# From the record's own ALT for a BND (VCF 4.2 bracket notation: `t[p[`/`t]p]` keep the
# left side here, `]p]t`/`[p[t` the right; the partner keeps its left side under `]`, its
# right under `[`) and for a single breakend (`t.` / `.t`); from the type for a symbolic
# DEL (+1/-1) and DUP (-1/+1); from INFO/STRANDS for an INV, only when it states one
# orientation for both sides (`++`/`--`). NA otherwise (INS, unknown ALT forms)
.sv_orientation = function(alt, svtype, strands = NULL) {
  alt = as.character(alt); n = length(alt)
  own = rep(NA_integer_, n); partner = rep(NA_integer_, n)
  bracket = grepl("[][]", alt)
  own[bracket & grepl("^[^][]", alt)]  = 1L
  own[bracket & grepl("[^][]$", alt)]  = -1L
  partner[bracket & grepl("]", alt, fixed = TRUE)] = 1L
  partner[bracket & grepl("[", alt, fixed = TRUE)] = -1L
  single = !bracket & grepl("^\\.[A-Za-z]+$|^[A-Za-z]+\\.$", alt)
  own[single & endsWith(alt, ".")]   = 1L
  own[single & startsWith(alt, ".")] = -1L
  sym = !bracket & !single
  own[sym & svtype == "DEL"] = 1L;  partner[sym & svtype == "DEL"] = -1L
  own[sym & svtype == "DUP"] = -1L; partner[sym & svtype == "DUP"] = 1L
  if (!is.null(strands)) {
    inv = sym & svtype == "INV" & strands %in% c("++", "--")
    own[inv] = partner[inv] = ifelse(strands[inv] == "++", 1L, -1L)
  }
  list(own = own, partner = partner)
}

# Per-caller VCF dialect. Both callers write one record per breakend and link mates by an
# INFO key, but they disagree on its name, on where the allele fraction and read support live,
# and on whether a deletion is a DEL record at all (SAVANA writes every rearrangement as a BND
# pair and states the junction orientation in BP_NOTATION instead).
#   mate_key        INFO key naming the mate record
#   vaf_format      FORMAT key for the allele fraction, or NA
#   vaf_info        INFO key for it, or NA (SAVANA's TUMOUR_AF is Number=2: this record's own
#                   breakpoint first, its mate's second — the first value is taken)
#   support_format  FORMAT key for supporting reads (Severus DV), or NA
#   support_info    INFO key for it (SAVANA TUMOUR_READ_SUPPORT), or NA
#   orientation     INFO key typing same-chromosome BND pairs (see SV_ORIENTATION_TYPES), or NA
SV_CALLER_DIALECTS = list(
  severus = list(mate_key = "MATE_ID", vaf_format = "VAF", vaf_info = NA_character_,
                 support_format = "DV", support_info = NA_character_,
                 orientation = NA_character_),
  savana  = list(mate_key = "MATEID", vaf_format = NA_character_, vaf_info = "TUMOUR_AF",
                 support_format = NA_character_, support_info = "TUMOUR_READ_SUPPORT",
                 orientation = "BP_NOTATION")
)

# BP_NOTATION of a same-chromosome BND pair -> the contiguous type it is. Without this every
# SAVANA deletion is an "intra-chr breakend": a junction with 1 Mb panel windows, an arc on the
# breakend circos and no size (~1,150 of P15's 1,246 SAVANA pairs are small +- deletions).
SV_ORIENTATION_TYPES = c("+-" = "DEL", "-+" = "DUP", "++" = "INV", "--" = "INV")

# FORMAT/<key> of the first sample column, format-group by format-group as in parse_caller_vcf()
.sv_format_val = function(dt, key) {
  out = rep(NA_real_, nrow(dt))
  if (is.na(key)) return(out)
  for (fmt in unique(dt$FORMAT[!is.na(dt$FORMAT)])) {
    idx_rows = which(dt$FORMAT == fmt)
    k = match(key, strsplit(fmt, ":", fixed = TRUE)[[1]])
    if (is.na(k)) next
    split_s = strsplit(as.character(dt$SAMPLE1[idx_rows]), ":", fixed = TRUE)
    out[idx_rows] = vapply(split_s, function(x)
      if (length(x) >= k) suppressWarnings(as.numeric(x[k])) else NA_real_,
      numeric(1))
  }
  out
}

# First value of a (possibly Number=2) numeric INFO key
.sv_info_num = function(info, key) {
  if (is.na(key)) return(rep(NA_real_, length(info)))
  suppressWarnings(as.numeric(sub(",.*$", "", .info_val(info, key))))
}

# Severus' parser, kept by name: the section, the tests and the overview all call it
parse_severus_somatic_records = function(vcf_file) parse_sv_vcf(vcf_file, "severus")

# Parse a somatic SV VCF into one row per rearrangement:
#   id, id_b, svtype, svtype_raw, svclass, chrom_a, pos_a, chrom_b, pos_b, sv_len, vaf, support,
#   filter, wl_rescue
# `caller` picks the dialect (SV_CALLER_DIALECTS). attr(, "n_records") is the VCF's record count
# before mate collapsing, for the footnote. A pair is rescued when either mate carries WL_RESCUE,
# and passes only when both do; a non-PASS record that was not rescued is dropped, as small
# variants are, and counted in attr(, "n_failed").
parse_sv_vcf = function(vcf_file, caller = "severus") {
  if (is.null(vcf_file) || !file.exists(vcf_file)) return(data.table())
  d = SV_CALLER_DIALECTS[[caller]]
  if (is.null(d)) stop("parse_sv_vcf(): no VCF dialect for caller '", caller, "'")

  dt = .severus_read_vcf(vcf_file, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                                     "FILTER", "INFO", "FORMAT", "SAMPLE1"))
  if (nrow(dt) == 0) return(data.table())
  n_records = nrow(dt)

  dt[, ID := as.character(ID)]
  dt[, CHROM := ensure_chr_prefix(CHROM)]
  dt[, SVTYPE := .info_val(INFO, "SVTYPE")]
  dt[, is_bnd := SVTYPE %in% BND_SVTYPES]

  # END and SVLEN are optional; create the columns unconditionally
  dt[, `:=`(END = NA_integer_, SVLEN = NA_integer_)]
  dt[grepl("END=",   INFO, fixed = TRUE), END   := as.integer(.info_val(INFO, "END"))]
  dt[grepl("SVLEN=", INFO, fixed = TRUE), SVLEN := as.integer(.info_val(INFO, "SVLEN"))]
  dt[is.na(END), END := POS]

  # Partner locus: from ALT for a BND, the SV's own end otherwise
  partner = .bnd_partner(dt$ALT)
  dt[, `:=`(chrom_b = fifelse(is_bnd, partner$chrom, CHROM),
            pos_b   = fifelse(is_bnd, partner$pos,   as.integer(END)))]

  dt[, mate_id := .info_val(INFO, d$mate_key)]

  # Orientation of this record's side and its partner's (.sv_orientation())
  ori = .sv_orientation(dt$ALT, dt$SVTYPE, .info_val(dt$INFO, "STRANDS"))
  dt[, `:=`(orient_a = ori$own, orient_b = ori$partner)]

  # Group mates on the unordered {ID, MATE_ID} pair (fallback for a breakend: shared `_1`/`_2`
  # stem); non-pairs stay as singletons. The stem fallback is BND-only because SAVANA names its
  # unpaired INS records `ID_n_1` too.
  dt[, pair_key := fifelse(!is.na(mate_id),
                           paste(pmin(ID, mate_id), pmax(ID, mate_id), sep = "|"),
                           fifelse(is_bnd, sub("_[12]$", "", ID), ID))]
  dt[, side_rank := fifelse(grepl("_1$", ID), 1L, 2L)]
  dt[, n_in_pair := .N, by = pair_key]
  if (any(dt$n_in_pair > 2L)) {
    message(caller, " VCF: ", sum(dt$n_in_pair > 2L),
            " records share a mate group with more than two members; kept unpaired.")
    dt[n_in_pair > 2L, pair_key := ID]
    dt[, n_in_pair := .N, by = pair_key]
  }

  dt[, FILTER := fifelse(is.na(FILTER) | !nzchar(FILTER), ".", as.character(FILTER))]
  dt[, wl := .info_val(INFO, "WL_RESCUE")]
  dt[, `:=`(pair_filter = { f = unique(FILTER[!FILTER %in% SV_PASS_FILTERS])
                            if (length(f)) paste(f, collapse = ";") else FILTER[1] },
            pair_wl = { w = wl[!is.na(wl)]; if (length(w)) w[1] else NA_character_ }),
     by = pair_key]
  failed = !dt$pair_filter %in% SV_PASS_FILTERS & is.na(dt$pair_wl)
  n_failed = sum(failed)
  if (n_failed > 0) dt = dt[!failed]
  if (nrow(dt) == 0) {
    out = data.table()
    setattr(out, "n_records", n_records)
    setattr(out, "n_failed", n_failed)
    return(out)
  }
  setorder(dt, pair_key, side_rank, ID)
  dt = dt[dt[, .I[1L], by = pair_key]$V1]

  # The mate's record ID, which is what the per-side VEP annotation is keyed on.
  dt[, id_b := fifelse(n_in_pair == 2L & !is.na(mate_id), mate_id, NA_character_)]
  dt[, svtype_raw := SVTYPE]

  # Orientation typing (SAVANA): a same-chromosome mate pair becomes the contiguous type its
  # BP_NOTATION describes and from here on is a span exactly like a Severus DEL/DUP/INV —
  # svclass, panel windows, size, circos ring. svtype_raw keeps "BND". Interchromosomal pairs
  # stay translocations; an orphaned record (mate filtered out) is not typed.
  dt[, typed := FALSE]
  if (!is.na(d$orientation)) {
    dt[, orient := .info_val(INFO, d$orientation)]
    dt[, typed := is_bnd & n_in_pair == 2L & !is.na(id_b) & !is.na(chrom_b) &
                  chrom_b == CHROM & !is.na(pos_b) & orient %in% names(SV_ORIENTATION_TYPES)]
    if (any(dt$typed)) {
      dt[typed == TRUE, SVTYPE := unname(SV_ORIENTATION_TYPES[orient])]
      # A span runs left to right (sv_display_columns() shows pos_b <= pos_a as a point);
      # the ids swap with the loci so each side keeps its own VEP record
      dt[typed == TRUE & pos_b < POS,
         c("POS", "pos_b", "ID", "id_b", "orient_a", "orient_b") :=
           .(pos_b, POS, id_b, ID, orient_b, orient_a)]
      dt[typed == TRUE, is_bnd := FALSE]
    }
  }

  dt[, svclass := fcase(
    is_bnd & is.na(chrom_b),      "single breakend",
    is_bnd & chrom_b != CHROM,    "translocation",
    is_bnd,                       "intra-chr breakend",
    default = SVTYPE
  )]

  # A typed pair's length is SVLEN where the caller wrote a positive one, else its span; a
  # junction has none (SAVANA writes SVLEN=0 on a translocation)
  dt[typed == TRUE & (is.na(SVLEN) | SVLEN <= 0L), SVLEN := as.integer(abs(pos_b - POS))]
  dt[svclass %in% SV_JUNCTION_CLASSES & !is.na(SVLEN) & SVLEN == 0L, SVLEN := NA_integer_]

  vaf = .sv_format_val(dt, d$vaf_format)
  if (!is.na(d$vaf_info)) vaf = fifelse(is.na(vaf), .sv_info_num(dt$INFO, d$vaf_info), vaf)
  support = .sv_format_val(dt, d$support_format)
  if (!is.na(d$support_info))
    support = fifelse(is.na(support), .sv_info_num(dt$INFO, d$support_info), support)
  dt[, `:=`(VAF = vaf, SUPPORT = as.integer(round(support)))]

  out = dt[, .(id = ID, id_b, svtype = SVTYPE, svtype_raw, svclass,
               chrom_a = CHROM, pos_a = POS, chrom_b, pos_b,
               sv_len = SVLEN, vaf = VAF, support = SUPPORT,
               filter = pair_filter, wl_rescue = pair_wl, orient_a, orient_b)]
  setattr(out, "n_records", n_records)
  setattr(out, "n_failed", n_failed)
  out
}

# Circos tracks from the collapsed records: BND links and a non-BND track with colour/y-position
severus_circos_tracks = function(records) {
  empty = list(translocations = data.table(), nontrans = data.table())
  if (is.null(records) || nrow(records) == 0) return(empty)

  is_bnd = records$svtype %in% BND_SVTYPES

  trans = records[is_bnd & !is.na(chrom_b) & !is.na(pos_b),
                  .(chrom = chrom_a, pos = pos_a, chrom2 = chrom_b, pos2 = pos_b)]
  # Dedup on the unordered endpoint pair so a VCF without MATE_ID cannot draw an arc twice
  if (nrow(trans) > 0) {
    a = paste(trans$chrom, trans$pos, sep = ":")
    b = paste(trans$chrom2, trans$pos2, sep = ":")
    trans = unique(trans[, .link_key := paste(pmin(a, b), pmax(a, b))],
                   by = ".link_key")[, .link_key := NULL]
  }

  nontrans = records[!is_bnd,
    .(chrom = chrom_a, pos = pos_a, end = pos_b, svtype, svlen = sv_len,
      circos_pos = unname(SV_YPOS[svtype]), circos_col = unname(SV_COLOURS[svtype]))]

  list(translocations = trans, nontrans = nontrans)
}

# SV display table: mate-collapsed records plus per-breakend VEP symbols, joined on record ID (locus only as fallback); panel matching is on coordinates, see sv_panel_hits()
build_sv_table_from_vep = function(records, vep_sv_vcf) {
  if (is.character(records)) records = parse_severus_somatic_records(records)
  if (is.null(records) || nrow(records) == 0) return(data.table())

  sv = copy(records)
  out_cols = SV_TABLE_COLS
  # Hand-built records (tests, older callers) may predate the support/svtype_raw columns
  if (!"support" %in% names(sv))    sv[, support := NA_integer_]
  if (!"svtype_raw" %in% names(sv)) sv[, svtype_raw := svtype]
  if (!"id_b" %in% names(sv))       sv[, id_b := NA_character_]
  if (!"filter" %in% names(sv))     sv[, filter := NA_character_]
  if (!"wl_rescue" %in% names(sv))  sv[, wl_rescue := NA_character_]
  if (!"orient_a" %in% names(sv))   sv[, orient_a := NA_integer_]
  if (!"orient_b" %in% names(sv))   sv[, orient_b := NA_integer_]

  vep = parse_vep_vcf(vep_sv_vcf)
  if (is.null(vep) || nrow(vep) == 0) {
    sv[, `:=`(gene_a = NA_character_, gene_b = NA_character_,
              consequence = NA_character_, impact = NA_character_)]
    out = sv[, ..out_cols]
    setattr(out, "n_annotated", 0L)
    return(out)
  }

  # Highest-impact annotation first, so annot[1] below is the one that survives.
  impact_rank = c(HIGH = 1L, MODERATE = 2L, LOW = 3L, MODIFIER = 4L)
  vep[, impact_rank := impact_rank[impact]]
  vep[is.na(impact_rank), impact_rank := 5L]
  setorder(vep, impact_rank)

  collapse_genes = function(x) {
    g = unique(x[!is.na(x) & nzchar(x)])
    if (length(g) == 0) NA_character_ else paste(g, collapse = ",")
  }

  by_id = if ("id" %in% names(vep) && any(!is.na(vep$id))) {
    vep[!is.na(id) & nzchar(id), .(gene = collapse_genes(symbol),
                                   consequence = consequence[1], impact = impact[1]),
        by = .(id)]
  } else NULL

  by_locus = vep[, .(gene = collapse_genes(symbol),
                     consequence = consequence[1], impact = impact[1]),
                 by = .(chrom, pos)]

  # Key choice is made once for the whole table, not per side
  use_id = !is.null(by_id) && any(c(sv$id, sv$id_b) %in% by_id$id)

  # Records filtered out on the annotation side simply get no symbols
  side_annot = function(ids, chroms, positions) {
    src = if (use_id) by_id else by_locus
    m   = if (use_id) match(ids, by_id$id)
          else match(paste(chroms, positions), paste(by_locus$chrom, by_locus$pos))
    res = data.table(gene = NA_character_, consequence = NA_character_,
                     impact = NA_character_)[rep(1L, length(chroms))]
    hit = !is.na(m)
    if (any(hit)) res[hit, `:=`(gene = src$gene[m[hit]],
                                consequence = src$consequence[m[hit]],
                                impact = src$impact[m[hit]])]
    res
  }

  a = side_annot(sv$id,   sv$chrom_a, sv$pos_a)
  b = side_annot(sv$id_b, sv$chrom_b, sv$pos_b)

  sv[, `:=`(gene_a = a$gene, gene_b = b$gene)]
  # Consequence/impact describe the rearrangement, taken from the higher-impact side.
  rank_of = function(x) { r = unname(impact_rank[x]); r[is.na(r)] = 5L; r }
  use_b   = rank_of(b$impact) < rank_of(a$impact)
  sv[, `:=`(consequence = fifelse(use_b, b$consequence, a$consequence),
            impact      = fifelse(use_b, b$impact,      a$impact))]

  # Join coverage for the footnote: rows with an annotation record on either side
  src_ids = if (use_id) by_id$id else NULL
  n_annot = if (use_id) sum(sv$id %in% src_ids | (!is.na(sv$id_b) & sv$id_b %in% src_ids))
            else sum(paste(sv$chrom_a, sv$pos_a) %in% paste(by_locus$chrom, by_locus$pos))
  out = sv[, ..out_cols]
  setattr(out, "n_annotated", as.integer(n_annot))
  out
}

# Readable locus/size/genes columns: span (DEL/DUP/INV/INS) vs junction (BND) is decided by svclass, not by chromosome; returns locus, locus_sort, size, size_bp, genes
sv_display_columns = function(sv_table, chrom_levels = NULL) {
  empty = data.table(locus = character(), locus_sort = character(),
                     size = character(), size_bp = numeric(), genes = character())
  if (is.null(sv_table) || nrow(sv_table) == 0) return(empty)

  d = as.data.table(sv_table)
  need = c("svclass", "chrom_a", "pos_a", "chrom_b", "pos_b", "sv_len",
           "gene_a", "gene_b")
  absent = setdiff(need, names(d))
  if (length(absent) > 0)
    stop("sv_display_columns(): sv_table is missing ", paste(absent, collapse = ", "))

  is_junction = d$svclass %in% SV_JUNCTION_CLASSES
  pa = suppressWarnings(as.numeric(d$pos_a))
  pb = suppressWarnings(as.numeric(d$pos_b))

  # formatC(), not format(): format() pads to a common width
  bp = function(x) formatC(x, format = "d", big.mark = ",")
  at = function(chrom, pos) paste0(chrom, ":", bp(pos))

  locus = fcase(
    is.na(pa) | is.na(d$chrom_a),   NA_character_,
    d$svclass == "single breakend", paste0(at(d$chrom_a, pa), " (unpaired)"),
    is_junction,
      # An arrow, not a dash: the two loci are joined, not a span; both sides named so a per-column search matches either
      fifelse(is.na(pb) | is.na(d$chrom_b),
              paste0(at(d$chrom_a, pa), " (unpaired)"),
              paste0(at(d$chrom_a, pa),
                     fifelse(d$svclass == "translocation", " → ", " ↔ "),
                     at(d$chrom_b, pb))),
    # A contiguous type: one span. An INS has END == POS, so it is a single point.
    default = fifelse(is.na(pb) | pb <= pa,
                      at(d$chrom_a, pa),
                      paste0(at(d$chrom_a, pa), "–", bp(pb)))
  )

  # SVLEN where Severus wrote one — an INS's length is *not* its span — else the span.
  len = suppressWarnings(as.numeric(d$sv_len))
  size_bp = fifelse(is_junction, NA_real_, fifelse(!is.na(len), abs(len), pb - pa))
  # as.character(): fmt_bp() is ifelse()-based and returns logical on an all-NA size_bp
  size = as.character(fmt_bp(size_bp))

  split_genes = function(x) {
    if (is.na(x) || !nzchar(x)) return(character())
    g = trimws(strsplit(x, ",", fixed = TRUE)[[1]])
    unique(g[nzchar(g)])
  }
  # Middot rather than spaces: DT escapes &nbsp; and HTML collapses whitespace
  genes = vapply(seq_len(nrow(d)), function(i) {
    a = split_genes(d$gene_a[i]); b = split_genes(d$gene_b[i])
    if (!is_junction[i]) return(paste(unique(c(a, b)), collapse = ", "))
    sides = c(if (length(a)) paste0("A: ", paste(a, collapse = ",")),
              if (length(b)) paste0("B: ", paste(b, collapse = ",")))
    paste(sides, collapse = " · ")
  }, character(1))

  rank = if (length(chrom_levels)) match(d$chrom_a, chrom_levels)
         else rep(NA_integer_, nrow(d))
  rank[is.na(rank)] = 999L
  locus_sort = fifelse(is.na(pa), "", sprintf("%03d:%011.0f", rank, pa))

  data.table(locus      = fifelse(is.na(locus), "", locus),
             locus_sort = locus_sort,
             size       = fifelse(is.na(size), "", size),
             size_bp    = size_bp,
             genes      = genes)
}

# Panel hits per SV row as "GENE (side, how)" labels ("" = none): coordinate panels match within the windows, symbol-only panels on VEP symbols; several panels union and add a " [name]" suffix. Mirrored by svPanelHits() in templates/per_sample.qmd. Built on sv_panel_hit_table(), the long form
sv_panel_hits = function(sv_table, panels,
                         bnd_window = SV_PANEL_WINDOW_BND,
                         other_window = SV_PANEL_WINDOW_OTHER) {
  n = if (is.null(sv_table)) 0L else nrow(sv_table)
  if (n == 0) return(character(0))

  ps = .as_panel_list(panels)
  if (length(ps) == 0) return(rep("", n))

  ht = sv_panel_hit_table(sv_table, ps, bnd_window, other_window)
  out = rep("", n)
  if (nrow(ht) == 0) return(out)

  lab = .sv_hit_labels(ht, tagged = length(ps) > 1)
  per_row = lab[, .(hit = paste(unique(label), collapse = ", ")), by = row]
  out[per_row$row] = per_row$hit
  out
}

# A panel object and a list of panels are both plain lists; tell them apart by load_gene_panel()'s fields
.as_panel_list = function(panels) {
  if (is.null(panels) || length(panels) == 0) return(list())
  if (!is.null(panels$genes) || !is.null(panels$has_coords)) {
    nm = if (!is.null(panels$name)) as.character(panels$name)[1] else "panel"
    return(setNames(list(panels), nm))
  }
  ps = panels[!vapply(panels, is.null, logical(1))]
  if (is.null(names(ps))) names(ps) = paste0("panel", seq_along(ps))
  ps
}

# Long-form SV panel hits: one row per (sv row, gene, side, gap, panel), in label order.
#   row     index into sv_table
#   id      sv_table$id[row] (NA when the table has no id column)
#   gene    panel gene (coordinate mode: the panel's interval_gene as written; symbol mode:
#           the upper-cased VEP symbol)
#   side    "A" | "B" (either breakend of a BND, or the annotated side in symbol mode) | "span"
#   gap     0 when the locus overlaps the gene, else bp distance from the unpadded locus
#   direct  gap == 0 (always TRUE in symbol mode: a VEP symbol sits on the breakend)
#   panel   panel name (names(panels); a lone panel object is named by its $name)
#   match   "coord" | "symbol" — which matching mode produced the row
# A gene with several loci can appear more than once for one row/side, at different gaps —
# exactly as the labels show it. `panels` is one load_gene_panel() object or a named list.
sv_panel_hit_table = function(sv_table, panels,
                              bnd_window = SV_PANEL_WINDOW_BND,
                              other_window = SV_PANEL_WINDOW_OTHER) {
  empty = data.table(row = integer(), id = character(), gene = character(),
                     side = character(), gap = numeric(), direct = logical(),
                     panel = character(), match = character())
  n = if (is.null(sv_table)) 0L else nrow(sv_table)
  ps = .as_panel_list(panels)
  if (n == 0 || length(ps) == 0) return(empty)

  parts = lapply(names(ps), function(nm) {
    h = .sv_panel_hits_one(sv_table, ps[[nm]], bnd_window, other_window)
    if (nrow(h) > 0) h[, panel := nm]
    h
  })
  ht = rbindlist(c(list(empty[, !"id"]), parts), use.names = TRUE, fill = TRUE)
  if (nrow(ht) == 0) return(empty)
  ids = if ("id" %in% names(sv_table)) as.character(sv_table$id) else rep(NA_character_, n)
  ht[, id := ids[row]]
  ht[, direct := gap == 0]
  ht[, names(empty), with = FALSE]
}

# Labels from the long form, in its order; consecutive symbol-mode genes on one side of one
# row share a label ("A,B (A, direct)"), as the cell they were read from listed them.
.sv_hit_labels = function(ht, tagged) {
  ht = copy(ht)
  ht[, tag := if (tagged) paste0(" [", panel, "]") else ""]
  # Symbol-mode rows of one (row, side, panel) run collapse; coordinate rows stay one each
  ht[, grp := rleid(row, side, panel, fifelse(match == "symbol", 0L, seq_len(.N)))]
  ht[, .(row = row[1],
         label = if (match[1] == "symbol")
                   paste0(paste(gene, collapse = ","), " (", side[1], ", direct)", tag[1])
                 else
                   paste0(gene, " (", side, ", ",
                          fifelse(gap == 0, "direct", fmt_bp(gap)), ")", tag)),
     by = grp][, .(row, label)]
}

# One panel's hits, long form (row, gene, side, gap, match), in label order
.sv_panel_hits_one = function(sv_table, panel,
                              bnd_window = SV_PANEL_WINDOW_BND,
                              other_window = SV_PANEL_WINDOW_OTHER) {
  n = nrow(sv_table)
  none = data.table(row = integer(), gene = character(), side = character(),
                    gap = numeric(), match = character())
  if (is.null(panel)) return(none)

  if (!isTRUE(panel$has_coords)) {
    # The SV-scoped subset, so an `applies_to = snv` gene cannot match a breakend. Falls
    # back to `genes` for a panel object built before the scope column (test fixtures).
    symbols = toupper(unlist(if (!is.null(panel$genes_sv)) panel$genes_sv else panel$genes))
    hit_side = function(col, side) {
      if (!col %in% names(sv_table)) return(none)
      g = lapply(sv_table[[col]], function(cell) {
        if (is.na(cell) || !nzchar(cell)) return(character(0))
        g = trimws(unlist(strsplit(toupper(as.character(cell)), "[;,]+")))
        unique(g[g %in% symbols])
      })
      data.table(row = rep.int(seq_len(n), lengths(g)),
                 gene = as.character(unlist(g, use.names = FALSE)),
                 side = side, gap = 0, match = "symbol")
    }
    # A VEP symbol sits on the breakend itself, so a symbol hit is always direct
    h = rbindlist(list(hit_side("gene_a", "A"), hit_side("gene_b", "B")))
    if (nrow(h) == 0) return(none)
    # Row-major, side A before B, genes in cell order (order() is stable)
    return(h[order(row, side)])
  }

  iv = panel_intervals(panel)
  if (is.null(iv) || nrow(iv) == 0) return(none)

  is_bnd = sv_table$svtype %in% BND_SVTYPES
  # `start`/`end` are window-padded for foverlaps(); `qlo`/`qhi` carry the unpadded locus for the distance
  q = rbindlist(list(
    # Each BND breakend gets its own window (the sides may be on different contigs)
    data.table(row = which(is_bnd), side = "A",
               chrom = sv_table$chrom_a[is_bnd],
               start = sv_table$pos_a[is_bnd] - bnd_window,
               end   = sv_table$pos_a[is_bnd] + bnd_window,
               qlo   = sv_table$pos_a[is_bnd],
               qhi   = sv_table$pos_a[is_bnd]),
    data.table(row = which(is_bnd), side = "B",
               chrom = sv_table$chrom_b[is_bnd],
               start = sv_table$pos_b[is_bnd] - bnd_window,
               end   = sv_table$pos_b[is_bnd] + bnd_window,
               qlo   = sv_table$pos_b[is_bnd],
               qhi   = sv_table$pos_b[is_bnd]),
    # Other types are contiguous: one window around the whole span.
    data.table(row = which(!is_bnd), side = "span",
               chrom = sv_table$chrom_a[!is_bnd],
               start = pmin(sv_table$pos_a[!is_bnd], sv_table$pos_b[!is_bnd]) - other_window,
               end   = pmax(sv_table$pos_a[!is_bnd], sv_table$pos_b[!is_bnd]) + other_window,
               qlo   = pmin(sv_table$pos_a[!is_bnd], sv_table$pos_b[!is_bnd]),
               qhi   = pmax(sv_table$pos_a[!is_bnd], sv_table$pos_b[!is_bnd]))
  ))
  q = q[!is.na(chrom) & !is.na(start) & !is.na(end)]
  if (nrow(q) == 0) return(none)
  q[, start := pmax(as.numeric(start), 0)]
  q[, end   := as.numeric(end)]
  q[, `:=`(qlo = as.numeric(qlo), qhi = as.numeric(qhi))]
  iv = copy(iv)
  iv[, `:=`(start = as.numeric(start), end = as.numeric(end))]
  setkey(iv, chrom, start, end)

  ov = data.table::foverlaps(q, iv, by.x = c("chrom", "start", "end"),
                             type = "any", nomatch = NULL)
  if (nrow(ov) == 0) return(none)

  # `start`/`end` are the gene's interval here (foverlaps() renamed the query to i.start/i.end)
  ov[, gap := pmax(0, pmax(start - qhi, qlo - end))]
  # Rows in foverlaps() order, grouped by row (stable): within a row, side A's hits, then B's
  ov = ov[order(row)]
  # The same gene/side/gap twice (one symbol, two overlapping loci) is one label
  unique(ov[, .(row, gene, side, gap = as.numeric(gap), match = "coord")])
}

# Parse the gene-annotated Severus TSV (filtered_SV2/SV_filtered_with_gene_annotations.tsv)
parse_severus_gene_tsv = function(tsv_file) {
  if (is.null(tsv_file) || !file.exists(tsv_file)) return(NULL)
  dt = fread(tsv_file, sep = "\t", header = TRUE, fill = TRUE)
  setnames(dt, toupper(names(dt)))

  if ("START_CHROM" %in% names(dt)) dt[, START_CHROM := ensure_chr_prefix(START_CHROM)]
  if ("END_CHROM"   %in% names(dt)) dt[, END_CHROM   := ensure_chr_prefix(END_CHROM)]

  # Gene column: prefer NHL hits
  gene_col = if ("NHL_GENE_HITS"    %in% names(dt)) "NHL_GENE_HITS"
             else if ("COSMIC_GENE_HITS" %in% names(dt)) "COSMIC_GENE_HITS"
             else NULL
  dt[, gene_hits := if (!is.null(gene_col)) get(gene_col) else NA_character_]
  dt
}

# SV display table from the gene-annotated Severus TSV, mapped onto the same column contract as the VEP path
build_sv_table = function(sv_tsv) {
  if (is.null(sv_tsv) || nrow(sv_tsv) == 0) return(data.table())

  col = function(nm, default = NA) if (nm %in% names(sv_tsv)) sv_tsv[[nm]] else default
  chrom_a = ensure_chr_prefix(as.character(col("START_CHROM", NA_character_)))
  chrom_b = ensure_chr_prefix(as.character(col("END_CHROM",   NA_character_)))
  svtype  = as.character(col("SVTYPE", NA_character_))

  data.table(
    id      = as.character(col("ID", NA_character_)),
    svclass = fcase(
      is.na(chrom_a) | is.na(chrom_b), svtype,
      svtype %in% BND_SVTYPES & chrom_a != chrom_b, "translocation",
      svtype %in% BND_SVTYPES, "intra-chr breakend",
      default = svtype
    ),
    svtype  = svtype,
    chrom_a = chrom_a,
    pos_a   = suppressWarnings(as.integer(col("START_POS", NA_integer_))),
    # The TSV's gene hits are span-level, with no per-breakend split — hence side A only.
    gene_a  = as.character(if ("NHL_GENE_HITS" %in% names(sv_tsv)) sv_tsv$NHL_GENE_HITS
                           else col("COSMIC_GENE_HITS", NA_character_)),
    chrom_b = chrom_b,
    pos_b   = suppressWarnings(as.integer(col("END_POS", NA_integer_))),
    gene_b  = NA_character_,
    sv_len  = suppressWarnings(as.integer(col("SV_LEN", NA_integer_))),
    vaf     = suppressWarnings(as.numeric(col("VAF", NA_real_))),
    consequence = as.character(col("DETAILED_TYPE", NA_character_)),
    impact  = NA_character_,
    support = NA_integer_,
    svtype_raw = svtype,
    filter  = NA_character_,
    wl_rescue = NA_character_,
    orient_a = NA_integer_,
    orient_b = NA_integer_
  )
}
