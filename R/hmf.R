# HMF (Hartwig Medical Foundation / WiGiTS) driver, fusion and hotspot annotation: an
# optional section, on only when bin/render_report.R is given --hmf-driver-panel (and
# --hmf-ensembl-data-dir). Reads the user's own copy of the oncoanalyser resource bundle
# (v2.x or v3.x), GRCh38 only; nothing from it is bundled here.
#
# What it reproduces, and what it does not (the README's "HMF driver annotation" section
# says the same in plain language — a constant changed here must be changed there too):
#   - Small variants: PURPLE's reportable rule, flag-based. A variant in a DriverGenePanel
#     gene is classified on the gene's HMF canonical transcript (MISSENSE incl. inframe
#     indels, NONSENSE_OR_FRAMESHIFT, SPLICE) and reported when the panel flag for that
#     class is true, or when it is a known hotspot and reportSomaticHotspot is true. No dNdS
#     driver likelihood, no biallelic status, no germline reporting.
#   - Copy number: PURPLE's AMP / PARTIAL_AMP / DEL on the gene's canonical exons, gated by
#     reportAmplification / reportDeletion, using the report's own CN segments and ploidy.
#   - Fusions: LINX's known-fusion classification (KNOWN_PAIR, EXON_DEL_DUP, PROMISCUOUS_5/3,
#     IG/ENHANCER, PROMISCUOUS_ENHANCER_TARGET) and exon ranges on canonical transcripts,
#     for single SVs only (no chaining), with a phase check on the fused exons.
#   - ActionableCodingPanel.38.bed is a region flag only: HMF uses it for SAGE's PANEL
#     calling tier, it carries no actionability.

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Constants -------------------------------------------------------------------------
HMF_DEL_MAX_CN          = 0.5     # DEL: gene minimum exonic CN below this
HMF_AMP_RATIO_DEFAULT   = 3       # AMP: gene minimum exonic CN above ratio x ploidy
HMF_UPSTREAM_DEFAULT    = 10000   # 3' breakend upstream of the gene: promiscuous partners
HMF_UPSTREAM_KNOWN      = 100000  # ... known pairs and high-impact promiscuous partners
HMF_UPSTREAM_ENHANCER   = 1000000 # 3' gene distance from an IG/enhancer-range breakend
HMF_UNPHASED_KNOWN_MIN  = 1000000 # an out-of-frame known pair is reported only further apart (LINX SHORT_UNPHASED_DISTANCE_KNOWN)
HMF_SPLICE_CONSEQUENCES = c("splice_donor_variant", "splice_acceptor_variant",
                            "splice_donor_5th_base_variant")

# VEP consequence -> HMF coding effect; the first class a variant's consequences hit wins,
# in HMF's own severity order
HMF_CODING_EFFECTS = list(
  NONSENSE_OR_FRAMESHIFT = c("frameshift_variant", "stop_gained", "stop_lost", "start_lost"),
  SPLICE                 = HMF_SPLICE_CONSEQUENCES,
  MISSENSE               = c("missense_variant", "inframe_insertion", "inframe_deletion",
                             "protein_altering_variant")
)

# DriverGenePanel columns read, by header name; any absent column reads as false (v2.x has
# 15 columns, v3.x 22, a targeted panel may drop one)
HMF_PANEL_FLAGS = c(report_missense = "reportMissense", report_nonsense = "reportNonsense",
                    report_splice = "reportSplice", report_deletion = "reportDeletion",
                    report_disruption = "reportDisruption",
                    report_amplification = "reportAmplification",
                    report_hotspot = "reportSomaticHotspot")

# Which panel flag gates each coding effect
HMF_EFFECT_FLAG = c(MISSENSE = "report_missense", NONSENSE_OR_FRAMESHIFT = "report_nonsense",
                    SPLICE = "report_splice")

# Fusion types after normalising v2's IG_* names to v3's ENHANCER_*
HMF_ENHANCER_TYPES = c("ENHANCER_KNOWN_PAIR", "ENHANCER_PROMISCUOUS")

# ---- Loaders ---------------------------------------------------------------------------

hmf_bool = function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")

# DriverGenePanel.38.tsv: one row per gene with role and report flags
hmf_read_panel = function(path) {
  dt = fread(path, sep = "\t", colClasses = "character", na.strings = c("", "NA"))
  if (!"gene" %in% names(dt)) stop("HMF driver panel ", path, " has no 'gene' column")
  out = data.table(gene = trimws(dt$gene))
  for (nm in names(HMF_PANEL_FLAGS)) {
    col = HMF_PANEL_FLAGS[[nm]]
    out[, (nm) := if (col %in% names(dt)) hmf_bool(dt[[col]]) else FALSE]
  }
  lt = if ("likelihoodType" %in% names(dt)) toupper(trimws(dt$likelihoodType)) else NA_character_
  out[, likelihood_type := lt]
  out[, role := fcase(likelihood_type == "ONCO", "oncogene", likelihood_type == "TSG", "tsg",
                      default = NA_character_)]
  ratio = if ("amplificationRatio" %in% names(dt)) suppressWarnings(as.numeric(dt$amplificationRatio))
          else NA_real_
  out[, amp_ratio := fifelse(is.na(ratio), HMF_AMP_RATIO_DEFAULT, ratio)]
  out = out[nzchar(gene) & !is.na(gene)]
  out = unique(out, by = "gene")
  setattr(out, "n_columns", ncol(dt))
  out
}

# known_fusion_data.38.csv. Overrides is space-separated KEY=a;b;c; IG_RANGE (v2) and
# ENHANCER_RANGE (v3) are the same thing, as are the IG_* and ENHANCER_* types
hmf_read_fusions = function(path) {
  dt = fread(path, sep = ",", colClasses = "character", na.strings = c("", "NA"))
  need = c("Type", "FiveGene", "ThreeGene")
  if (!all(need %in% names(dt)))
    stop("HMF known fusion file ", path, " lacks column(s) ",
         paste(setdiff(need, names(dt)), collapse = ", "))
  col = function(nm) if (nm %in% names(dt)) as.character(dt[[nm]]) else rep(NA_character_, nrow(dt))
  raw_type = trimws(dt$Type)
  range_of = function(x, i) suppressWarnings(as.integer(vapply(strsplit(ifelse(is.na(x), "", x), ";", fixed = TRUE),
                                                               function(p) if (length(p) >= i) p[i] else NA_character_,
                                                               character(1))))
  up = col("KnownExonUpRange"); down = col("KnownExonDownRange")
  out = data.table(
    type        = sub("^IG_", "ENHANCER_", raw_type),
    five_gene   = trimws(col("FiveGene")),
    three_gene  = trimws(col("ThreeGene")),
    transcript  = sub("\\..*$", "", trimws(col("KnownExonTranscript"))),
    up_min = range_of(up, 1), up_max = range_of(up, 2),
    down_min = range_of(down, 1), down_max = range_of(down, 2),
    high_impact = hmf_bool(col("HighImpactPromiscuous")),
    cancer_types = col("CancerTypes")
  )
  out[five_gene == "", five_gene := NA_character_]
  out[three_gene == "", three_gene := NA_character_]

  # Overrides: the region keys parse to orient/chrom/start/end
  ov = col("Overrides")
  ov_val = function(keys) {
    vapply(ifelse(is.na(ov), "", ov), function(s) {
      kv = strsplit(s, "\\s+")[[1]]
      for (k in keys) {
        hit = kv[startsWith(kv, paste0(k, "="))]
        if (length(hit)) return(sub("^[^=]+=", "", hit[1]))
      }
      NA_character_
    }, character(1), USE.NAMES = FALSE)
  }
  region = function(v) {
    p = strsplit(ifelse(is.na(v), "", v), ";", fixed = TRUE)
    get = function(i) vapply(p, function(x) if (length(x) >= i) x[i] else NA_character_, character(1))
    list(orient = suppressWarnings(as.integer(get(1))),
         chrom  = ifelse(is.na(get(2)), NA_character_, ensure_chr_prefix(get(2))),
         start  = suppressWarnings(as.integer(get(3))),
         end    = suppressWarnings(as.integer(get(4))))
  }
  enh = region(ov_val(c("ENHANCER_RANGE", "IG_RANGE")))
  tp  = region(ov_val("THREE_PRIME_RANGE"))
  out[, `:=`(enh_orient = enh$orient, enh_chrom = enh$chrom, enh_start = enh$start, enh_end = enh$end,
             tp_orient = tp$orient, tp_chrom = tp$chrom, tp_start = tp$start, tp_end = tp$end,
             down_distance = suppressWarnings(as.integer(ov_val("DOWN_GENE_DOWNSTREAM_DISTANCE"))))]
  setattr(out, "format", if (any(startsWith(raw_type, "IG_"))) "v2 (IG_*)" else
                         if (any(startsWith(raw_type, "ENHANCER_"))) "v3 (ENHANCER_*)" else "unknown")
  out
}

# HMF Ensembl data cache: gene spans and the canonical transcript's exons, for every gene
# (a promiscuous fusion's partner can be any gene) or only `genes`. The exon file is ~200 MB;
# only the columns used are read, and only canonical rows kept (~250k)
hmf_read_ensembl = function(dir, genes = NULL) {
  gene_file = file.path(dir, "ensembl_gene_data.csv")
  exon_file = file.path(dir, "ensembl_trans_exon_data.csv")
  for (f in c(gene_file, exon_file))
    if (!file.exists(f)) stop("HMF Ensembl data cache file not found: ", f)

  g = fread(gene_file, sep = ",", colClasses = "character",
            select = c("GeneId", "GeneName", "Chromosome", "Strand", "GeneStart", "GeneEnd"))
  if (!is.null(genes)) g = g[GeneName %in% genes]
  g = g[, .(gene_id = GeneId, gene = GeneName, chrom = ensure_chr_prefix(Chromosome),
            strand = as.integer(Strand), start = as.integer(GeneStart), end = as.integer(GeneEnd))]
  # A symbol on two contigs (PAR genes): keep the first, as HMF's own lookup does
  g = unique(g, by = "gene")

  ex = fread(exon_file, sep = ",", na.strings = c("NULL", ""),
             select = c("GeneId", "CanonicalTranscriptId", "TransId", "TransName", "BioType",
                        "TransStart", "TransEnd", "ExonRank", "ExonStart", "ExonEnd",
                        "ExonPhase", "ExonEndPhase", "CodingStart", "CodingEnd"))
  ex = ex[GeneId %in% g$gene_id & CanonicalTranscriptId == TransId]
  ex = ex[, .(gene_id = GeneId, transcript = TransName, biotype = BioType,
              trans_start = as.integer(TransStart), trans_end = as.integer(TransEnd),
              rank = as.integer(ExonRank), start = as.integer(ExonStart), end = as.integer(ExonEnd),
              phase = as.integer(ExonPhase), end_phase = as.integer(ExonEndPhase),
              coding_start = as.integer(CodingStart), coding_end = as.integer(CodingEnd))]
  setorder(ex, gene_id, rank)
  ex = merge(ex, g[, .(gene_id, gene, chrom, strand)], by = "gene_id")
  setorder(ex, gene, rank)

  # Canonical transcript id and span per gene
  tr = unique(ex[, .(gene_id, gene, chrom, strand, transcript, biotype, trans_start, trans_end,
                     coding_start, coding_end)], by = "gene_id")
  # Per-gene exon tables, for the per-breakend lookups in hmf_annotate_breakend()
  list(genes = g, exons = ex, transcripts = tr, exons_by_gene = split(ex, by = "gene_id"))
}

# KnownHotspots.somatic.38.vcf.gz (SERVE-derived; `input=` carries the protein change in 1-letter form, p.V600E): coordinates plus the protein change and
# sources from INFO
hmf_read_hotspots = function(path) {
  con = if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  lines = readLines(con, warn = FALSE)
  close(con)
  body = lines[!startsWith(lines, "#")]
  if (length(body) == 0) return(data.table(chrom = character(), pos = integer(), ref = character(),
                                           alt = character(), hs_gene = character(),
                                           hs_protein = character(), sources = character()))
  f = tstrsplit(body, "\t", fixed = TRUE, keep = c(1L, 2L, 4L, 5L, 8L))
  info = f[[5]]
  info_val = function(key) {
    m = regexpr(paste0("(?:^|;)", key, "=([^;]+)"), info, perl = TRUE)
    out = rep(NA_character_, length(info))
    out[m > 0] = sub("^;?[^=]+=", "", regmatches(info, m))
    out
  }
  inp = strsplit(ifelse(is.na(info_val("input")), "", info_val("input")), "|", fixed = TRUE)
  part = function(i) vapply(inp, function(x) if (length(x) >= i) x[i] else NA_character_, character(1))
  dt = data.table(chrom = ensure_chr_prefix(f[[1]]), pos = as.integer(f[[2]]),
                  ref = f[[3]], alt = f[[4]], hs_gene = part(1), hs_protein = part(3),
                  sources = info_val("sources"))
  # Multi-allelic records: one row per ALT
  dt = dt[, .(alt = unlist(strsplit(alt, ",", fixed = TRUE))),
          by = .(chrom, pos, ref, hs_gene, hs_protein, sources)]
  dt[hs_protein %in% c("-", "null", ""), hs_protein := NA_character_]
  dt[]
}

# ActionableCodingPanel.38.bed(.gz): 0-based BED, name GENE_exonRank
hmf_read_bed = function(path) {
  dt = fread(path, sep = "\t", header = FALSE, colClasses = "character")
  if (ncol(dt) < 3) stop("HMF actionable coding BED ", path, " has fewer than 3 columns")
  data.table(chrom = ensure_chr_prefix(dt[[1]]), start = as.integer(dt[[2]]) + 1L,
             end = as.integer(dt[[3]]),
             gene = if (ncol(dt) >= 4) sub("_[0-9]+$", "", dt[[4]]) else NA_character_)
}

# The build the HMF resources are on. The published bundle is GRCh38 and declares nothing, so
# undeclared reads as "hg38"; a copy lifted to CHM13 (the pipeline's HMF_LIFTOVER) declares it
# in a one-line REFERENCE file inside the Ensembl data cache directory, the coordinate anchor
# of the whole set (each resource is staged on its own, so the declaration has to travel with
# one of them). A lifted hotspot VCF also carries ##assembly=, which must agree.
# Returns list(reference, hotspots) with hotspots NA when undeclared
hmf_declared_reference = function(files) {
  ref = "hg38"
  f = if (!is.null(files$ensembl_data_dir)) file.path(files$ensembl_data_dir, "REFERENCE") else ""
  if (nzchar(f) && file.exists(f)) {
    x = normalise_reference_name(readLines(f, n = 1, warn = FALSE)[1])
    if (!is.na(x)) ref = x
  }
  hs = NA_character_
  if (!is.null(files$hotspots) && file.exists(files$hotspots)) {
    con = if (grepl("\\.gz$", files$hotspots)) gzfile(files$hotspots, "rt") else file(files$hotspots, "rt")
    hdr = character(0)
    repeat {
      l = readLines(con, n = 200, warn = FALSE)
      hdr = c(hdr, l[startsWith(l, "##")])
      if (length(l) < 200 || any(!startsWith(l, "##"))) break
    }
    close(con)
    a = grep("^##assembly=", hdr, value = TRUE)
    if (length(a)) hs = normalise_reference_name(sub("^##assembly=", "", a[1]))
  }
  list(reference = ref, hotspots = hs)
}

# NULL when the HMF resources are on `reference`, else why not ("mismatch": declared another
# build; "undeclared": GRCh38 resources on another reference)
hmf_check_reference = function(files, reference) {
  d = hmf_declared_reference(files)
  if (!is.na(d$hotspots) && d$hotspots != d$reference)
    return(list(kind = "mismatch", msg = paste0("the HMF hotspot VCF declares '", d$hotspots,
      "' but the Ensembl data cache '", d$reference, "'; use resources lifted together")))
  if (d$reference == reference) return(NULL)
  list(kind = if (d$reference == "hg38") "undeclared" else "mismatch",
       msg = paste0("the HMF resources are on '", d$reference, "' and the report on '", reference, "'"))
}

# Every resource named in `files` (list of paths, NULL = not given). Panel and Ensembl
# cache are required for anything at all; the others switch their part on
hmf_load_resources = function(files) {
  if (is.null(files$driver_panel) || is.null(files$ensembl_data_dir)) return(NULL)
  panel = hmf_read_panel(files$driver_panel)
  fusions = if (!is.null(files$known_fusions)) hmf_read_fusions(files$known_fusions) else NULL
  # Every gene when fusions are matched (any gene can partner a promiscuous one), else the panel's
  ens = hmf_read_ensembl(files$ensembl_data_dir, if (is.null(fusions)) panel$gene else NULL)
  n_columns = attr(panel, "n_columns")   # merge() drops attributes
  panel = merge(panel, ens$transcripts[, .(gene, gene_id, chrom, strand, transcript)],
                by = "gene", all.x = TRUE)
  setattr(panel, "n_columns", n_columns)
  list(
    files    = files,
    panel    = panel,
    ensembl  = ens,
    fusions  = fusions,
    hotspots = if (!is.null(files$hotspots)) hmf_read_hotspots(files$hotspots) else NULL,
    bed      = if (!is.null(files$actionable_bed)) hmf_read_bed(files$actionable_bed) else NULL,
    # Not in the cache, or (in a lifted cache) dropped because a transcript did not lift whole
    missing_genes = setdiff(panel$gene, ens$genes$gene)
  )
}

# ---- Small variants --------------------------------------------------------------------

# First HMF coding effect a comma-joined VEP consequence string hits; NA for none
hmf_coding_effect = function(consequence) {
  toks = strsplit(ifelse(is.na(consequence), "", consequence), "[,&]")
  vapply(toks, function(t) {
    for (eff in names(HMF_CODING_EFFECTS)) if (any(t %in% HMF_CODING_EFFECTS[[eff]])) return(eff)
    NA_character_
  }, character(1))
}

# Protein change without the ENSP prefix: "ENSP0001.2:p.Val600Glu" -> "p.Val600Glu"
hmf_protein = function(hgvsp) {
  x = sub("^[^:]*:", "", as.character(hgvsp))
  x[is.na(hgvsp) | !startsWith(x, "p.")] = NA_character_
  x
}

# VEP writes HGVSp in 3-letter codes (p.Val600Glu), SERVE's hotspot `input=` in 1-letter
# (p.V600E): both go to 1-letter before comparing. Ter -> *, as HGVS does; a 1-letter string
# has no capital followed by lower case, so it passes through unchanged
HMF_AA3 = c(Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C", Gln = "Q", Glu = "E",
            Gly = "G", His = "H", Ile = "I", Leu = "L", Lys = "K", Met = "M", Phe = "F",
            Pro = "P", Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V", Ter = "*",
            Sec = "U", Pyl = "O", Xaa = "X")
hmf_protein_1letter = function(p) {
  for (k in names(HMF_AA3)) p = gsub(k, HMF_AA3[[k]], p, fixed = TRUE)
  p
}

# Is each position inside any BED interval of the same contig?
hmf_in_regions = function(chrom, pos, bed) {
  out = rep(FALSE, length(chrom))
  if (is.null(bed) || nrow(bed) == 0 || length(chrom) == 0) return(out)
  q = data.table(i = seq_along(chrom), chrom = chrom, start = as.integer(pos), end = as.integer(pos))
  b = bed[, .(chrom, start, end)]
  setkey(b, chrom, start, end)
  hit = foverlaps(q[!is.na(start)], b, by.x = c("chrom", "start", "end"), nomatch = NULL)
  out[unique(hit$i)] = TRUE
  out
}

# Driver small variants: every VEP row on a panel gene's HMF canonical transcript that has a
# coding effect, plus every variant at a known hotspot. `vep_data` is parse_vep()'s output
# (one row per variant x transcript annotation); `variant_table` supplies VAF
hmf_small_variants = function(vep_data, variant_table, res) {
  empty = data.table(gene = character(), chrom = character(), pos = integer(), ref = character(),
                     alt = character(), hgvsp = character(), consequence = character(),
                     effect = character(), hotspot = character(), role = character(),
                     reported = logical(), reason = character(), in_panel_region = logical(),
                     vaf = numeric(), in_panel = logical())
  if (is.null(vep_data) || nrow(vep_data) == 0) return(empty)
  panel = res$panel
  space = if ("coord_space" %in% names(vep_data)) vep_data$coord_space[1] else "vep"
  v = vep_data[, intersect(c("chrom", "pos", "ref", "alt", "symbol", "gene_id", "consequence",
                             "hgvsp", "feature", "canonical"), names(vep_data)), with = FALSE]
  if (!"feature" %in% names(v)) v[, feature := NA_character_]
  if (!"canonical" %in% names(v)) v[, canonical := NA_character_]
  v[, key := variant_key(chrom, pos, ref, alt, space = space)]

  # Panel gene of each row: Ensembl gene id first (symbols drift between Ensembl releases),
  # the symbol as fallback
  gi = match(v$gene_id, panel$gene_id)
  gs = match(v$symbol, panel$gene)
  v[, pidx := fifelse(!is.na(gi), gi, gs)]

  # Coding rows: the HMF canonical transcript, or VEP's CANONICAL flag where no transcript id
  # was annotated. When VEP's transcript ids are in another namespace altogether, none can equal
  # an HMF one: the CHM13 rapid-release cache (ENST0522...), a --refseq run (NM_...). Every
  # coding row would then be dropped in silence, so VEP's own canonical transcript stands in
  coding = v[!is.na(pidx)]
  coding[, hmf_tx := panel$transcript[pidx]]
  coding[, tx := sub("\\..*$", "", feature)]
  ns_ok = nrow(coding) == 0 || all(is.na(coding$tx)) || any(coding$tx %in% panel$transcript)
  coding = if (ns_ok) coding[(!is.na(tx) & !is.na(hmf_tx) & tx == hmf_tx) |
                             (is.na(tx) & canonical %in% "YES")]
           else coding[canonical %in% "YES"]
  transcript_rule = if (ns_ok) "hmf_canonical" else "vep_canonical"
  coding[, effect := hmf_coding_effect(consequence)]
  coding = coding[!is.na(effect)]

  # Hotspots, keyed in the same normalised space as the VEP rows
  hs = res$hotspots
  hot_exact = character(0)
  if (!is.null(hs) && nrow(hs) > 0) {
    hs = copy(hs)[, key := variant_key(chrom, pos, ref, alt, space = "vcf")]
    hot_exact = hs$key
  }
  v[, hot := key %in% hot_exact]

  # Same protein change on the same gene as a hotspot whose alleles differ (an MNV hotspot
  # called as separate SNVs, or a differently normalised indel)
  hot_protein = if (!is.null(hs) && nrow(hs) > 0)
    unique(paste(hs$hs_gene[!is.na(hs$hs_protein)], hmf_protein_1letter(hs$hs_protein[!is.na(hs$hs_protein)])))
    else character(0)

  rows = rbindlist(list(
    coding[, .(key, chrom, pos, ref, alt, gene = panel$gene[pidx], pidx, hgvsp, consequence, effect)],
    # Hotspot hits on any row: a panel gene, a non-panel gene or none (TERT promoter)
    v[hot == TRUE, .(key, chrom, pos, ref, alt,
                     gene = fifelse(!is.na(pidx), panel$gene[pidx], symbol), pidx,
                     hgvsp, consequence, effect = NA_character_,
                     # rank the annotation shown: HMF canonical transcript, VEP canonical, any
                     tx_rank = fifelse(!is.na(pidx) & sub("\\..*$", "", feature) %in% panel$transcript[pidx], 1L,
                                       fifelse(canonical %in% "YES", 2L, 3L)))]
  ), use.names = TRUE, fill = TRUE)
  if (nrow(rows) == 0) return(setattr(empty, "transcript_rule", transcript_rule))
  rows[is.na(tx_rank), tx_rank := 0L]
  # One row per variant and gene: the coding row (it carries the effect) over a hotspot-only one
  setorder(rows, key, gene, effect, tx_rank, na.last = TRUE)
  rows = unique(rows, by = c("key", "gene"))

  rows[, hs_match := fcase(key %in% hot_exact, "exact",
                           paste(gene, hmf_protein_1letter(hmf_protein(hgvsp))) %in% hot_protein, "protein",
                           default = "")]
  rows[, in_panel := !is.na(pidx)]
  # A hotspot outside the panel is one row per variant, named by the hotspot file's own gene,
  # and none at all where a panel row already carries it (TERT promoter: TERT plus intergenic)
  if (!is.null(hs) && nrow(hs) > 0) {
    off = !rows$in_panel
    rows[off, gene := hs$hs_gene[match(key, hs$key)]]
    rows = rows[in_panel | !key %in% rows[in_panel == TRUE, key]]
    rows = rbind(rows[in_panel == TRUE], unique(rows[in_panel == FALSE], by = "key"))
  }
  flag = function(nm) { x = panel[[nm]][rows$pidx]; !is.na(x) & x }
  rows[, flag_effect := FALSE]
  for (eff in names(HMF_EFFECT_FLAG)) {
    idx = which(rows$effect %in% eff)
    if (length(idx)) set(rows, idx, "flag_effect", flag(HMF_EFFECT_FLAG[[eff]])[idx])
  }
  rows[, flag_hotspot := flag("report_hotspot")]
  rows[, hotspot := nzchar(hs_match)]
  rows[, reported := in_panel & (flag_effect | (hotspot & flag_hotspot))]
  rows[, reason := fcase(
    reported & hotspot & flag_hotspot, "hotspot",
    reported,                         paste0(tolower(sub("_OR_", " or ", effect)), " reportable"),
    !in_panel & hotspot,              "hotspot; gene not in driver panel",
    in_panel & hotspot,               "hotspot; reportSomaticHotspot is false",
    default = paste0("report", c(MISSENSE = "Missense", NONSENSE_OR_FRAMESHIFT = "Nonsense",
                                 SPLICE = "Splice")[effect], " is false"))]
  rows[reported & hs_match == "protein" & flag_hotspot, reason := "hotspot (same protein change)"]
  rows[, role := panel$role[pidx]]
  rows[, in_panel_region := hmf_in_regions(chrom, pos, res$bed)]

  # VAF from the display table (same normalised key, by gene)
  rows[, vaf := NA_real_]
  if (!is.null(variant_table) && nrow(variant_table) > 0 && "vaf" %in% names(variant_table)) {
    vt = variant_table[, .(key = variant_key(chrom, pos, ref, alt, space = space), vaf)]
    vt = unique(vt[!is.na(vaf)], by = "key")
    rows[, vaf := vt$vaf[match(key, vt$key)]]
  }
  out = rows[, .(gene, chrom, pos, ref, alt, hgvsp = hmf_protein(hgvsp), consequence, effect,
                 hotspot = hs_match, role, reported, reason, in_panel_region, vaf, in_panel)]
  setorder(out, -reported, -in_panel, gene, pos)
  # Which transcript the effects were read from: "hmf_canonical" or "vep_canonical"
  setattr(out, "transcript_rule", transcript_rule)
  out[]
}

# ---- Copy number -----------------------------------------------------------------------

# PURPLE's gene calls from one segment set: min/max total CN over the canonical exons of every
# panel gene flagged for amplification or deletion. Wakhan's zero-coverage segments over
# reference gaps were already left out by cna_parse() (cna_gap_segments()), so a gene inside an
# unmappable gap reads as uncovered; one off the gaps is kept, as the homozygous deletion it is
hmf_cn_calls = function(segments, ploidy, label, res, sex = "female") {
  panel = res$panel[(report_amplification | report_deletion) & !is.na(gene_id)]
  ex = res$ensembl$exons[gene %in% panel$gene]
  if (is.null(segments) || nrow(segments) == 0 || nrow(ex) == 0) return(NULL)
  seg = copy(segments)
  seg = seg[!is.na(total), .(chrom = ensure_chr_prefix(chrom), start = as.integer(start),
                            end = as.integer(end), total = as.numeric(total))]
  setkey(seg, chrom, start, end)
  e = ex[, .(gene, chrom, start, end)]
  hit = foverlaps(e, seg, by.x = c("chrom", "start", "end"), nomatch = NULL)
  hit[, ov := pmin(end, i.end) - pmax(start, i.start) + 1L]

  span = ex[, .(chrom = chrom[1], start = min(start), end = max(end),
                exonic_bp = sum(end - start + 1L)), by = gene]
  agg = if (nrow(hit) > 0) hit[, .(min_cn = min(total), max_cn = max(total), covered_bp = sum(ov)), by = gene]
        else data.table(gene = character(), min_cn = numeric(), max_cn = numeric(), covered_bp = integer())
  g = merge(span, agg, by = "gene", all.x = TRUE)
  g[, coverage := fifelse(is.na(covered_bp), 0, pmin(1, covered_bp / exonic_bp))]
  g = merge(g, panel[, .(gene, report_amplification, report_deletion, amp_ratio, role)], by = "gene")

  male = tolower(as.character(sex)) %in% c("male", "xy", "m")
  ploidy = suppressWarnings(as.numeric(ploidy))
  # chrX in a male: PURPLE halves the amplification threshold
  g[, threshold := amp_ratio * ploidy * fifelse(male & chrom == "chrX", 0.5, 1)]
  g[, call := fcase(
    report_deletion & !is.na(min_cn) & min_cn < HMF_DEL_MAX_CN,                         "DEL",
    report_amplification & !is.na(threshold) & !is.na(min_cn) & min_cn > threshold,  "AMP",
    report_amplification & !is.na(threshold) & !is.na(max_cn) & max_cn > threshold,  "PARTIAL_AMP",
    default = NA_character_)]
  g[, `:=`(source = label, ploidy = ploidy)]
  g[]
}

# Panel genes flagged by reportDeletion whose canonical exons a deletion called by an SV caller
# removes (Severus DEL records, SAVANA +- pairs typed DEL), whatever the copy number says
hmf_sv_deletions = function(sv_table, res) {
  if (is.null(sv_table) || nrow(sv_table) == 0 || !"svtype" %in% names(sv_table)) return(NULL)
  d = sv_table[svtype == "DEL" & !is.na(chrom_a) & !is.na(pos_a) & !is.na(pos_b) &
               (is.na(chrom_b) | chrom_b == chrom_a)]
  if (nrow(d) == 0) return(NULL)
  panel = res$panel[report_deletion == TRUE]
  ex = res$ensembl$exons[gene %in% panel$gene, .(gene, chrom, start, end, rank)]
  if (nrow(ex) == 0) return(NULL)
  q = d[, .(sv_id = id, caller = if ("caller" %in% names(d)) caller else NA_character_,
            chrom = chrom_a, start = as.integer(pmin(pos_a, pos_b)),
            end = as.integer(pmax(pos_a, pos_b)), vaf = vaf,
            filter = if ("filter" %in% names(d)) filter else NA_character_)]
  setkey(q, chrom, start, end)
  hit = foverlaps(ex, q, by.x = c("chrom", "start", "end"), nomatch = NULL)
  if (nrow(hit) == 0) return(NULL)
  n_ex = ex[, .(n_exons = .N, g_start = min(start), g_end = max(end)), by = gene]
  out = hit[, .(chrom = chrom[1], first = min(rank), last = max(rank), n_hit = .N,
                sv_start = start[1], sv_end = end[1], vaf = vaf[1], filter = filter[1]),
            by = .(gene, caller, sv_id)]
  out = merge(out, n_ex, by = "gene")
  out = merge(out, panel[, .(gene, role)], by = "gene")
  out[, `:=`(start = g_start, end = g_end, call = "DEL",
             source = paste0(fcase(caller %in% "severus", "Severus", caller %in% "savana", "SAVANA",
                                   is.na(caller), "SV", default = caller), " SV"),
             note = paste0(hmf_exon_label(first, last, n_hit, n_exons), " · ", fmt_bp(sv_end - sv_start + 1L),
                           " deletion", fifelse(!is.na(filter) & !filter %in% c("PASS", "."),
                                                paste0(" · FILTER ", filter), "")))]
  out[]
}
hmf_exon_label = function(first, last, n_hit, n_exons)
  fifelse(n_hit == n_exons, paste0("all ", n_exons, " exons"),
          fifelse(first == last, paste0("exon ", first, " of ", n_exons),
                  paste0("exons ", first, "–", last, " of ", n_exons)))

# Copy-number driver calls: ASCAT and Wakhan each called on their own segments and ploidy,
# plus deletions called by the SV callers. One row per gene x source with a call; `support`
# lists every source calling the same type for the gene. `cna` is SECTION_DATA$cna
hmf_copy_number = function(cna, sv_table, res, sex = "female") {
  empty = data.table(gene = character(), chrom = character(), start = integer(), end = integer(),
                     source = character(), call = character(), min_cn = numeric(),
                     max_cn = numeric(), threshold = numeric(), coverage = numeric(),
                     note = character(), role = character(), support = character())
  per = list()
  if (!is.null(cna)) {
    wk = if (length(cna$profiles) > 0) cna$profiles[[1]] else NULL
    wk_ploidy = if (!is.null(wk) && !is.na(wk$ploidy)) wk$ploidy else
                if (identical(cna$source, "wakhan")) cna$ploidy else NA_real_
    if (!is.null(wk)) per$wakhan = hmf_cn_calls(wk$segments, wk_ploidy, "Wakhan", res, sex)
    as_seg = if (!is.null(cna$ascat_segments)) cna$ascat_segments else
             if (identical(cna$source, "ascat")) cna$segments else NULL
    as_ploidy = if (!is.null(cna$ascat_ploidy) && !is.na(cna$ascat_ploidy)) cna$ascat_ploidy else
                if (identical(cna$source, "ascat")) cna$ploidy else NA_real_
    per$ascat = hmf_cn_calls(as_seg, as_ploidy, "ASCAT", res, sex)
  }
  per = Filter(Negate(is.null), per)
  rows = rbindlist(lapply(per, function(g) g[!is.na(call), .(gene, chrom, start, end, source, call,
    min_cn, max_cn, threshold, coverage,
    note = fifelse(coverage < 1, paste0(round(100 * coverage), "% of exonic bp covered"), ""),
    role)]), fill = TRUE)
  svd = hmf_sv_deletions(sv_table, res)
  if (!is.null(svd))
    rows = rbind(rows, svd[, .(gene, chrom, start, end, source, call, min_cn = NA_real_,
                               max_cn = NA_real_, threshold = NA_real_, coverage = NA_real_,
                               note, role)], fill = TRUE)
  out = if (nrow(rows) > 0) rows else empty
  if (nrow(out) > 0) {
    # PARTIAL_AMP and AMP agree on direction; DEL stands alone
    out[, dir := fifelse(call == "DEL", "DEL", "AMP")]
    out[, support := paste(unique(source), collapse = ", "), by = .(gene, dir)]
    out[, dir := NULL]
    src_rank = function(x) fifelse(x == "Wakhan", 1L, fifelse(x == "ASCAT", 2L, 3L))
    out = out[order(call != "DEL", match(chrom, paste0("chr", c(1:22, "X", "Y"))), start, gene, src_rank(source))]
  }
  # Genes no segment covers, per source (an ASCAT gap can hide a homozygous deletion)
  setattr(out, "uncovered", lapply(per, function(g) g[coverage == 0, gene]))
  setattr(out, "sources", vapply(per, function(g) paste0(g$source[1], " (ploidy ",
                                  trimws(formatC(g$ploidy[1], format = "fg", digits = 3)), ")"), character(1)))
  out[]
}

# ---- Fusions ---------------------------------------------------------------------------

# LINX's exonic base phase (hmf-common TranscriptUtils.calcExonicCodingPhase): the HMF phase of
# the base at `pos` in exon `e` (coding bases up to and including it), -1 outside the CDS. The
# first coding base of an exon is tick(ExonPhase), so a base's phase is (ExonPhase + n) mod 3,
# with ExonPhase -1 read as 0 where the CDS starts in this exon
hmf_exonic_base_phase = function(e, cs, ce, strand, pos) {
  if (is.na(cs) || pos < cs || pos > ce) return(-1L)
  plus = strand == 1L
  n = if (plus) pos - max(e$start, cs) + 1L else min(e$end, ce) - pos + 1L
  before = if (plus) cs < e$start else ce > e$end      # the CDS started in an earlier exon
  as.integer(((if (before) e$phase else 0L) + n) %% 3L)
}

# Every canonical transcript a breakend can belong to, as 5' or 3' partner, the way LINX places
# it (linx BreakendGenePrep / BreakendTransData):
#   - region: exon k, intron k (after exon k), or upstream of the gene (3' partner only)
#   - context: coding / 5utr / 3utr / noncoding, from the position against the CDS; a 3' partner
#     must be coding, have a splice acceptor (>= 2 exons) and not be broken in its 3'UTR
#   - fused exon and phase, for a partner that is not itself exonic-to-exonic: the last exon kept
#     for a 5' partner (k after intron k; k - 1 for exon k, which is spliced round), with its end
#     phase; the first exon kept for a 3' partner (k + 1 after intron k or exon k; exon 2 when
#     upstream, since exon 1 has no splice acceptor), with its start phase, -1 where that exon
#     starts at the start codon. A 5' breakend in its 3'UTR gets phase -2, which matches nothing
#   - base_phase: for an exonic breakend, the phase of the breakend base itself, used when both
#     breakends are exonic (exact base phasing)
#   - alt_phases: phases reachable by skipping exons (LINX setAlternativeTranscriptPhasings): end
#     phases of the 5' exons before the fused one, start phases of the 3' exons after it (-2 in
#     the 3'UTR)
hmf_annotate_breakend = function(bk_chrom, pos, orient, res, upstream) {
  tr = res$ensembl$transcripts
  ex = res$ensembl$exons
  if (is.na(bk_chrom) || is.na(pos) || is.na(orient)) return(NULL)
  cand = tr[chrom == bk_chrom]
  if (nrow(cand) == 0) return(NULL)
  # `upstream`: a window per gene name, or one unnamed window for every gene
  up = if (is.null(names(upstream))) rep(upstream[1], nrow(cand)) else unname(upstream[cand$gene])
  if (length(up) != nrow(cand)) up = rep(NA_real_, nrow(cand))
  up[is.na(up)] = HMF_UPSTREAM_DEFAULT
  lo = fifelse(cand$strand == 1L, cand$trans_start - up, cand$trans_start)
  hi = fifelse(cand$strand == 1L, cand$trans_end, cand$trans_end + up)
  cand = cand[pos >= lo & pos <= hi]
  if (nrow(cand) == 0) return(NULL)

  rbindlist(lapply(seq_len(nrow(cand)), function(i) {
    t = cand[i]
    e = if (!is.null(res$ensembl$exons_by_gene)) res$ensembl$exons_by_gene[[t$gene_id]] else ex[gene_id == t$gene_id]
    is_5p = t$strand == orient
    if (!is_5p && (is.na(t$coding_start) || nrow(e) < 2)) return(NULL)
    plus = t$strand == 1L
    cs = t$coding_start; ce = t$coding_end
    context = if (is.na(cs)) "noncoding" else if (pos >= cs && pos <= ce) "coding" else
              if (plus == (pos > ce)) "3utr" else "5utr"
    if (!is_5p && context == "3utr") return(NULL)
    starts_cds = function(x) nrow(x) == 1 && !is.na(cs) && (if (plus) x$start == cs else x$end == ce)
    start_phase = function(x) if (nrow(x) == 0) NA_integer_ else if (starts_cds(x)) -1L else x$phase
    inside = pos >= t$trans_start && pos <= t$trans_end
    in_exon = which(e$start <= pos & e$end >= pos)
    base_phase = NA_integer_
    if (!inside) {
      if (is_5p) return(NULL)   # a 5' partner keeps its own promoter: never upstream of it
      region = "upstream"; k = 0L; fused = 2L
      dist = if (plus) t$trans_start - pos else pos - t$trans_end
      e2 = e[rank == 2L]
      # coding when the CDS already starts in exon 1 (LINX isCodingTypeOverride)
      context = if (if (plus) e2$start > cs else e2$end < ce) "coding" else "5utr"
      phase = start_phase(e2)
    } else if (length(in_exon)) {
      region = "exonic"; k = e$rank[in_exon[1]]; dist = 0L
      ek = e[in_exon[1]]
      base_phase = hmf_exonic_base_phase(ek, cs, ce, t$strand, pos)
      # Partner not exonic: exon k is spliced round (LINX nextSpliceExonRank)
      fused = if (is_5p) k - 1L else k + 1L
      phase = if (is_5p) ek$phase else start_phase(e[rank == fused])
    } else {
      # Intron k: between exon k and k + 1 in transcript order
      before = if (plus) e$end < pos else e$start > pos
      if (!any(before)) return(NULL)
      region = "intronic"; k = max(e$rank[before]); dist = 0L
      fused = if (is_5p) k else k + 1L
      fe = e[rank == fused]
      phase = if (is_5p) (if (nrow(fe) == 0) NA_integer_ else fe$end_phase) else start_phase(fe)
    }
    if (is_5p && context == "3utr") phase = -2L
    alt = if (is_5p) e[rank < fused, end_phase] else {
      d = e[rank > fused]
      fifelse(d$phase == -1L & (if (plus) d$end > ce else d$start < cs), -2L, d$phase)
    }
    data.table(gene = t$gene, strand = t$strand, side = if (is_5p) "5p" else "3p",
               region = region, exon = k, fused_exon = fused, phase = as.integer(phase),
               base_phase = base_phase, context = context,
               alt_phases = list(unique(setdiff(alt, phase))),
               upstream_distance = as.integer(dist), n_exons = nrow(e), pos = pos, chrom = bk_chrom)
  }), fill = TRUE)
}

# LINX's frame and reportability for one 5'/3' transcript pair (linx FusionFinder
# transcriptPairHasFusion, FusionReportability). NULL when the pair is no fusion at all:
# 5'UTR/non-coding -> coding and coding -> 5'UTR fuse only for a known pair. Phased: exact base
# phase when both breakends are exonic (exon -> exon); otherwise fused-exon phases equal
# (-1 == -1 is a 5'UTR -> 5'UTR fusion using the 3' partner's own start codon), or equal after
# skipping exons. Reported: a known pair in frame, or out of frame only from a 5'UTR/non-coding
# 5' breakend more than HMF_UNPHASED_KNOWN_MIN apart (or on another chromosome); an exon
# deletion/duplication in frame, or exonic at both ends; a promiscuous fusion only in frame,
# and through skipped exons only for a high-impact or known-exon partner
hmf_frame_rule = function(p, type, known_exons, hi_imp) {
  utr5 = p$context_5 %in% c("5utr", "noncoding")
  if (type != "KNOWN_PAIR" && ((utr5 && p$context_3 == "coding") ||
                               (p$context_5 == "coding" && p$context_3 == "5utr")))
    return(NULL)
  both_ex = p$region_5 == "exonic" && p$region_3 == "exonic"
  skipped = FALSE
  if (both_ex) {
    phased = if (p$context_5 == "coding" && p$context_3 == "coding")
               (p$base_phase_5 + 1L) %% 3L == p$base_phase_3   # LINX exonToExonInPhase: tick(5') == 3'
             else utr5 && p$context_3 == "5utr"
    phased = isTRUE(phased)
  } else {
    phased = isTRUE(p$phase_5 == p$phase_3)
    if (!phased && (p$gene_5 != p$gene_3 || type == "EXON_DEL_DUP")) {
      a5 = unlist(p$alt_phases_5); a3 = unlist(p$alt_phases_3)
      phased = isTRUE(p$phase_3 %in% a5 || p$phase_5 %in% a3 || length(intersect(a5, a3)) > 0)
      skipped = phased
    }
  }
  frame = if (phased && skipped) "skipped exons" else
          if (phased && !both_ex && p$phase_5 == -1L) "inframe (5′UTR)" else
          if (phased) "inframe" else
          if (utr5) "5′UTR to coding" else "out-of-frame"
  short = p$chrom_5 == p$chrom_3 && abs(p$pos_5 - p$pos_3) <= HMF_UNPHASED_KNOWN_MIN
  reported = switch(type,
    KNOWN_PAIR   = phased || (utr5 && !short),
    EXON_DEL_DUP = phased || both_ex,
    phased && (!skipped || hi_imp || known_exons))
  list(frame = frame, reported = reported)
}

# Every 5' candidate against every 3' one, columns suffixed _5 / _3; NULL when either is empty
hmf_cross = function(five, three) {
  if (nrow(five) == 0 || nrow(three) == 0) return(NULL)
  i = rep(seq_len(nrow(five)), each = nrow(three))
  j = rep(seq_len(nrow(three)), times = nrow(five))
  x = copy(five)[i]; y = copy(three)[j]
  setnames(x, paste0(names(x), "_5")); setnames(y, paste0(names(y), "_3"))
  cbind(x, y)
}

# Is pos within [lo, hi] of a region on the same contig, and the orientation as configured
# (0 = not checked)?
hmf_in_range = function(chrom, pos, orient, r_chrom, r_start, r_end, r_orient) {
  !is.na(r_chrom) & !is.na(chrom) & chrom == r_chrom & pos >= r_start & pos <= r_end &
    (is.na(r_orient) | r_orient == 0L | (!is.na(orient) & orient == r_orient))
}

# Known-fusion calls over the report's SV table (both callers, one row per rearrangement)
hmf_fusions = function(sv_table, res) {
  empty = data.table(caller = character(), sv_id = character(), svtype = character(),
                     type = character(), five_gene = character(), three_gene = character(),
                     five_exon = character(), three_exon = character(), frame = character(),
                     likelihood = character(), reported = logical(), high_impact = logical(),
                     note = character(), locus = character())
  fu = res$fusions
  if (is.null(fu) || is.null(sv_table) || nrow(sv_table) == 0) return(empty)
  if (!all(c("orient_a", "orient_b") %in% names(sv_table))) return(empty)

  known  = fu[type == "KNOWN_PAIR"]
  p5     = fu[type == "PROMISCUOUS_5"]
  p3     = fu[type == "PROMISCUOUS_3"]
  deldup = fu[type == "EXON_DEL_DUP"]
  enh    = fu[type %in% HMF_ENHANCER_TYPES & !is.na(enh_chrom)]
  target = fu[type == "PROMISCUOUS_ENHANCER_TARGET" & !is.na(tp_chrom)]

  # Upstream window per gene: wider for known pairs and high-impact partners
  wide = unique(c(known$three_gene, fu[high_impact == TRUE, c(five_gene, three_gene)]))
  genes = unique(na.omit(c(fu$five_gene, fu$three_gene)))
  upstream = setNames(ifelse(genes %in% wide, HMF_UPSTREAM_KNOWN, HMF_UPSTREAM_DEFAULT), genes)

  sv = sv_table[!is.na(chrom_a) & !is.na(pos_a)]
  if (!"caller" %in% names(sv)) sv[, caller := NA_character_]
  out = list()
  for (i in seq_len(nrow(sv))) {
    r = sv[i]
    locus = paste0(r$chrom_a, ":", r$pos_a, if (!is.na(r$chrom_b)) paste0(" / ", r$chrom_b, ":", r$pos_b) else "")
    base = list(caller = r$caller, sv_id = r$id, svtype = r$svtype, locus = locus)

    # Enhancer-driven: one side in an IG/enhancer range, the other near the 3' gene
    for (side in c("a", "b")) {
      ch = r[[paste0("chrom_", side)]]; ps = r[[paste0("pos_", side)]]; or = r[[paste0("orient_", side)]]
      oth = if (side == "a") "b" else "a"
      och = r[[paste0("chrom_", oth)]]; ops = r[[paste0("pos_", oth)]]
      if (is.na(ch) || is.na(ps)) next
      for (j in which(hmf_in_range(ch, ps, or, enh$enh_chrom, enh$enh_start, enh$enh_end, enh$enh_orient))) {
        e = enh[j]
        # Both breakends in the locus: a rearrangement inside it (V(D)J recombination at IG)
        if (isTRUE(hmf_in_range(och, ops, NA_integer_, e$enh_chrom, e$enh_start, e$enh_end, 0L))) next
        partner = NA_character_; ok = TRUE
        if (!is.na(e$three_gene)) {
          g = res$ensembl$genes[gene == e$three_gene]
          ok = nrow(g) == 1 && !is.na(och) && och == g$chrom && {
            dd = if (is.na(e$down_distance)) 0L else e$down_distance
            lo = if (g$strand == 1L) g$start - HMF_UPSTREAM_ENHANCER else g$start - dd
            hi = if (g$strand == 1L) g$end + dd else g$end + HMF_UPSTREAM_ENHANCER
            ops >= lo && ops <= hi
          }
          partner = e$three_gene
        } else {
          # Promiscuous: the partner is a gene whose 5' end the other breakend lies in or upstream of
          pb = hmf_annotate_breakend(och, ops, r[[paste0("orient_", oth)]], res,
                                     HMF_UPSTREAM_ENHANCER)
          pb = if (!is.null(pb)) pb[side == "3p"] else NULL
          ok = !is.null(pb) && nrow(pb) > 0
          partner = if (ok) paste(unique(pb$gene), collapse = ",") else NA_character_
        }
        if (!ok) next
        out[[length(out) + 1L]] = c(base, list(type = e$type, five_gene = e$five_gene,
          three_gene = partner, five_exon = NA_character_, three_exon = NA_character_,
          frame = "n/a", likelihood = if (e$type == "ENHANCER_KNOWN_PAIR") "HIGH" else "LOW",
          reported = e$type == "ENHANCER_KNOWN_PAIR", high_impact = FALSE,
          note = paste0(e$five_gene, " enhancer region")))
      }
      for (j in which(hmf_in_range(ch, ps, or, target$tp_chrom, target$tp_start, target$tp_end, target$tp_orient))) {
        tg = target[j]
        out[[length(out) + 1L]] = c(base, list(type = tg$type, five_gene = r[[paste0("gene_", oth)]],
          three_gene = tg$three_gene, five_exon = NA_character_, three_exon = NA_character_,
          frame = "n/a", likelihood = "LOW", reported = TRUE, high_impact = FALSE,
          note = paste0("breakend in the ", tg$three_gene, " enhancer-target region")))
      }
    }

    # Gene fusions: a 5' partner on one side and a 3' partner on the other (or the same gene)
    if (is.na(r$chrom_b) || is.na(r$pos_b)) next
    a = hmf_annotate_breakend(r$chrom_a, r$pos_a, r$orient_a, res, upstream)
    b = hmf_annotate_breakend(r$chrom_b, r$pos_b, r$orient_b, res, upstream)
    if (is.null(a) || is.null(b) || nrow(a) == 0 || nrow(b) == 0) next
    pairs = rbind(hmf_cross(a[side == "5p"], b[side == "3p"]),
                  hmf_cross(b[side == "5p"], a[side == "3p"]))
    if (is.null(pairs) || nrow(pairs) == 0) next
    for (k in seq_len(nrow(pairs))) {
      p = pairs[k]
      g5 = p$gene_5; g3 = p$gene_3
      # Exon -> exon: both exons are kept (and the 3' one may not be the last unless still 5'UTR);
      # otherwise an exonic side is spliced round, which needs an exon to splice to
      if (p$region_5 == "exonic" && p$region_3 == "exonic") {
        p[, `:=`(fused_exon_5 = exon_5, fused_exon_3 = exon_3)]
        if (p$exon_3 == p$n_exons_3 && p$context_3 != "5utr") next
      } else if (p$fused_exon_5 < 1L || p$fused_exon_3 > p$n_exons_3) next
      # EXON_DEL_DUP needs both breakends in one gene; any other same-gene pair is no fusion
      if (g5 == g3) {
        d = deldup[five_gene == g5 | three_gene == g5]
        if (nrow(d) == 0 || p$region_5 == "upstream" || p$region_3 == "upstream") next
        ok = any(p$fused_exon_5 >= d$up_min & p$fused_exon_5 <= d$up_max &
                 p$fused_exon_3 >= d$down_min & p$fused_exon_3 <= d$down_max, na.rm = TRUE)
        if (!ok) next
        type = "EXON_DEL_DUP"; known_exons = TRUE; hi_imp = FALSE
      } else if (nrow(known[five_gene == g5 & three_gene == g3]) > 0) {
        type = "KNOWN_PAIR"; known_exons = FALSE; hi_imp = FALSE
      } else {
        r3 = p3[three_gene == g3]; r5 = p5[five_gene == g5]
        if (nrow(r3) == 0 && nrow(r5) == 0) next
        in3 = nrow(r3) > 0 && (is.na(r3$down_min[1]) ||
                               (p$fused_exon_3 >= r3$down_min[1] && p$fused_exon_3 <= r3$down_max[1]))
        in5 = nrow(r5) > 0 && (is.na(r5$up_min[1]) ||
                               (p$fused_exon_5 >= r5$up_min[1] && p$fused_exon_5 <= r5$up_max[1]))
        type = if (nrow(r3) > 0 && nrow(r5) > 0) "PROMISCUOUS_BOTH" else
               if (nrow(r5) > 0) "PROMISCUOUS_5" else "PROMISCUOUS_3"
        known_exons = (nrow(r3) > 0 && !is.na(r3$down_min[1]) && in3) ||
                      (nrow(r5) > 0 && !is.na(r5$up_min[1]) && in5)
        hi_imp = (nrow(r3) > 0 && isTRUE(r3$high_impact[1]) && in3) ||
                 (nrow(r5) > 0 && isTRUE(r5$high_impact[1]) && in5)
      }
      # 3' breakend upstream of its gene: within the window this fusion type allows
      max_up = if (type == "KNOWN_PAIR" || hi_imp) HMF_UPSTREAM_KNOWN else HMF_UPSTREAM_DEFAULT
      if (p$region_3 == "upstream" && p$upstream_distance_3 > max_up) next

      fr = hmf_frame_rule(p, type, known_exons, hi_imp)
      if (is.null(fr)) next
      frame = fr$frame; reported = fr$reported
      likelihood = if (type %in% c("KNOWN_PAIR", "EXON_DEL_DUP") || known_exons) "HIGH" else "LOW"
      exon_lbl = function(region, exon, fused) switch(region,
        upstream = paste0("upstream → exon ", fused),
        exonic = if (exon == fused) paste0("exon ", exon) else paste0("exon ", exon, " → exon ", fused),
        intronic = paste0("intron ", exon, " → exon ", fused))
      out[[length(out) + 1L]] = c(base, list(type = type, five_gene = g5, three_gene = g3,
        five_exon = exon_lbl(p$region_5, p$exon_5, p$fused_exon_5),
        three_exon = exon_lbl(p$region_3, p$exon_3, p$fused_exon_3),
        frame = frame, likelihood = likelihood, reported = reported, high_impact = hi_imp,
        note = if (p$region_3 == "upstream") paste0("3' breakend ", fmt_bp(p$upstream_distance_3), " upstream") else ""))
    }
  }
  if (length(out) == 0) return(empty)
  res_dt = rbindlist(out, use.names = TRUE, fill = TRUE)
  setorder(res_dt, -reported, likelihood, five_gene, three_gene, caller)
  res_dt = unique(res_dt, by = c("caller", "sv_id", "type", "five_gene", "three_gene"))
  res_dt[, names(empty), with = FALSE]
}

# ---- Entry point -----------------------------------------------------------------------

# The whole HMF section's data, or NULL when it is off. `hmf_files` is params$hmf: a named
# list of paths (driver_panel, ensembl_data_dir, known_fusions, hotspots, actionable_bed)
build_hmf = function(hmf_files, vep_data, variant_table, section_data, ctx) {
  files = lapply(hmf_files, function(x) { x = as.character(unlist(x)); if (length(x) && nzchar(x[1])) x[1] else NULL })
  files = Filter(Negate(is.null), files)
  if (is.null(files$driver_panel)) return(NULL)
  res = tryCatch(hmf_load_resources(files), error = function(e) {
    message("HMF annotation disabled: ", conditionMessage(e))
    list(error = conditionMessage(e))
  })
  if (!is.null(res$error)) return(list(error = res$error, files = files))

  sex = if (!is.null(ctx$sex)) as.character(unlist(ctx$sex))[1] else "female"
  # A part that errors is NULL and its message goes into `failed` (named by part), so _hmf.qmd
  # says the step failed instead of rendering its "nothing found" notice. message() alone never
  # reaches the HTML: the setup chunk runs with message = FALSE
  failed = list()
  safely = function(expr, part, what) tryCatch(expr, error = function(e) {
    message("HMF ", what, " failed: ", conditionMessage(e))
    failed[[part]] <<- conditionMessage(e)
    NULL
  })
  small   = safely(hmf_small_variants(vep_data, variant_table, res), "small", "small-variant matching")
  cn      = safely(hmf_copy_number(section_data$cna, section_data$sv$table, res, sex), "cn", "copy-number calls")
  fusions = if (!is.null(res$fusions)) safely(hmf_fusions(section_data$sv$table, res), "fusions", "fusion matching")
  list(
    files     = files,
    reference = hmf_declared_reference(files)$reference,
    small     = small,
    cn        = cn,
    fusions   = fusions,
    failed    = failed,
    cn_source = if (!is.null(section_data$cna)) section_data$cna$source else NA_character_,
    ploidy    = if (!is.null(section_data$cna)) section_data$cna$ploidy else NA_real_,
    n_genes   = nrow(res$panel),
    n_panel_columns = attr(res$panel, "n_columns"),
    fusion_format = if (!is.null(res$fusions)) attr(res$fusions, "format") else NA_character_,
    n_hotspots = if (!is.null(res$hotspots)) nrow(res$hotspots) else NA_integer_,
    n_bed     = if (!is.null(res$bed)) nrow(res$bed) else NA_integer_,
    missing_genes = res$missing_genes
  )
}

# The driver panel as a report gene list (Key Findings chip): gene, canonical gene span from
# the Ensembl cache, role from likelihoodType. Written as a TSV so load_gene_panel() — the one
# loader every list goes through — reads it. Returns the path written
hmf_gene_list_tsv = function(panel_path, ensembl_dir, out_path, reference = "hg38") {
  panel = hmf_read_panel(panel_path)
  g = fread(file.path(ensembl_dir, "ensembl_gene_data.csv"), sep = ",", colClasses = "character",
            select = c("GeneName", "Chromosome", "GeneStart", "GeneEnd"))
  g = unique(g[GeneName %in% panel$gene], by = "GeneName")
  dt = merge(panel[, .(gene, role)], g, by.x = "gene", by.y = "GeneName", all.x = TRUE)
  dt = dt[!is.na(Chromosome)]
  out = dt[, .(gene, chrom = ensure_chr_prefix(Chromosome), start = GeneStart, end = GeneEnd,
               role = fifelse(is.na(role), "", role))]
  con = file(out_path, "w")
  writeLines(c(paste0("# reference: ", reference),
               "# HMF DriverGenePanel genes; spans from the HMF Ensembl data cache"), con)
  close(con)
  fwrite(out, out_path, sep = "\t", append = TRUE, col.names = TRUE)
  out_path
}
