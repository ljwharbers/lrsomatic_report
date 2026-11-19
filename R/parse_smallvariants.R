suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# INFO/CALLER values that denote a somatic caller (the merged VEP VCF is mostly germline records); ClairS is "clairs" or "clairs-to"
SOMATIC_CALLERS = c("clairs", "clairs-to", "clairsto", "deepsomatic")

# Derive dbsnp/cosmic columns from VEP's Existing_variation list
derive_dbsnp_cosmic = function(dt) {
  dt[, dbsnp := sub("(rs[0-9]+).*", "\\1", existing)]
  dt[!grepl("^rs", dbsnp, perl = TRUE), dbsnp := NA_character_]

  dt[, cosmic := sub(".*(COS[VM][0-9]+).*", "\\1", existing)]
  dt[!grepl("^COS", cosmic, perl = TRUE), cosmic := NA_character_]
  dt
}

# Report columns fed by VEP plugins or optional cache data, and the CSQ/Extra fields that can supply each (first present wins). Built-in SIFT/PolyPhen carry class and score in one field, `deleterious(0.01)`; the PolyPhen_SIFT plugin (T2T, whose cache has neither) splits them. humVar, not humDiv: it is what VEP's built-in PolyPhen reports
VEP_ANNOTATION_FIELDS = list(
  sift           = c("SIFT", "SIFT_pred"),
  sift_score     = c("SIFT", "SIFT_score"),
  polyphen       = c("PolyPhen", "PolyPhen_humvar_pred"),
  polyphen_score = c("PolyPhen", "PolyPhen_humvar_score"),
  am_class       = "am_class",
  am_score       = "am_pathogenicity",
  clinvar        = c("ClinVar_CLNSIG", "CLIN_SIG"),
  clinvar_id     = "ClinVar",
  cadd           = "CADD_PHRED",
  revel          = "REVEL",
  eve_class      = "EVE_CLASS",
  eve_score      = "EVE_SCORE"
)
VEP_SCORE_COLUMNS = c("sift_score", "polyphen_score", "am_score", "cadd", "revel", "eve_score")

# Annotation sources as the footnote names them, with the report columns each feeds
VEP_ANNOTATION_GROUPS = list(
  SIFT          = c("sift", "sift_score"),
  PolyPhen      = c("polyphen", "polyphen_score"),
  AlphaMissense = c("am_class", "am_score"),
  ClinVar       = c("clinvar", "clinvar_id"),
  CADD          = "cadd",
  REVEL         = "revel",
  EVE           = c("eve_class", "eve_score"),
  dbSNP         = "dbsnp",
  COSMIC        = "cosmic"
)

# Prediction classes tinted in the table; low-confidence and "possibly" calls stay plain
PATHOGENIC_CLASSES = c("deleterious", "probably_damaging", "pathogenic", "likely_pathogenic",
                       "Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")

# Which source field supplies each annotation column given the fields a file declares; NA where none does
resolve_vep_fields = function(available) {
  vapply(VEP_ANNOTATION_FIELDS, function(cands) {
    hit = cands[cands %in% available]
    if (length(hit)) hit[1] else NA_character_
  }, character(1))
}

# Whether the VEP cache carried dbSNP/COSMIC IDs, from the header (`##VEP=... dbSNP="156"` in a VCF, `## dbSNP version 156` in text output); NA when the header says nothing either way. The text scan is anchored at line start: every VEP text header describes its own Extra keys, and `## CLIN_SIG : ClinVar clinical significance of the dbSNP variant` matched an unanchored "dbSNP", which kept a permanently empty dbsnp column on every T2T text report
vep_cache_sources = function(header_lines) {
  vep_line = grep("^##VEP=", header_lines, value = TRUE)
  is_text  = any(grepl("^## ENSEMBL VARIANT EFFECT PREDICTOR", header_lines))
  if (length(vep_line) == 0 && !is_text) return(c(dbsnp = NA, cosmic = NA))
  said = function(key) {
    if (length(vep_line)) grepl(paste0("\\b", key, "[= ]"), vep_line[1])
    else any(grepl(paste0("^## ", key, "[ =]"), header_lines))
  }
  c(dbsnp = said("dbSNP"), cosmic = said("COSMIC"))
}

# Attach the annotation columns a file declares. `get_field(name)` returns that CSQ/Extra field as a character vector. A column whose source is not declared is not created (the plugin did not run); one declared but empty stays all NA (it ran and found nothing). The resolution is kept as attr "vep_fields" for the footnote
add_vep_annotations = function(dt, get_field, available, cache = c(dbsnp = NA, cosmic = NA)) {
  resolved = resolve_vep_fields(available)
  for (col in names(resolved)) {
    if (is.na(resolved[[col]])) next
    dt[, (col) := as.character(get_field(resolved[[col]]))]
  }
  # Built-in `deleterious(0.01)`: one field feeds both the class and the score column
  for (base in c("sift", "polyphen")) {
    sc = paste0(base, "_score")
    if (!is.na(resolved[[base]]) && identical(resolved[[base]], resolved[[sc]])) {
      raw = dt[[base]]
      dt[, (sc)   := sub("^.*\\(([^()]*)\\)\\s*$", "\\1", raw)]
      dt[, (base) := sub("\\(.*$", "", raw)]
      dt[!grepl("(", raw, fixed = TRUE), (sc) := NA_character_]
    }
  }
  for (col in intersect(names(resolved), names(dt))) {
    dt[get(col) == "", (col) := NA_character_]
  }
  # Multi-valued custom fields arrive "&"-joined, like Consequence; the facet splits on "," and the ClinVar link render splits the id cell the same way
  for (col in intersect(c("clinvar", "clinvar_id"), names(dt))) {
    dt[, (col) := gsub("&", ",", get(col), fixed = TRUE)]
  }
  for (col in intersect(VEP_SCORE_COLUMNS, names(dt))) {
    dt[, (col) := suppressWarnings(as.numeric(get(col)))]
  }
  # dbSNP/COSMIC come from the cache, not a plugin: dropped only when the header says the cache has none
  for (col in c("dbsnp", "cosmic")) {
    if (isFALSE(cache[[col]])) dt[, (col) := NULL]
    resolved[col] = if (isFALSE(cache[[col]])) NA_character_ else "Existing_variation"
  }
  setattr(dt, "vep_fields", resolved)
  dt
}

# One row per annotation source for the footnote: the columns it feeds, the CSQ fields that fed them, and whether this run had it at all
vep_annotation_summary = function(vep_fields) {
  if (is.null(vep_fields)) return(NULL)
  rbindlist(lapply(names(VEP_ANNOTATION_GROUPS), function(g) {
    cols = VEP_ANNOTATION_GROUPS[[g]]
    src  = unique(unname(vep_fields[cols]))
    src  = src[!is.na(src)]
    data.table(source = g, columns = paste(cols, collapse = ", "),
               fields = paste(src, collapse = ", "), present = length(src) > 0)
  }))
}

# Dispatch on file contents: VEP text output (#Uploaded_variation header) vs VCF with CSQ (#CHROM header)
parse_vep = function(vep_file) {
  if (is.null(vep_file) || !file.exists(vep_file)) return(NULL)

  con = gzfile(vep_file, "rb")
  is_vcf = FALSE
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#Uploaded_variation")) break
    if (startsWith(line, "#CHROM")) { is_vcf = TRUE; break }
  }
  close(con)

  if (is_vcf) parse_vep_vcf(vep_file) else parse_vep_text(vep_file)
}

# Parse VEP default text output (tab-delimited, not a VCF); one row per consequence per variant
parse_vep_text = function(vep_file) {
  if (is.null(vep_file) || !file.exists(vep_file)) return(NULL)

  # Count meta-lines (start with ##) to find the column-header line
  con = gzfile(vep_file, "rb")
  skip_n = 0L
  hdr = character(0)
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#Uploaded_variation")) break
    hdr = c(hdr, line)
    skip_n = skip_n + 1L
  }
  close(con)

  dt = tryCatch(
    fread(vep_file, skip = skip_n, sep = "\t", header = TRUE,
          col.names = function(x) gsub("^#", "", x)),
    error = function(e) {
      message("Failed to parse VEP file: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  setnames(dt, old = "Uploaded_variation", new = "variant_id", skip_absent = TRUE)
  setnames(dt, old = "Gene",              new = "gene_id",     skip_absent = TRUE)
  setnames(dt, old = "Consequence",       new = "consequence", skip_absent = TRUE)

  # Coordinates and alleles both from variant_id where it has the canonical shape; Location is off by one for dash-form insertions, so it is only the fallback
  vid = "^.+_[0-9]+_[^_]+/[^_]+$"
  dt[, from_vid := grepl(vid, variant_id)]

  dt[, chrom := ifelse(from_vid, sub("_[0-9]+_[^_]+$", "", variant_id),
                                 sub(":.*", "", Location))]
  dt[, pos   := as.integer(ifelse(from_vid, sub(".*_([0-9]+)_[^_]+$", "\\1", variant_id),
                                            sub(".*:(\\d+).*", "\\1", Location)))]
  dt[, chrom := ensure_chr_prefix(chrom)]

  dt[, ref := sub(".*_([^/]+)/.*", "\\1", variant_id)]
  dt[, alt := sub(".*/", "", variant_id)]

  # Extra keys the file declares (`## KEY : description` header lines), plus any annotation key actually present in a row, for files whose header omits them
  cands    = unique(unlist(VEP_ANNOTATION_FIELDS, use.names = FALSE))
  declared = sub("^## (\\S+) : .*$", "\\1", grep("^## \\S+ : ", hdr, value = TRUE))

  # One alternation pass over Extra, only for keys the header did not already declare, and the matches are extracted only from the rows that hit: one grepl() per key was 8.8 s of every render at 171k rows, this is 0.8 s
  unseen   = setdiff(cands, declared)
  observed = character(0)
  if (length(unseen) > 0) {
    xtr = as.character(dt$Extra)
    xtr[is.na(xtr)] = ""
    pat = paste0("(?:^|;)(", paste(unseen, collapse = "|"), ")=")
    hit = grepl(pat, xtr, perl = TRUE)
    if (any(hit)) {
      hits = unlist(regmatches(xtr[hit], gregexpr(pat, xtr[hit], perl = TRUE)), use.names = FALSE)
      observed = intersect(unseen, unique(sub("^;?(.*)=$", "\\1", hits)))
    }
  }
  available = union(declared, observed)

  # Parse the VEP Extra key=value field in one pass
  ex = extract_extra_keys(dt$Extra, unique(c("SYMBOL", "IMPACT", "HGVSp",
                                             intersect(cands, available))))
  dt[, symbol := ex$SYMBOL]
  dt[, impact := ex$IMPACT]
  dt[, hgvsp  := ex$HGVSp]

  # Existing_variation is a column of the text format, never an Extra key; reading it from Extra left it - and so dbsnp/cosmic - empty on every text report
  dt[, existing := if ("Existing_variation" %in% names(dt)) as.character(Existing_variation)
                   else NA_character_]
  dt[existing %in% c("", "-"), existing := NA_character_]

  # dbSNP / COSMIC IDs, derived from Existing_variation
  dt = derive_dbsnp_cosmic(dt)
  dt = add_vep_annotations(dt, function(nm) ex[[nm]], available, vep_cache_sources(hdr))

  # No per-variant caller in the text format; kept for contract parity with parse_vep_vcf()
  dt[, caller := NA_character_]

  # Record identity for parity with parse_vep_vcf(): VEP's own Uploaded_variation name, never a caller ID
  dt[, id := variant_id]

  # chrom/pos/ref/alt are in VEP notation (indels shifted, possibly dash-form); see coord_space in build_variant_table()
  dt[, coord_space := "vep"]

  dt
}

# Parse a VCF with VEP CSQ annotation; same column contract as parse_vep_text()
parse_vep_vcf = function(vep_file) {
  if (is.null(vep_file) || !file.exists(vep_file)) return(NULL)

  # Skip header to #CHROM, capturing the CSQ field order and whether a CALLER tag exists
  con = gzfile(vep_file, "rb")
  skip_n = 0L
  csq_format = NULL
  has_caller_info = FALSE
  hdr = character(0)
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    hdr = c(hdr, line)
    if (startsWith(line, "##INFO=<ID=CSQ")) {
      m = regmatches(line, regexpr("Format: [^\"]+", line))
      if (length(m) > 0) csq_format = strsplit(sub("^Format: ", "", m), "|", fixed = TRUE)[[1]]
    }
    if (startsWith(line, "##INFO=<ID=CALLER")) has_caller_info = TRUE
    if (startsWith(line, "#CHROM")) break
    skip_n = skip_n + 1L
  }
  close(con)

  if (is.null(csq_format)) {
    message("Failed to parse VEP VCF: no CSQ Format found in header")
    return(NULL)
  }

  dt = tryCatch(
    fread(vep_file, skip = skip_n + 1L, sep = "\t", header = FALSE, select = 1:8,
          col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")),
    error = function(e) {
      message("Failed to parse VEP VCF: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  # Drop records the callers themselves rejected (RefCall / LowQual / GERMLINE).
  dt = dt[FILTER %in% c("PASS", ".")]
  if (nrow(dt) == 0) return(NULL)

  # Fixed-string splits: much faster than sub() over millions of rows
  dt[, CSQ := tstrsplit(INFO, "CSQ=", fixed = TRUE, keep = 2L)[[1]]]
  dt[, CSQ := tstrsplit(CSQ,  ";",    fixed = TRUE, keep = 1L)[[1]]]

  # Keep somatic-caller records only; skipped for older VCFs with no CALLER tag at all
  if (has_caller_info) {
    dt[, caller := tstrsplit(INFO, "CALLER=", fixed = TRUE, keep = 2L)[[1]]]
    dt[, caller := tstrsplit(caller, ";", fixed = TRUE, keep = 1L)[[1]]]
    dt = dt[caller %in% SOMATIC_CALLERS]
    if (nrow(dt) == 0) return(NULL)
  } else {
    dt[, caller := NA_character_]
  }

  # Collapse multi-caller duplicates to one row per variant before CSQ expansion; ID rides along for the SV join
  dt = dt[, .(CSQ = CSQ[1], id = ID[1],
              caller = paste(sort(unique(caller)), collapse = ",")),
          by = .(CHROM, POS, REF, ALT)]
  dt[id == ".", id := NA_character_]

  # One row per gene/transcript annotation (comma-separated CSQ entries)
  dt_long = dt[, .(csq_entry = unlist(strsplit(CSQ, ",", fixed = TRUE))),
               by = .(CHROM, POS, REF, ALT, id, caller)]
  if (nrow(dt_long) == 0) return(NULL)

  # Split each entry on "|", keeping the fields used below plus whichever annotation fields this header declares. `need` is not in CSQ order - a real GRCh38 header puts ClinVar_CLNSIG after CADD_PHRED - and names(parts) relies on tstrsplit() returning `keep` in the order asked for, not sorted
  resolved = resolve_vep_fields(csq_format)
  need = unique(c("Consequence", "IMPACT", "SYMBOL", "Gene", "HGVSp", "Existing_variation",
                  resolved[!is.na(resolved)]))
  keep_idx = match(need, csq_format)
  ok = !is.na(keep_idx)
  parts = tstrsplit(dt_long$csq_entry, "|", fixed = TRUE, fill = NA_character_,
                    keep = keep_idx[ok])
  names(parts) = need[ok]

  get_field = function(name) {
    if (is.null(parts[[name]])) rep(NA_character_, nrow(dt_long)) else parts[[name]]
  }

  dt_long[, consequence := gsub("&", ",", get_field("Consequence"))]
  dt_long[, impact      := get_field("IMPACT")]
  dt_long[, symbol      := get_field("SYMBOL")]
  dt_long[, gene_id     := get_field("Gene")]
  dt_long[, hgvsp       := get_field("HGVSp")]
  dt_long[, existing    := get_field("Existing_variation")]

  # Blank CSQ fields are "", not NA; normalise for parity with parse_vep_text()
  for (col in c("symbol", "impact", "hgvsp", "existing", "gene_id", "caller")) {
    dt_long[get(col) == "", (col) := NA_character_]
  }

  dt_long[, chrom := ensure_chr_prefix(CHROM)]
  dt_long[, pos   := as.integer(POS)]
  dt_long[, ref   := REF]
  dt_long[, alt   := ALT]

  # These are the VCF's own POS/REF/ALT (not VEP-shifted); see build_variant_table()
  dt_long[, coord_space := "vcf"]

  # dbSNP / COSMIC IDs, derived from Existing_variation
  dt_long = derive_dbsnp_cosmic(dt_long)
  dt_long = add_vep_annotations(dt_long, get_field, csq_format, vep_cache_sources(hdr))

  keep = c("chrom", "pos", "ref", "alt", "id", "symbol", "gene_id", "consequence", "impact",
           "hgvsp", "existing", "dbsnp", "cosmic", names(VEP_ANNOTATION_FIELDS),
           "caller", "coord_space")
  out = dt_long[, intersect(keep, names(dt_long)), with = FALSE]
  setattr(out, "vep_fields", attr(dt_long, "vep_fields"))
  out
}

# Pick the sample column: one is unambiguous; else the report's sample, then the first non-normal; the first column is a flagged guess
pick_sample_column = function(samples, sample_id = NULL, vcf_file = "") {
  if (length(samples) == 1L) return(1L)

  if (!is.null(sample_id) && nzchar(sample_id)) {
    i = match(sample_id, samples)
    if (!is.na(i)) return(i)
  }
  not_normal = which(!grepl("^(normal|blood|germline|ref)", samples, ignore.case = TRUE))
  if (length(not_normal) > 0) {
    if (length(not_normal) > 1L || is.null(sample_id))
      warning("VCF '", vcf_file, "' has samples [", paste(samples, collapse = ", "),
              "]; reading '", samples[not_normal[1]], "'.")
    return(not_normal[1])
  }
  warning("VCF '", vcf_file, "' has samples [", paste(samples, collapse = ", "),
          "] and none looks like a tumour; reading '", samples[1], "'.")
  1L
}

# Allele fraction from whichever FORMAT tag is present: LRSomatic renames AF/VAF per caller (STANDARDIZE_AF); AD is the last resort
allele_fraction = function(field, fields) {
  for (tag in c("AF", "VAF")) {
    if (tag %in% fields) {
      # Multi-allelic: the first value pairs with the first ALT, which is all this table joins on
      v = suppressWarnings(as.numeric(sub(",.*$", "", field(tag))))
      if (any(!is.na(v))) return(v)
    }
  }
  if ("AD" %in% fields) {
    ad = strsplit(field("AD"), ",", fixed = TRUE)
    return(vapply(ad, function(x) {
      n = suppressWarnings(as.numeric(x))
      if (length(n) < 2L || anyNA(n[1:2]) || sum(n[1:2]) == 0) return(NA_real_)
      n[2] / sum(n[1:2])
    }, numeric(1)))
  }
  rep(NA_real_, length(field("GT")))
}

# Parse raw caller VCF(s) into chrom, pos, ref, alt, vaf, dp, gt, ps, caller; several paths are stacked (ClairS splits snvs/indels). `sample_id` selects the column by name so a matched VCF never reports the normal's values
parse_caller_vcf = function(vcf_file, caller_name = "unknown", sample_id = NULL) {
  if (is.null(vcf_file)) return(NULL)
  if (length(vcf_file) > 1) {
    parts = lapply(vcf_file, parse_caller_vcf, caller_name = caller_name,
                   sample_id = sample_id)
    parts = parts[!vapply(parts, is.null, logical(1))]
    return(if (length(parts) > 0) rbindlist(parts) else NULL)
  }
  if (!file.exists(vcf_file)) return(NULL)

  # Count header lines, keeping the #CHROM line itself — it carries the sample names.
  con = gzfile(vcf_file, "rb")
  skip_n = 0L
  chrom_line = NULL
  has_caller_tag = FALSE; has_dv = FALSE; sources = character(0)
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#CHROM")) { chrom_line = line; break }
    if (startsWith(line, "##INFO=<ID=CALLER,")) has_caller_tag = TRUE
    if (startsWith(line, "##DeepVariant_version")) has_dv = TRUE
    if (startsWith(line, "##source=")) sources = c(sources, line)
    skip_n = skip_n + 1L
  }
  close(con)
  if (is.null(chrom_line)) return(NULL)

  header  = strsplit(chrom_line, "\t", fixed = TRUE)[[1]]
  samples = if (length(header) > 9L) header[10:length(header)] else character(0)
  # Sites-only VCF: no genotypes, return NULL
  if (length(samples) == 0) return(NULL)

  sample_col = pick_sample_column(samples, sample_id, vcf_file)

  # 1:9 are the fixed VCF columns; the chosen sample sits at 9 + its index.
  col_names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                "SAMPLE1")
  dt = tryCatch(
    fread(vcf_file, skip = skip_n + 1L, sep = "\t", header = FALSE,
          select = c(1:9, 9L + sample_col), col.names = col_names),
    error = function(e) {
      warning("Could not read VCF '", vcf_file, "': ", conditionMessage(e))
      NULL
    })
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  dt[, CHROM := ensure_chr_prefix(CHROM)]

  # Extract AF/DP/GT/PS from FORMAT + SAMPLE, format-group by format-group
  fmt_groups = unique(dt$FORMAT)
  vaf_list = rep(NA_real_,      nrow(dt))
  dp_list  = rep(NA_integer_,   nrow(dt))
  gt_list  = rep(NA_character_, nrow(dt))
  ps_list  = rep(NA_character_, nrow(dt))

  for (fmt in fmt_groups) {
    idx_rows = which(dt$FORMAT == fmt)
    fields  = strsplit(fmt, ":", fixed = TRUE)[[1]]
    split_s = strsplit(dt$SAMPLE1[idx_rows], ":", fixed = TRUE)

    # One FORMAT field, by name, across this group's rows
    field = function(name) {
      i = match(name, fields)
      if (is.na(i)) return(rep(NA_character_, length(split_s)))
      vapply(split_s, function(x) if (length(x) >= i) x[i] else NA_character_,
             character(1))
    }

    vaf_list[idx_rows] = allele_fraction(field, fields)
    dp_list[idx_rows]  = suppressWarnings(as.integer(field("DP")))
    gt_list[idx_rows]  = field("GT")
    ps_list[idx_rows]  = field("PS")
  }

  # Blank the unphased placeholders ("." PS, "/" GT) so cells render empty
  ps_list[!is.na(ps_list) & ps_list == "."] = NA_character_
  gt_list[!is.na(gt_list) & gt_list %in% c(".", "./.")] = NA_character_

  # Which caller wrote each record: INFO/CALLER when the header declares it, else a
  # DeepSomatic-only VCF (DeepVariant header, no other caller's ##source)
  vcf_caller = if (has_caller_tag) {
    m = regmatches(dt$INFO, regexec("(^|;)CALLER=([^;]+)", dt$INFO))
    tolower(vapply(m, function(x) if (length(x) == 3) x[3] else NA_character_, character(1)))
  } else if (has_dv && !any(grepl("clair", sources, ignore.case = TRUE))) {
    rep("deepsomatic", nrow(dt))
  } else rep(NA_character_, nrow(dt))
  # DeepSomatic writes GT 1/1 on every somatic call by convention: not a genotype, so not shown
  gt_list[vcf_caller %in% "deepsomatic"] = NA_character_

  data.table(chrom = dt$CHROM, pos = dt$POS, ref = dt$REF, alt = dt$ALT,
             vaf = vaf_list, dp = dp_list, gt = gt_list, ps = ps_list,
             caller = caller_name, vcf_caller = vcf_caller)
}

# Header-only provenance for the VCF(s) supplying VAF/DP/GT/PS (after a consensus merge the FORMAT fields need not match the `callers` column); NULL when there is nothing to describe
vaf_provenance = function(vcf_files, sample_dir = NULL, sample_id = NULL) {
  if (is.null(vcf_files) || length(vcf_files) == 0) return(NULL)
  vcf_files = vcf_files[!is.na(vcf_files) & nzchar(vcf_files)]
  if (length(vcf_files) == 0) return(NULL)

  # Paths read better relative to the sample directory the caller passed in.
  rel = function(p) {
    if (is.null(sample_dir) || !nzchar(sample_dir)) return(p)
    root = sub("/+$", "", normalizePath(sample_dir, mustWork = FALSE))
    full = normalizePath(p, mustWork = FALSE)
    if (startsWith(full, paste0(root, "/"))) substring(full, nchar(root) + 2L) else p
  }

  # "##source=..." is optional colour; many VCFs carry none
  read_sources = function(p) {
    if (!file.exists(p)) return(character(0))
    con = gzfile(p, "rb")
    on.exit(close(con), add = TRUE)
    out = character(0)
    repeat {
      line = readLines(con, n = 1, warn = FALSE)
      if (length(line) == 0) break
      if (!startsWith(line, "##")) break   # #CHROM or a malformed header ends the scan
      if (startsWith(line, "##source=")) out = c(out, sub("^##source=", "", line))
    }
    out
  }

  # Which sample column parse_caller_vcf() read; named in the footnote because a wrong pick is otherwise invisible
  read_sample = function(p) {
    if (!file.exists(p)) return(NA_character_)
    con = gzfile(p, "rb")
    on.exit(close(con), add = TRUE)
    repeat {
      line = readLines(con, n = 1, warn = FALSE)
      if (length(line) == 0) return(NA_character_)
      if (startsWith(line, "#CHROM")) {
        header = strsplit(line, "\t", fixed = TRUE)[[1]]
        if (length(header) <= 9L) return(NA_character_)
        samples = header[10:length(header)]
        # Only worth reporting when there was a choice to get wrong.
        if (length(samples) == 1L) return(NA_character_)
        return(samples[suppressWarnings(pick_sample_column(samples, sample_id, p))])
      }
      if (!startsWith(line, "##")) return(NA_character_)
    }
  }
  chosen = unique(unlist(lapply(vcf_files, read_sample)))
  chosen = chosen[!is.na(chosen)]

  list(
    paths   = unname(vapply(vcf_files, rel, character(1))),
    sources = unique(unlist(lapply(vcf_files, read_sources))),
    sample  = chosen
  )
}

# ClairS-TO Verdict status file (<id>.verdict_status.tsv: status, purity, ploidy, purity_source,
# message; one data row). NULL when absent or unreadable; columns are read by name
parse_verdict_status = function(path) {
  if (is.null(path) || length(path) == 0 || !file.exists(path) || isTRUE(file.size(path) == 0)) return(NULL)
  dt = tryCatch(fread(path, sep = "\t", header = TRUE, colClasses = "character", na.strings = c("NA", ""),
                      quote = ""),
                error = function(e) { message("Verdict status: cannot read ", path, ": ", conditionMessage(e)); NULL })
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  setnames(dt, tolower(trimws(names(dt))))
  get1 = function(k) if (k %in% names(dt)) as.character(dt[[k]][1]) else NA_character_
  status = tolower(trimws(get1("status")))
  if (is.na(status) || !nzchar(status)) return(NULL)
  list(status = status,
       purity = suppressWarnings(as.numeric(get1("purity"))),
       ploidy = suppressWarnings(as.numeric(get1("ploidy"))),
       purity_source = tolower(trimws(get1("purity_source"))),
       message = trimws(get1("message")))
}

VERDICT_SOURCE_LABELS = c(ascat = "ASCAT", wakhan = "Wakhan", savana = "SAVANA",
                          clairsto_internal = "ClairS-TO's own estimate", none = "no estimate")

# Reader-facing sentence for a Verdict status; NULL when tagging was applied (or the status is unknown)
verdict_note = function(v) {
  if (is.null(v) || is.null(v$status)) return(NULL)
  if (identical(v$status, "skipped")) {
    msg = if (!is.na(v$message) && nzchar(v$message))
            sub("[.!]+$", "", sub("^\\[[A-Z]+\\]\\s*", "", v$message)) else "no reason recorded"
    src = if (!is.na(v$purity_source) && v$purity_source %in% names(VERDICT_SOURCE_LABELS))
            VERDICT_SOURCE_LABELS[[v$purity_source]] else if (!is.na(v$purity_source)) v$purity_source else "unknown source"
    pur = if (!is.na(v$purity)) formatC(v$purity, format = "f", digits = 2) else "unknown"
    return(paste0("ClairS-TO germline tagging (Verdict) was skipped: ", msg, "; purity ", pur, " from ", src,
                  ". Tumour-only calls likely contain many germline variants."))
  }
  if (identical(v$status, "not_run"))
    return("ClairS-TO germline tagging (Verdict) did not run for this sample, so tumour-only calls were not screened for germline variants.")
  NULL
}

# Canonical variant key joining VEP rows to VCF records: VEP reports indels at anchor + 1 with raw or dash-trimmed alleles; all forms reconcile as trimmed alleles at anchor + 1, SNVs/MNVs verbatim. `space` is "vcf" or "vep"
variant_key = function(chrom, pos, ref, alt, space = c("vcf", "vep")) {
  space = match.arg(space)
  ref = toupper(as.character(ref)); alt = toupper(as.character(alt))
  pos = as.integer(pos)

  trim = function(x) { t = substr(x, 2L, nchar(x)); ifelse(t == "", "-", t) }

  dash     = ref == "-" | alt == "-"   # already trimmed by VEP
  is_indel = dash | nchar(ref) != nchar(alt)

  # Raw allele pairs still need the anchor base dropped; dash forms are already trimmed.
  need_trim = is_indel & !dash

  # Only the VCF side needs shifting — VEP has already done it.
  key_pos = ifelse(space == "vcf" & is_indel, pos + 1L, pos)
  key_ref = ifelse(need_trim, trim(ref), ref)
  key_alt = ifelse(need_trim, trim(alt), alt)

  paste(chrom, key_pos, key_ref, key_alt, sep = "|")
}

# Classify SNV into 6 SBS mutation categories (C/T-ref normalised)
classify_mut = function(ref, alt) {
  comp = c(A = "T", T = "A", C = "G", G = "C")
  ref = toupper(ref); alt = toupper(alt)
  use_comp = !(ref %in% c("C", "T"))
  norm_ref = ifelse(use_comp, comp[ref], ref)
  norm_alt = ifelse(use_comp, comp[alt], alt)
  paste0(norm_ref, ">", norm_alt)
}

# Small-variant display table: VEP rows with VAF/depth/phasing joined from the annotated VCF; gene_panel filters symbols (NULL = all)
build_variant_table = function(vep_data, vaf_data, gene_panel = NULL) {
  if (is.null(vep_data) || nrow(vep_data) == 0) return(NULL)
  vep_fields = attr(vep_data, "vep_fields")   # merge()/unique() below drop attributes

  # Impact ranking for deduplication
  impact_rank = c(HIGH = 1L, MODERATE = 2L, LOW = 3L, MODIFIER = 4L)
  vep_data[, impact_rank := impact_rank[impact]]
  vep_data[is.na(impact_rank), impact_rank := 5L]

  # Filter to gene panel (by gene symbol or Ensembl ID fallback)
  if (!is.null(gene_panel)) {
    if (length(gene_panel) > 0) {
      vep_data = vep_data[symbol %in% gene_panel | gene_id %in% gene_panel]
    } else {
      vep_data = vep_data[FALSE]  # Empty panel → empty result
    }
  }
  if (nrow(vep_data) == 0) return(data.table())

  # Keep best consequence per variant×gene (lowest impact rank)
  key_cols = c("chrom", "pos", "ref", "alt", "symbol")
  setorder(vep_data, impact_rank)
  vep_data = unique(vep_data, by = key_cols)

  # Join on the canonical key; the coordinate space is declared by each parser (text = vep, CSQ = vcf), never assumed
  vep_space = if ("coord_space" %in% names(vep_data)) unique(vep_data$coord_space) else "vep"
  stopifnot(length(vep_space) == 1L, vep_space %in% c("vep", "vcf"))
  vep_data[, join_key := variant_key(chrom, pos, ref, alt, space = vep_space)]

  if (!is.null(vaf_data) && nrow(vaf_data) > 0) {
    if (!"vcf_caller" %in% names(vaf_data)) vaf_data = copy(vaf_data)[, vcf_caller := NA_character_]
    vdt = vaf_data[, .(join_key = variant_key(chrom, pos, ref, alt, space = "vcf"),
                       vaf, dp, gt, ps, vcf_caller)]
    vdt = unique(vdt, by = "join_key")
    vep_data = merge(vep_data, vdt, by = "join_key", all.x = TRUE)
  } else {
    vep_data[, `:=`(vaf = NA_real_, dp = NA_integer_,
                    gt  = NA_character_, ps = NA_character_)]
  }

  # Which caller reported each variant (INFO/CALLER); empty on the VEP text path
  if ("caller" %in% names(vep_data) && any(!is.na(vep_data$caller))) {
    vep_data[, callers := caller]
  } else {
    vep_data[, callers := ""]
  }

  # Mutation category for SNVs
  vep_data[nchar(ref) == 1 & nchar(alt) == 1,
           mut_cat := classify_mut(ref, alt)]

  # Annotation columns exist only when the VEP run declared their source (see add_vep_annotations())
  display_cols = c("symbol", "chrom", "pos", "ref", "alt",
                   "consequence", "impact", "hgvsp",
                   "vaf", "dp", "gt", "ps",
                   "callers", "cosmic", "dbsnp", names(VEP_ANNOTATION_FIELDS))
  display_cols = display_cols[display_cols %in% names(vep_data)]
  out = vep_data[, ..display_cols]
  setattr(out, "vep_fields", vep_fields)
  # Rows whose GT was blanked as DeepSomatic's (parse_caller_vcf()), for the provenance footnote
  setattr(out, "n_gt_deepsomatic",
          if ("vcf_caller" %in% names(vep_data)) sum(vep_data$vcf_caller %in% "deepsomatic") else 0L)
  out
}

# ClinVar significance tokens that count as pathogenic, lower-cased. Exact tokens, never a
# substring test: "Conflicting_classifications_of_pathogenicity" contains "pathogenic"
CLINVAR_PATHOGENIC_TOKENS = c("pathogenic", "likely_pathogenic")

# Is each ClinVar significance cell pathogenic / likely pathogenic? Vectorised, NA-safe
# (NA/"" -> FALSE). Cells are token lists: add_vep_annotations() rewrites VEP's "&" to ",",
# ClinVar itself joins a combined call with "/" ("Pathogenic/Likely_pathogenic") and its
# secondary terms with "|" or ",_" ("Pathogenic,_low_penetrance"), so all four separate
# tokens and a leading "_" is dropped. Case-insensitive: the cache's CLIN_SIG is lower-case
# ("likely_pathogenic"), the ClinVar plugin's ClinVar_CLNSIG capitalised
clinvar_is_pathogenic = function(x) {
  if (length(x) == 0) return(logical(0))
  x = tolower(as.character(x))
  x[is.na(x)] = ""
  vapply(strsplit(x, "[,&/|]+"), function(tok) {
    tok = sub("^[_[:space:]]+", "", trimws(tok))
    any(tok %in% CLINVAR_PATHOGENIC_TOKENS)
  }, logical(1))
}

# Upper-cased SNV-scope symbols of every panel in `all_panels` (a named list of panel lists
# after the YAML round trip, so every vector is unlist()ed). A panel with an applies_to
# column (has_scopes) contributes genes_snv only, even when that is empty (YAML NULL):
# falling back to `genes` there would pull its SV-only genes into the SNV slim set
panel_snv_genes = function(all_panels) {
  if (is.null(all_panels) || length(all_panels) == 0) return(character(0))
  g = unlist(lapply(all_panels, function(p) {
    if (!is.list(p)) return(NULL)
    if (isTRUE(as.logical(unlist(p$has_scopes)))) unlist(p$genes_snv)
    else if (!is.null(p$genes_snv)) unlist(p$genes_snv)
    else unlist(p$genes)
  }), use.names = FALSE)
  g = toupper(trimws(as.character(g)))
  unique(g[!is.na(g) & nzchar(g)])
}

# Rows of the small-variant table embedded in the HTML.
#
# mode "full" returns every row; mode "slim" (also the default for NULL) keeps rows that are
# impact HIGH/MODERATE, OR ClinVar pathogenic/likely pathogenic (clinvar_is_pathogenic() on
# the `clinvar` column, when the run carried one), OR whose symbol is in any panel's SNV
# gene set (panel_snv_genes(); case-insensitive). Columns are never dropped or renamed (JS
# resolves them by name via SNV_COLS). The result keeps `vep_fields` and every other custom
# attribute of `vt` (row subsetting drops them, so they are copied back), plus
# attr(, "n_full") = nrow(vt) and attr(, "n_embedded") = rows kept. `vt` is never modified
# by reference. NULL in, NULL out.
slim_variant_table = function(vt, all_panels, mode = "slim") {
  if (is.null(vt)) return(NULL)
  if (is.null(mode) || length(mode) == 0 || is.na(mode[1])) mode = "slim"
  mode = mode[1]
  if (!mode %in% c("slim", "full")) stop("slim_variant_table(): unknown mode '", mode, "'")
  src = vt
  vt = as.data.table(vt)

  if (mode == "full" || nrow(vt) == 0) {
    out = copy(vt)
  } else {
    keep = rep(FALSE, nrow(vt))
    if ("impact" %in% names(vt))
      keep = keep | (!is.na(vt$impact) & vt$impact %in% c("HIGH", "MODERATE"))
    if ("clinvar" %in% names(vt))
      keep = keep | clinvar_is_pathogenic(vt$clinvar)
    genes = panel_snv_genes(all_panels)
    if (length(genes) > 0 && "symbol" %in% names(vt))
      keep = keep | (!is.na(vt$symbol) & toupper(vt$symbol) %in% genes)
    out = vt[keep]   # a new table: `vt` is untouched
  }

  # Custom attributes (vep_fields, ...) do not survive `[`; data.table's own are left alone
  skip = c("names", "row.names", "class", ".internal.selfref", "sorted", "index",
           "n_full", "n_embedded")
  for (a in setdiff(names(attributes(src)), skip)) setattr(out, a, attr(src, a))
  setattr(out, "n_full", nrow(vt))
  setattr(out, "n_embedded", nrow(out))
  out
}

# Window of the genome-wide circos mutation-density track, in bp
CIRCOS_SNV_BIN = 1e6

# SNV counts per `bin`-bp window for the circos density track. Pure (no circlize) so it is
# testable. `dt` needs chrom/pos/ref/alt; only single-base substitutions are counted, and a
# duplicated chrom/pos/ref/alt once. `chrom_lens` is a named vector of chromosome lengths
# and fixes which chromosomes appear and their order. Every window is returned, empty ones
# with n = 0, so the track draws as one continuous area per chromosome; the last window of a
# chromosome ends at its length. Returns data.table(chrom, start, end, n, per_mb), 0-based
# half-open [start, end); per_mb = n * 1e6 / bin (scaled to a full window, so a short final
# window is not inflated into a spike).
bin_snv_density = function(dt, chrom_lens, bin = CIRCOS_SNV_BIN) {
  empty = data.table(chrom = character(), start = numeric(), end = numeric(),
                     n = integer(), per_mb = numeric())
  if (is.null(chrom_lens) || length(chrom_lens) == 0) return(empty)
  bins = rbindlist(lapply(names(chrom_lens), function(chr) {
    len = as.numeric(chrom_lens[[chr]])
    if (is.na(len) || len <= 0) return(NULL)
    st = seq(0, len - 1, by = bin)
    data.table(chrom = chr, start = st, end = pmin(st + bin, len))
  }))
  if (nrow(bins) == 0) return(empty)
  bins[, idx := as.integer(start %/% bin)]
  bins[, n := 0L]

  if (!is.null(dt) && nrow(dt) > 0 && all(c("chrom", "pos", "ref", "alt") %in% names(dt))) {
    snv = unique(data.table(chrom = as.character(dt$chrom), pos = as.numeric(dt$pos),
                            ref = as.character(dt$ref), alt = as.character(dt$alt)))
    snv = snv[nchar(ref) == 1L & nchar(alt) == 1L & ref != "-" & alt != "-" &
              chrom %in% names(chrom_lens) & !is.na(pos)]
    # 1-based VCF positions: window k holds positions k*bin+1 .. (k+1)*bin
    if (nrow(snv) > 0) {
      counts = snv[, .(cnt = .N), by = .(chrom, idx = as.integer((pos - 1) %/% bin))]
      bins[counts, n := i.cnt, on = .(chrom, idx)]
    }
  }
  bins[, per_mb := n * 1e6 / bin]
  bins[, idx := NULL]
  bins[]
}

# Y-axis ceiling of the density track: the `q` quantile of the non-empty windows, rounded
# up to an integer, never below 1 and never above the real maximum (rounded up). Taken over
# non-empty windows so a low-burden sample, where most windows are 0, still gets a scale
# set by the windows that have mutations; a kataegis window above it is clipped
circos_density_cap = function(per_mb, q = 0.99) {
  v = per_mb[!is.na(per_mb) & per_mb > 0]
  if (length(v) == 0) return(1)
  cap = ceiling(as.numeric(stats::quantile(v, q, names = FALSE)))
  max(1, min(cap, ceiling(max(v))))
}

# Full small-variant table as a gzipped TSV sidecar next to the report; list columns are
# dropped (fwrite cannot write them flat). Returns the path invisibly, NULL when nothing written.
write_full_variant_table = function(vt, path) {
  if (is.null(vt)) return(invisible(NULL))
  vt = as.data.table(vt)
  list_cols = names(vt)[vapply(vt, is.list, logical(1))]
  if (length(list_cols) > 0) vt = vt[, setdiff(names(vt), list_cols), with = FALSE]
  data.table::fwrite(vt, path, sep = "\t", compress = "gzip")
  invisible(path)
}
