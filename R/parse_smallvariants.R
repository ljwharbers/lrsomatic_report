suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# Derive dbsnp/cosmic columns from a VEP "Existing_variation" column (semicolon- or
# comma-joined list of IDs, e.g. "rs123&COSV456"). Shared by parse_vep_text/parse_vep_vcf.
derive_dbsnp_cosmic = function(dt) {
  dt[, dbsnp := sub("(rs[0-9]+).*", "\\1", existing)]
  dt[!grepl("^rs", dbsnp, perl = TRUE), dbsnp := NA_character_]

  dt[, cosmic := sub(".*(COS[VM][0-9]+).*", "\\1", existing)]
  dt[!grepl("^COS", cosmic, perl = TRUE), cosmic := NA_character_]
  dt
}

# Dispatch to the right VEP parser based on actual file contents — both forms ship as
# "*_SOMATIC_VEP.vcf.gz" so the filename alone doesn't tell you which one you have.
# - VEP default text output: "##"-commented header, column line starts with "#Uploaded_variation"
# - genuine VCF w/ CSQ INFO field: "##fileformat=VCFv4.2", column line starts with "#CHROM"
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

# Parse the VEP default text output (tab-delimited, ##-commented header, NOT a VCF).
# Returns a data.table with one row per consequence per variant.
parse_vep_text = function(vep_file) {
  if (is.null(vep_file) || !file.exists(vep_file)) return(NULL)

  # Count meta-lines (start with ##) to find the column-header line
  con = gzfile(vep_file, "rb")
  skip_n = 0L
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#Uploaded_variation")) break
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

  # Parse Location → chrom, pos (format: "chr1:3506" or "chr1:3506-3510")
  dt[, chrom := sub(":.*", "", Location)]
  dt[, pos   := as.integer(sub(".*:(\\d+).*", "\\1", Location))]
  dt[, chrom := ensure_chr_prefix(chrom)]

  # Parse variant_id → ref, alt ("chr1_3506_A/G")
  dt[, ref := sub(".*_([^/]+)/.*", "\\1", variant_id)]
  dt[, alt := sub(".*/", "", variant_id)]

  # Parse VEP Extra key=value field
  dt[, symbol   := extract_extra_key(Extra, "SYMBOL")]
  dt[, impact   := extract_extra_key(Extra, "IMPACT")]
  dt[, existing := extract_extra_key(Extra, "Existing_variation")]
  dt[, sift     := extract_extra_key(Extra, "SIFT")]
  dt[, polyphen := extract_extra_key(Extra, "PolyPhen")]
  dt[, hgvsp    := extract_extra_key(Extra, "HGVSp")]

  # dbSNP / COSMIC IDs, derived from Existing_variation
  dt = derive_dbsnp_cosmic(dt)

  dt
}

# Parse a genuine VCF carrying VEP annotation in a CSQ INFO field (VEP run with --vcf,
# as opposed to the default text output handled by parse_vep_text()).
# Returns a data.table with one row per gene/transcript annotation per variant, using the
# same column contract as parse_vep_text(): chrom, pos, ref, alt, symbol, gene_id,
# consequence, impact, hgvsp, existing, dbsnp, cosmic, sift, polyphen.
parse_vep_vcf = function(vep_file) {
  if (is.null(vep_file) || !file.exists(vep_file)) return(NULL)

  # Skip header to #CHROM, capturing the CSQ field order from its INFO meta-line
  # (e.g. "...Format: Allele|Consequence|IMPACT|SYMBOL|Gene|...")
  con = gzfile(vep_file, "rb")
  skip_n = 0L
  csq_format = NULL
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "##INFO=<ID=CSQ")) {
      m = regmatches(line, regexpr("Format: [^\"]+", line))
      if (length(m) > 0) csq_format = strsplit(sub("^Format: ", "", m), "|", fixed = TRUE)[[1]]
    }
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

  # Fixed-string splits (no regex backtracking) — much faster than sub() on the full
  # semicolon-delimited INFO field across millions of rows.
  dt[, CSQ := tstrsplit(INFO, "CSQ=", fixed = TRUE, keep = 2L)[[1]]]
  dt[, CSQ := tstrsplit(CSQ,  ";",    fixed = TRUE, keep = 1L)[[1]]]

  # One row per gene/transcript annotation (comma-separated CSQ entries)
  dt_long = dt[, .(csq_entry = unlist(strsplit(CSQ, ",", fixed = TRUE))),
               by = .(CHROM, POS, REF, ALT)]
  if (nrow(dt_long) == 0) return(NULL)

  # Split each entry on "|", materialising only the fields actually used below
  # (tstrsplit returns columns in the order requested via `keep`).
  need = c("Consequence", "IMPACT", "SYMBOL", "Gene", "HGVSp",
           "Existing_variation", "SIFT", "PolyPhen")
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
  dt_long[, sift        := get_field("SIFT")]
  dt_long[, polyphen    := get_field("PolyPhen")]

  # Blank fields are "" in CSQ, not NA — normalise for consistency with parse_vep_text()
  for (col in c("symbol", "impact", "hgvsp", "existing", "gene_id", "sift", "polyphen")) {
    dt_long[get(col) == "", (col) := NA_character_]
  }

  dt_long[, chrom := ensure_chr_prefix(CHROM)]
  dt_long[, pos   := as.integer(POS)]
  dt_long[, ref   := REF]
  dt_long[, alt   := ALT]

  # dbSNP / COSMIC IDs, derived from Existing_variation
  dt_long = derive_dbsnp_cosmic(dt_long)

  dt_long[, .(chrom, pos, ref, alt, symbol, gene_id, consequence, impact, hgvsp,
              existing, dbsnp, cosmic, sift, polyphen)]
}

# Parse a raw caller VCF for variant coordinates + VAF.
# Returns data.table with: chrom, pos, ref, alt, vaf, dp, caller
parse_caller_vcf = function(vcf_file, caller_name = "unknown") {
  if (is.null(vcf_file) || !file.exists(vcf_file)) return(NULL)

  # Count header lines
  con = gzfile(vcf_file, "rb")
  skip_n = 0L
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#CHROM")) break
    skip_n = skip_n + 1L
  }
  close(con)

  # Read up to 10 columns (standard VCF single-sample layout)
  col_names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE1")
  dt = fread(vcf_file, skip = skip_n + 1L, sep = "\t", header = FALSE,
             select = 1:10, col.names = col_names)
  if (nrow(dt) == 0) return(NULL)

  dt[, CHROM := ensure_chr_prefix(CHROM)]

  # Extract AF and DP from FORMAT + SAMPLE1 columns
  # Work format-group by format-group to avoid splitting every single row redundantly
  fmt_groups = unique(dt$FORMAT)
  vaf_list = vector("numeric", nrow(dt))
  dp_list  = vector("integer", nrow(dt))

  for (fmt in fmt_groups) {
    idx_rows = which(dt$FORMAT == fmt)
    fields = strsplit(fmt, ":", fixed = TRUE)[[1]]
    af_idx = match("AF", fields)
    dp_idx = match("DP", fields)
    samples = dt$SAMPLE1[idx_rows]
    split_s = strsplit(samples, ":", fixed = TRUE)

    if (!is.na(af_idx)) {
      vaf_list[idx_rows] = vapply(split_s, function(x)
        if (length(x) >= af_idx) suppressWarnings(as.numeric(x[af_idx])) else NA_real_,
        numeric(1))
    } else {
      vaf_list[idx_rows] = NA_real_
    }
    if (!is.na(dp_idx)) {
      dp_list[idx_rows] = vapply(split_s, function(x)
        if (length(x) >= dp_idx) suppressWarnings(as.integer(x[dp_idx])) else NA_integer_,
        integer(1))
    } else {
      dp_list[idx_rows] = NA_integer_
    }
  }

  dt[, vaf    := vaf_list]
  dt[, dp     := dp_list]
  dt[, caller := caller_name]
  dt[, .(chrom = CHROM, pos = POS, ref = REF, alt = ALT, vaf, dp, caller)]
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

# Build the small-variant display table:
#   canonical rows from VEP text, per-caller VAFs joined on position.
# gene_panel: character vector of HGNC symbols to keep, or NULL to return all variants.
build_variant_table = function(vep_data, caller_list, gene_panel = NULL) {
  if (is.null(vep_data) || nrow(vep_data) == 0) return(NULL)

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

  # Join key for caller data
  vep_data[, join_key := paste(chrom, pos, ref, alt, sep = "|")]

  # Left-join per-caller VAFs
  for (nm in names(caller_list)) {
    cdt = caller_list[[nm]]
    vaf_col = paste0("vaf_", nm)
    if (is.null(cdt) || nrow(cdt) == 0) {
      vep_data[, (vaf_col) := NA_real_]
      next
    }
    cdt_sub = cdt[, .(join_key = paste(chrom, pos, ref, alt, sep = "|"), vaf)]
    cdt_sub = unique(cdt_sub, by = "join_key")
    vep_data = merge(vep_data, cdt_sub, by = "join_key", all.x = TRUE)
    setnames(vep_data, "vaf", vaf_col)
  }

  # Summarise which callers reported each variant
  vaf_cols = paste0("vaf_", names(caller_list))
  vaf_cols = vaf_cols[vaf_cols %in% names(vep_data)]
  if (length(vaf_cols) > 0) {
    vep_data[, callers := {
      .sd = .SD
      apply(.sd, 1, function(r) paste(names(caller_list)[!is.na(r)], collapse = ","))
    }, .SDcols = vaf_cols]
  } else {
    vep_data[, callers := ""]
  }

  # Mutation category for SNVs
  vep_data[nchar(ref) == 1 & nchar(alt) == 1,
           mut_cat := classify_mut(ref, alt)]

  display_cols = c("symbol", "chrom", "pos", "ref", "alt",
                   "consequence", "impact", "hgvsp",
                   vaf_cols, "callers", "cosmic", "dbsnp", "sift", "polyphen")
  display_cols = display_cols[display_cols %in% names(vep_data)]
  vep_data[, ..display_cols]
}
