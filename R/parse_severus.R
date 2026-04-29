suppressPackageStartupMessages({
  library(data.table)
})

# Parse Severus somatic VCF for circos plot data.
# Returns list: $translocations (BND records) and $nontrans (DEL/DUP/INV/INS)
parse_severus_vcf = function(vcf_file) {
  if (is.null(vcf_file) || !file.exists(vcf_file)) {
    return(list(translocations = data.table(), nontrans = data.table()))
  }

  con = gzcon(file(vcf_file, "rb"))
  skip_n = 0L
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#CHROM")) break
    skip_n = skip_n + 1L
  }
  close(con)

  dt = fread(vcf_file, skip = skip_n + 1L, sep = "\t", header = FALSE,
             select = 1:8,
             col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
  if (nrow(dt) == 0) {
    return(list(translocations = data.table(), nontrans = data.table()))
  }

  dt[, CHROM := ensure_chr_prefix(CHROM)]

  # Extract INFO sub-fields
  .info_val = function(info_vec, key) {
    pattern = paste0("(?:^|;)", key, "=([^;]+)")
    m = regmatches(info_vec, regexpr(pattern, info_vec, perl = TRUE))
    ifelse(nchar(m) > 0, sub(paste0(".*="), "", m), NA_character_)
  }

  dt[, SVTYPE := .info_val(INFO, "SVTYPE")]
  dt[grepl("END=",   INFO, fixed = TRUE), END   := as.integer(.info_val(INFO[grepl("END=",   INFO, fixed = TRUE)], "END"))]
  dt[grepl("SVLEN=", INFO, fixed = TRUE), SVLEN := as.integer(.info_val(INFO[grepl("SVLEN=", INFO, fixed = TRUE)], "SVLEN"))]

  # BND partner chromosome/position from ALT field
  # ALT format examples: "N[chr7:24547089[" or "]chr7:24547089]N"
  dt[SVTYPE == "BND", CHROM2 := {
    m = regmatches(ALT, regexpr("chr[^:]+", ALT, perl = TRUE))
    ifelse(nchar(m) > 0, m, NA_character_)
  }]
  dt[SVTYPE == "BND", POS2 := as.integer(regmatches(ALT, regexpr("(?<=:)\\d+", ALT, perl = TRUE)))]

  # Insertions have no END — use POS
  dt[SVTYPE == "INS" | is.na(END), END := POS]

  # Colours and y-positions for non-BND SV track
  SV_COL  = c(INS = "#f97e02", DEL = "#020272", INV = "#e7cc02", DUP = "#e41a1c")
  SV_YPOS = c(INS = 1.0,       DEL = 0.66,      INV = 0.33,      DUP = 0.05)
  dt[SVTYPE %in% names(SV_COL),  circos_col := SV_COL[SVTYPE]]
  dt[SVTYPE %in% names(SV_YPOS), circos_pos := SV_YPOS[SVTYPE]]

  translocations = dt[SVTYPE == "BND" & !is.na(CHROM2) & !is.na(POS2),
    .(chrom = CHROM, pos = POS, chrom2 = CHROM2, pos2 = POS2)]

  nontrans = dt[SVTYPE != "BND",
    .(chrom = CHROM, pos = POS, end = END, svtype = SVTYPE,
      svlen = SVLEN, circos_pos, circos_col)]

  list(translocations = translocations, nontrans = nontrans)
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

# Build the SV display table.
# gene_panel: character vector of HGNC symbols to keep, or NULL to return all SVs (one row each).
build_sv_table = function(sv_tsv, gene_panel = NULL) {
  if (is.null(sv_tsv) || nrow(sv_tsv) == 0) return(data.table())

  sv_tsv = copy(sv_tsv)

  # When no panel is supplied, return one row per SV without explosion
  if (is.null(gene_panel)) {
    display_cols = intersect(
      c("ID", "SVTYPE", "DETAILED_TYPE",
        "START_CHROM", "START_POS", "END_CHROM", "END_POS",
        "SV_LEN", "VAF", "NHL_GENE_HITS", "COSMIC_GENE_HITS",
        "NHL_NEAREST_GENE_HITS_1MBWINDOW"),
      names(sv_tsv)
    )
    return(sv_tsv[, ..display_cols])
  }

  # Panel-filtered path: explode multi-gene gene_hits, filter, return one row per gene×SV
  sv_tsv[, .ridx := .I]

  sv_long = sv_tsv[, {
    raw = as.character(gene_hits[1])
    genes = unique(trimws(unlist(strsplit(raw, "[;,]+"))))
    genes = genes[nchar(genes) > 0 & genes != "-" & toupper(genes) != "NA"]
    if (length(genes) == 0) genes = NA_character_
    list(gene = genes)
  }, by = .ridx]

  sv_long = merge(sv_long, sv_tsv, by = ".ridx")
  sv_long[, .ridx := NULL]
  sv_tsv[,  .ridx := NULL]

  # Filter by gene panel (always, even if panel is empty)
  sv_long = sv_long[!is.na(gene) & gene %in% gene_panel]
  if (nrow(sv_long) == 0) return(data.table())

  display_cols = intersect(
    c("gene", "ID", "SVTYPE", "DETAILED_TYPE",
      "START_CHROM", "START_POS", "END_CHROM", "END_POS",
      "SV_LEN", "VAF", "NHL_GENE_HITS", "COSMIC_GENE_HITS",
      "NHL_NEAREST_GENE_HITS_1MBWINDOW"),
    names(sv_long)
  )
  sv_long[, ..display_cols]
}
