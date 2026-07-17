suppressPackageStartupMessages({
  library(data.table)
})

# Parse ASCAT raw segments (segments_raw.txt)
# Columns: sample, chr, startpos, endpos, nMajor, nMinor, nAraw, nBraw
parse_ascat_segments = function(segments_file) {
  if (is.null(segments_file) || !file.exists(segments_file)) return(NULL)
  dt = fread(segments_file, sep = "\t", header = TRUE)

  # Normalise column names to lowercase
  setnames(dt, tolower(names(dt)))

  # Add chr prefix if missing
  dt[, chr := ensure_chr_prefix(as.character(chr))]

  # Column names after tolower(): naraw, nbraw
  dt[, total_cn := pmin(naraw + nbraw, 4)]
  dt[, major_cn := pmin(naraw, 4)]
  dt[, minor_cn := pmin(nbraw, 4)]

  dt
}

# Parse ASCAT purity/ploidy file
# Columns: AberrantCellFraction, Ploidy
parse_ascat_purityploidy = function(pp_file) {
  if (is.null(pp_file) || !file.exists(pp_file)) return(list(purity = NA_real_, ploidy = NA_real_))
  dt = fread(pp_file, sep = "\t", header = TRUE)
  setnames(dt, tolower(names(dt)))
  list(
    purity = round(as.numeric(dt$aberrantcellfraction[1]), 3),
    ploidy = round(as.numeric(dt$ploidy[1]), 3)
  )
}

# Parse Wakhan's ranked purity/ploidy solutions table (wakhan/solutions_ranks.tsv).
# Columns: repository_name, dna_purity, cell_purity, ploidy, confidence, solution_rank
parse_wakhan_solutions = function(tsv_file) {
  if (is.null(tsv_file) || !file.exists(tsv_file)) return(NULL)
  dt = fread(tsv_file, sep = "\t", header = TRUE)
  if (nrow(dt) == 0) return(NULL)
  setorder(dt, solution_rank)
  dt
}
