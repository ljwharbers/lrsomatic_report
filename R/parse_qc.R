suppressPackageStartupMessages({
  library(data.table)
})

# Parse mosdepth summary (*.mosdepth.summary.txt)
# Returns list: mean_depth, total_row (the "total" row from mosdepth)
parse_mosdepth_summary = function(summary_file) {
  if (is.null(summary_file) || !file.exists(summary_file)) {
    return(list(mean_depth = NA_real_, table = data.table()))
  }
  dt = fread(summary_file, sep = "\t", header = TRUE)
  setnames(dt, tolower(names(dt)))
  total_row = dt[chrom == "total"]
  mean_depth = if (nrow(total_row) > 0) total_row$mean[1] else NA_real_

  # Keep per-chromosome rows (exclude region-level and total)
  chr_rows = dt[grepl("^chr", chrom) & !grepl("_region", chrom)]
  list(mean_depth = round(mean_depth, 2), table = chr_rows)
}

# Parse mosdepth global distribution (*.mosdepth.global.dist.txt)
# Returns data.table with columns: chrom, coverage, fraction
parse_mosdepth_dist = function(dist_file) {
  if (is.null(dist_file) || !file.exists(dist_file)) return(NULL)
  dt = fread(dist_file, sep = "\t", header = FALSE,
             col.names = c("chrom", "coverage", "fraction"))
  dt
}

# Parse cramino alignment report
# Returns list: n50, yield_gb, mapped_pct, n_reads
parse_cramino = function(cramino_file) {
  if (is.null(cramino_file) || !file.exists(cramino_file)) {
    return(list(n50 = NA_real_, yield_gb = NA_real_,
                mapped_pct = NA_real_, n_reads = NA_integer_))
  }
  lines = readLines(cramino_file, warn = FALSE)
  get_val = function(pattern) {
    hit = grep(pattern, lines, value = TRUE, ignore.case = TRUE)
    if (length(hit) == 0) return(NA_character_)
    trimws(sub(paste0(".*", pattern, "\\s*"), "", hit[1], ignore.case = TRUE))
  }

  # Cramino outputs key\tvalue pairs
  dt = tryCatch(
    fread(cramino_file, sep = "\t", header = FALSE, col.names = c("key", "value"), fill = TRUE),
    error = function(e) NULL
  )
  if (is.null(dt)) return(list(n50 = NA_real_, yield_gb = NA_real_,
                               mapped_pct = NA_real_, n_reads = NA_integer_))

  get_field = function(pattern) {
    row = dt[grepl(pattern, key, ignore.case = TRUE)]
    if (nrow(row) == 0) NA_character_ else as.character(row$value[1])
  }

  list(
    n50        = suppressWarnings(as.numeric(get_field("N50"))),
    yield_gb   = suppressWarnings(as.numeric(get_field("Yield"))),
    mapped_pct = suppressWarnings(as.numeric(sub("%", "", get_field("% from total")))),
    n_reads    = suppressWarnings(as.integer(get_field("Number of reads")))
  )
}

# Parse samtools flagstat
# Returns a named list of counts (total, mapped, ...)
parse_flagstat = function(flagstat_file) {
  if (is.null(flagstat_file) || !file.exists(flagstat_file)) return(list())
  lines = readLines(flagstat_file, warn = FALSE)
  out = list()
  for (line in lines) {
    count = suppressWarnings(as.integer(sub(" .*", "", trimws(line))))
    if (grepl("in total", line))    out$total    = count
    if (grepl("mapped \\(", line))  out$mapped   = count
    if (grepl("paired in seq", line)) out$paired = count
    if (grepl("secondary", line))   out$secondary = count
    if (grepl("supplementary", line)) out$supplementary = count
    if (grepl("duplicate", line))   out$duplicate = count
  }
  out
}
