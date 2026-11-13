suppressPackageStartupMessages({
  library(data.table)
})

# Parse mosdepth summary; `keep_chroms` avoids a bare "^chr" match (admits decoys, empty on chr-less references)
parse_mosdepth_summary = function(summary_file, keep_chroms = NULL) {
  empty = list(mean_depth = NA_real_, total_length = NA_real_, total_bases = NA_real_, table = data.table())
  if (is.null(summary_file) || !file.exists(summary_file)) return(empty)
  # A 0-byte file (stub runs, interrupted publish) reads as a column-less table; guard before indexing `chrom`
  dt = tryCatch(suppressWarnings(fread(summary_file, sep = "\t", header = TRUE)), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(empty)
  setnames(dt, tolower(names(dt)))
  if (!all(c("chrom", "mean") %in% names(dt))) return(empty)
  total_row = dt[chrom == "total"]
  mean_depth = if (nrow(total_row) > 0) total_row$mean[1] else NA_real_

  # Keep per-chromosome rows (exclude region-level and total)
  chr_rows = dt[chrom != "total" & !grepl("_region", chrom)]
  if (!is.null(keep_chroms) && length(keep_chroms) > 0) {
    chr_rows = chr_rows[ensure_chr_prefix(chrom) %in% ensure_chr_prefix(keep_chroms)]
  } else {
    chr_rows = chr_rows[grepl("^chr", chrom)]
  }
  total_length = if (nrow(total_row) > 0 && "length" %in% names(dt)) total_row$length[1] else NA_real_
  total_bases  = if (nrow(total_row) > 0 && "bases" %in% names(dt))  total_row$bases[1]  else NA_real_
  list(mean_depth = round(mean_depth, 2), total_length = total_length, total_bases = total_bases, table = chr_rows)
}

# Parse mosdepth global distribution (<prefix>[.mosdepth].global.dist.txt)
parse_mosdepth_dist = function(dist_file) {
  if (is.null(dist_file) || !file.exists(dist_file)) return(NULL)
  dt = tryCatch(suppressWarnings(fread(dist_file, sep = "\t", header = FALSE)), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0 || ncol(dt) != 3) return(NULL)
  setnames(dt, c("chrom", "coverage", "fraction"))
  dt
}

# Cramino keys, looked up EXACTLY (key = text before the first tab, trimmed, case-insensitive).
# cramino 2.0 writes "<metric> aligned" twins (N50 aligned, N75 aligned, ...) after the full-read-length
# rows, so an unanchored grepl("N50") is right only by row order. Alternatives are tried in order:
# "% from total" is "alignments" in cramino 1.x and "reads" in 2.0.
CRAMINO_KEYS = list(
  n_alignments        = "Number of alignments",
  mapped_pct          = c("% from total reads", "% from total alignments"),
  n_reads             = "Number of reads",
  yield_gb            = "Yield [Gb]",
  mean_coverage       = "Mean coverage",
  n50                 = "N50",
  n75                 = "N75",
  median_length       = "Median length",
  mean_length         = "Mean length",
  n50_aligned         = "N50 aligned",
  n75_aligned         = "N75 aligned",
  median_length_aligned = "Median length aligned",
  mean_length_aligned = "Mean length aligned",
  median_identity     = c("Median identity", "Median estimated identity")
)

# Parse cramino alignment report. Every CRAMINO_KEYS field is present (NA when absent): the *_aligned
# fields are NA on cramino 1.x files, which is how a template tells the two formats apart.
parse_cramino = function(cramino_file) {
  out = lapply(CRAMINO_KEYS, function(k) NA_real_)
  out$n_reads = NA_integer_
  if (is.null(cramino_file) || !file.exists(cramino_file)) return(out)

  lines = readLines(cramino_file, warn = FALSE)
  lines = lines[grepl("\t", lines, fixed = TRUE)]
  if (length(lines) == 0) return(out)
  keys = tolower(trimws(sub("\t.*$", "", lines)))
  vals = trimws(sub("^[^\t]*\t", "", lines))

  get_field = function(alternatives) {
    for (k in tolower(alternatives)) {
      hit = which(keys == k)
      if (length(hit) > 0) return(vals[hit[1]])
    }
    NA_character_
  }
  for (f in names(CRAMINO_KEYS)) {
    out[[f]] = suppressWarnings(as.numeric(sub("%", "", get_field(CRAMINO_KEYS[[f]]), fixed = TRUE)))
  }
  out$n_reads = suppressWarnings(as.integer(out$n_reads))
  out
}

# samtools flagstat labels, matched EXACTLY on the text between "N + M " and any " (...)" suffix.
# A substring match is ambiguous on samtools >= 1.13: "mapped (" also hits "primary mapped (", and
# "duplicate" hits "primary duplicates", so the later line silently won. Names are the list keys;
# order is display order. Older samtools simply lacks the primary* lines.
FLAGSTAT_LABELS = c(
  total              = "in total",
  primary            = "primary",
  secondary          = "secondary",
  supplementary      = "supplementary",
  duplicate          = "duplicates",
  primary_duplicates = "primary duplicates",
  mapped             = "mapped",
  primary_mapped     = "primary mapped",
  paired             = "paired in sequencing"
)

# Parse samtools flagstat: QC-passed counts, keyed as in FLAGSTAT_LABELS (absent labels omitted)
parse_flagstat = function(flagstat_file) {
  if (is.null(flagstat_file) || !file.exists(flagstat_file)) return(list())
  lines = readLines(flagstat_file, warn = FALSE)
  re = "^\\s*([0-9]+) \\+ [0-9]+ (.*)$"
  lines = lines[grepl(re, lines)]
  if (length(lines) == 0) return(list())
  counts = suppressWarnings(as.numeric(sub(re, "\\1", lines)))
  labels = trimws(sub(" \\(.*$", "", sub(re, "\\2", lines)))
  out = list()
  for (k in names(FLAGSTAT_LABELS)) {
    hit = which(labels == FLAGSTAT_LABELS[[k]])
    if (length(hit) > 0) out[[k]] = counts[hit[1]]
  }
  out
}

# Parse samtools stats (*.stats), SN summary lines only
parse_samtools_stats = function(stats_file) {
  if (is.null(stats_file) || !file.exists(stats_file)) return(NULL)
  lines = readLines(stats_file, warn = FALSE)
  sn = lines[startsWith(lines, "SN\t")]
  if (length(sn) == 0) return(NULL)
  # Anchored on "SN\t<key>:\t", so "reads mapped" cannot collide with "reads mapped and paired"
  get_sn = function(key) {
    hit = grep(paste0("^SN\t", key, ":\t"), sn, value = TRUE)
    if (length(hit) == 0) return(NA_real_)
    suppressWarnings(as.numeric(trimws(sub(paste0("^SN\t", key, ":\t([^\t#]+).*"), "\\1", hit[1]))))
  }
  reads_total  = get_sn("raw total sequences")
  reads_mapped = get_sn("reads mapped")
  mapped_pct   = if (!is.na(reads_total) && reads_total > 0)
    round(reads_mapped / reads_total * 100, 2) else NA_real_
  list(
    reads_total  = reads_total,
    reads_mapped = reads_mapped,
    mapped_pct   = mapped_pct,
    total_length = get_sn("total length"),
    bases_mapped = get_sn("bases mapped \\(cigar\\)"),
    error_rate   = get_sn("error rate"),
    avg_length   = get_sn("average length"),
    max_length   = get_sn("maximum length"),
    avg_quality  = get_sn("average quality")
  )
}
