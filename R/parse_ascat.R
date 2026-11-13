suppressPackageStartupMessages({
  library(data.table)
})

# Parse ASCAT raw segments (segments_raw.txt)
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
parse_ascat_purityploidy = function(pp_file) {
  if (is.null(pp_file) || !file.exists(pp_file)) return(list(purity = NA_real_, ploidy = NA_real_))
  dt = fread(pp_file, sep = "\t", header = TRUE)
  setnames(dt, tolower(names(dt)))
  list(
    purity = round(as.numeric(dt$aberrantcellfraction[1]), 3),
    ploidy = round(as.numeric(dt$ploidy[1]), 3)
  )
}

# Parse Wakhan's ranked purity/ploidy solutions table (wakhan/solutions_ranks.tsv)
parse_wakhan_solutions = function(tsv_file) {
  if (is.null(tsv_file) || !file.exists(tsv_file)) return(NULL)
  dt = fread(tsv_file, sep = "\t", header = TRUE)
  if (nrow(dt) == 0) return(NULL)
  setorder(dt, solution_rank)
  dt
}

# ---------------------------------------------------------------------------
# Unified copy-number segments: data.table(chrom, start, end, cn_a, cn_b, total, minor,
# conf, sv_ids, conf_a, conf_b). cn_a/cn_b are HP1/HP2 for Wakhan and major/minor for
# ASCAT; coordinates as written by the tool (both are contiguous inclusive spans).
# Unlike parse_ascat_segments() nothing is capped — that one feeds the circos.
# ---------------------------------------------------------------------------

.empty_cn_segments = function() {
  data.table(chrom = character(), start = numeric(), end = numeric(),
             cn_a = numeric(), cn_b = numeric(), total = numeric(), minor = numeric(),
             conf = numeric(), sv_ids = character(), conf_a = numeric(), conf_b = numeric())
}

# Integer allele-specific segments from ASCAT's <id>.segments_raw.txt (nMajor/nMinor; the
# rounded nAraw/nBraw only when an older file lacks them). segments.txt is not staged.
parse_ascat_integer_segments = function(segments_file) {
  if (is.null(segments_file) || !file.exists(segments_file)) return(NULL)
  dt = fread(segments_file, sep = "\t", header = TRUE)
  if (nrow(dt) == 0) return(.empty_cn_segments())
  setnames(dt, tolower(names(dt)))
  maj = if ("nmajor" %in% names(dt)) as.numeric(dt$nmajor) else round(as.numeric(dt$naraw))
  mnr = if ("nminor" %in% names(dt)) as.numeric(dt$nminor) else round(as.numeric(dt$nbraw))
  out = data.table(chrom = ensure_chr_prefix(dt$chr),
                   start = as.numeric(dt$startpos), end = as.numeric(dt$endpos),
                   cn_a = maj, cn_b = mnr)
  out[, `:=`(total = cn_a + cn_b, minor = pmin(cn_a, cn_b), conf = NA_real_,
             sv_ids = NA_character_, conf_a = NA_real_, conf_b = NA_real_)]
  out[]
}

# ASCAT <id>.metrics.txt (one row): the genome-level fields the report uses, or NULL.
# WGD is ASCAT's own call (0/1, on the fitted profile), LOH its LOH genome fraction,
# GI its genomic-instability fraction.
parse_ascat_metrics = function(metrics_file) {
  if (is.null(metrics_file) || !file.exists(metrics_file)) return(NULL)
  dt = tryCatch(fread(metrics_file, sep = "\t", header = TRUE), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  num = function(k) if (k %in% names(dt)) suppressWarnings(as.numeric(dt[[k]][1])) else NA_real_
  list(wgd = num("WGD"), loh = num("LOH"), gi = num("GI"), purity = num("purity"),
       ploidy = num("ploidy"), goodness_of_fit = num("goodness_of_fit"),
       sex = if ("sex" %in% names(dt)) as.character(dt$sex[1]) else NA_character_)
}

# Read a Wakhan BED whose header is the LAST leading `#` line (the ones above it are a
# prose legend), with the `#` stripped. NULL when the file is absent.
read_hash_header_table = function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  con = file(path, "r")
  on.exit(close(con))
  hdr = character(0)
  n_hash = 0L
  repeat {
    l = readLines(con, n = 1, warn = FALSE)
    if (length(l) == 0 || !startsWith(l, "#")) break
    hdr = l
    n_hash = n_hash + 1L
  }
  if (length(hdr) == 0) return(NULL)
  cols = trimws(strsplit(sub("^#\\s*", "", hdr), "\t", fixed = TRUE)[[1]])
  dt = tryCatch(fread(path, skip = n_hash, header = FALSE, sep = "\t", fill = TRUE,
                      colClasses = "character"),
                error = function(e) data.table())
  if (nrow(dt) == 0) {
    dt = as.data.table(setNames(replicate(length(cols), character(), simplify = FALSE), cols))
    return(dt)
  }
  # A trailing empty svs_breakpoints_ids column may be dropped or kept by fread
  if (ncol(dt) > length(cols)) dt = dt[, seq_along(cols), with = FALSE]
  if (ncol(dt) < length(cols)) for (k in seq(ncol(dt) + 1L, length(cols))) dt[, (paste0("V", k)) := NA_character_]
  setnames(dt, cols)
  dt
}

# Wakhan solution_*/integer_profile.bed -> unified segments (cn_a = HP1, cn_b = HP2;
# conf = the lower of the two haplotype confidences)
parse_wakhan_integer_profile = function(path) {
  dt = read_hash_header_table(path)
  if (is.null(dt)) return(NULL)
  if (nrow(dt) == 0) return(.empty_cn_segments())
  num = function(k) if (k %in% names(dt)) suppressWarnings(as.numeric(dt[[k]])) else NA_real_
  out = data.table(chrom = ensure_chr_prefix(dt$chr), start = num("start"), end = num("end"),
                   cn_a = num("hp1_copynumber_state"), cn_b = num("hp2_copynumber_state"))
  out[, `:=`(total = cn_a + cn_b, minor = pmin(cn_a, cn_b))]
  out[, conf_a := num("hp1_confidence")]
  out[, conf_b := num("hp2_confidence")]
  out[, conf := pmin(conf_a, conf_b)]
  sv = if ("svs_breakpoints_ids" %in% names(dt)) as.character(dt$svs_breakpoints_ids) else NA_character_
  sv[!is.na(sv) & !nzchar(sv)] = NA_character_
  out[, sv_ids := sv]
  setcolorder(out, names(.empty_cn_segments()))
  out[]
}

# Wakhan solution_*/genes/genes_copynumber_states.bed -> list(genes, n_dropped). Rows whose
# haplotype state exceeds `max_cn` are artefacts (state 394 with ~3300x "adjusted
# coverage" on P15) and are dropped and counted. NULL when the file is absent.
parse_wakhan_genes = function(path, max_cn = 100) {
  dt = read_hash_header_table(path)
  if (is.null(dt)) return(NULL)
  empty = data.table(chrom = character(), start = numeric(), end = numeric(),
                     gene = character(), hp1_state = numeric(), hp2_state = numeric(),
                     total = numeric(), minor = numeric())
  if (nrow(dt) == 0) return(list(genes = empty, n_dropped = 0L))
  cn = resolve_columns(dt, list(hp1 = c("hp1_state", "copynumber_state_hp1"),
                                hp2 = c("hp2_state", "copynumber_state_hp2")))
  num = function(k) if (!is.na(k)) suppressWarnings(as.numeric(dt[[k]])) else rep(NA_real_, nrow(dt))
  out = data.table(chrom = ensure_chr_prefix(dt$chr), start = as.numeric(dt$start),
                   end = as.numeric(dt$end), gene = as.character(dt$gene),
                   hp1_state = num(cn[["hp1"]]), hp2_state = num(cn[["hp2"]]))
  bad = is.na(out$hp1_state) | is.na(out$hp2_state) |
        out$hp1_state > max_cn | out$hp2_state > max_cn
  out = out[!bad]
  out[, `:=`(total = hp1_state + hp2_state, minor = pmin(hp1_state, hp2_state))]
  list(genes = out, n_dropped = sum(bad))
}

# SAVANA <id>_fitted_purity_ploidy.tsv (purity ploidy distance rank) -> the rank-1 fit, or NULL
parse_savana_purity = function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  dt = tryCatch(fread(path, sep = "\t", header = TRUE), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0 || !all(c("purity", "ploidy") %in% names(dt))) return(NULL)
  if ("rank" %in% names(dt)) setorder(dt, rank)
  list(purity = as.numeric(dt$purity[1]), ploidy = as.numeric(dt$ploidy[1]))
}

# A file inside one ranked Wakhan solution. Literal file.path() lookups, never a listing:
# solution_rank_<n> may itself be a symlinked directory (which recursive listing skips),
# and in the pipeline-staged layout only individual symlinked files exist. Tries
# solution_rank_<n> (0.5), solution_<n> (0.4.x), then repository_name (the real dir).
locate_wakhan_solution_file = function(wakhan_dir, rank, repository_name = NULL, relpath) {
  if (is.null(wakhan_dir) || is.na(wakhan_dir)) return(NULL)
  dirs = c(paste0("solution_rank_", rank), paste0("solution_", rank),
           if (!is.null(repository_name) && !is.na(repository_name)) repository_name)
  cand = file.path(wakhan_dir, dirs, relpath)
  hit = cand[file.exists(cand)]
  if (length(hit) > 0) hit[1] else NULL
}

