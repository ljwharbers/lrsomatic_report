# ShatterSeek chromothripsis section module.
#
# Input: shatterseek/<combo>/<id>.shatterseek.tsv, written by the pipeline's
# assets/shatterseek/run_shatterseek.R — one row per chromosome (chr1-22, chrX), the
# chromosome's largest cluster of interleaved SVs, ShatterSeek's statistics, BH q-values,
# the PCAWG high-confidence tiers HC1/HC2/HC3, `chromothripsis_criteria` ("HC1,HC2") and
# `chromothripsis_confidence` ("high" when any tier is met, else "none"). <combo> is
# severus_ascat, severus_wakhan or savana (sv_source_cn_source, or the one source when both
# are the same caller); every row also carries `cn_source`/`sv_source`, which win over the
# directory name.
#
# Contract
#   id "shatterseek"; parse(inputs, section_data, ctx) returns
#   list(status = "not_run")                         no ShatterSeek file found
#   list(
#     status = "ok",
#     combos = named list of data.table,  one per combo, names = combo ids ("severus_ascat"),
#                                          canonical columns, see SHATTERSEEK_COLUMNS
#     labels = named character,            combo id -> display label ("Severus + ASCAT")
#     high   = data.table(combo, combo_label, chrom, confidence, criteria, start, end,
#                         n_sv, n_sv_incl_tra, n_tra, n_oscillating)
#                                          one row per (combo, chromosome) with a
#                                          confidence other than "none"
#     files  = named character             combo id -> source TSV
#   )
#   Rendered by templates/sections/_instability.qmd.

# Canonical column -> accepted spellings (resolve_columns(): exact, then normalised match)
SHATTERSEEK_COLUMNS = list(
  sample            = c("sample", "sample_id"),
  cn_source         = c("cn_source"),
  sv_source         = c("sv_source"),
  chrom             = c("chrom", "chr", "chromosome"),
  start             = c("start"),
  end               = c("end"),
  n_sv              = c("clusterSize", "cluster_size"),
  n_del             = c("number_DEL"),
  n_dup             = c("number_DUP"),
  n_h2hinv          = c("number_h2hINV"),
  n_t2tinv          = c("number_t2tINV"),
  n_tra             = c("number_TRA"),
  n_sv_incl_tra     = c("clusterSize_including_TRA", "cluster_size_including_tra"),
  n_oscillating     = c("max_number_oscillating_CN_segments_2_states", "oscillating_cn_2_states"),
  p_fragment_joins  = c("pval_fragment_joins"),
  p_chr_enrichment  = c("chr_breakpoint_enrichment", "pval_chr_breakpoint_enrichment"),
  p_exp_chr         = c("pval_exp_chr"),
  p_exp_cluster     = c("pval_exp_cluster"),
  q_fragment_joins  = c("q_fragment_joins"),
  q_chr_enrichment  = c("q_chr_breakpoint_enrichment"),
  q_exp_chr         = c("q_exp_chr"),
  q_exp_cluster     = c("q_exp_cluster"),
  hc1               = c("HC1"),
  hc2               = c("HC2"),
  hc3               = c("HC3"),
  criteria          = c("chromothripsis_criteria", "criteria"),
  confidence        = c("chromothripsis_confidence", "confidence")
)

SHATTERSEEK_TIERS = c("hc1", "hc2", "hc3")

# Display name of one caller id; unknown ids pass through unchanged
shatterseek_tool_label = function(x) {
  known = c(severus = "Severus", ascat = "ASCAT", wakhan = "Wakhan", savana = "SAVANA")
  x = as.character(x)
  out = unname(known[tolower(x)])
  ifelse(is.na(out), x, out)
}

# "severus_ascat" -> "Severus + ASCAT"; "savana" -> "SAVANA"
shatterseek_combo_label = function(combo) {
  vapply(as.character(combo), function(cb) {
    parts = strsplit(cb, "_", fixed = TRUE)[[1]]
    paste(shatterseek_tool_label(parts), collapse = " + ")
  }, character(1), USE.NAMES = FALSE)
}

# TRUE/FALSE/"True"/"1" -> logical; anything unreadable is NA
.ss_logical = function(x) {
  if (is.logical(x)) return(x)
  v = tolower(trimws(as.character(x)))
  out = rep(NA, length(v))
  out[v %in% c("true", "t", "1", "yes")] = TRUE
  out[v %in% c("false", "f", "0", "no")] = FALSE
  out
}

# Natural chromosome rank (chr1..chr22, X, Y, then the rest)
.ss_chrom_rank = function(chrom) {
  s = sub("^chr", "", as.character(chrom), ignore.case = TRUE)
  r = suppressWarnings(as.integer(s))
  r[toupper(s) == "X"] = 23L
  r[toupper(s) == "Y"] = 24L
  r[is.na(r)] = 99L
  r
}

# Read one <id>.shatterseek.tsv into canonical columns. Returns NULL (with a message) when
# the file cannot be read or has no chromosome column; a 0-byte or header-only file returns
# an empty table. `dir_combo` is the directory name, used when the rows carry no source columns.
parse_shatterseek_file = function(path, dir_combo = basename(dirname(path))) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  if (isTRUE(file.size(path) == 0)) {
    return(list(combo = dir_combo, dt = data.table(chrom = character(0))))
  }
  raw = tryCatch(fread(path, sep = "\t", header = TRUE, na.strings = c("NA", ""),
                       colClasses = "character"),
                 error = function(e) {
                   message("ShatterSeek: cannot read ", path, ": ", conditionMessage(e))
                   NULL
                 })
  if (is.null(raw)) return(NULL)
  cols = resolve_columns(raw, SHATTERSEEK_COLUMNS)
  if (is.na(cols[["chrom"]])) {
    message("ShatterSeek: no chromosome column in ", path, "; skipped")
    return(NULL)
  }

  get_col = function(k) {
    if (is.na(cols[[k]])) rep(NA_character_, nrow(raw)) else as.character(raw[[cols[[k]]]])
  }
  num = function(k) suppressWarnings(as.numeric(get_col(k)))

  dt = data.table(chrom = ensure_chr_prefix(get_col("chrom")))
  for (k in c("start", "end", "n_sv", "n_del", "n_dup", "n_h2hinv", "n_t2tinv", "n_tra",
              "n_sv_incl_tra", "n_oscillating", "p_fragment_joins", "p_chr_enrichment",
              "p_exp_chr", "p_exp_cluster", "q_fragment_joins", "q_chr_enrichment",
              "q_exp_chr", "q_exp_cluster")) {
    set(dt, j = k, value = num(k))
  }
  for (k in SHATTERSEEK_TIERS) set(dt, j = k, value = .ss_logical(get_col(k)))

  # Criteria: as written, else rebuilt from the tier flags
  crit = trimws(get_col("criteria"))
  tier_cols = SHATTERSEEK_TIERS[!is.na(cols[SHATTERSEEK_TIERS])]
  if (is.na(cols[["criteria"]]) && length(tier_cols) > 0) {
    m = as.matrix(dt[, ..tier_cols])
    crit = apply(m, 1, function(r) paste(toupper(tier_cols)[!is.na(r) & r], collapse = ","))
  }
  crit[is.na(crit)] = ""
  dt[, criteria := crit]

  # Confidence: as written (lower-cased; blank/NA reads as "none"), else derived — any tier
  # met or any criteria listed means "high", which is what the wrapper writes
  conf = tolower(trimws(get_col("confidence")))
  if (is.na(cols[["confidence"]])) {
    if (length(tier_cols) > 0 || !is.na(cols[["criteria"]])) {
      any_tier = if (length(tier_cols) > 0)
        Reduce(`|`, lapply(tier_cols, function(k) !is.na(dt[[k]]) & dt[[k]])) else FALSE
      conf = ifelse(any_tier | nzchar(dt$criteria), "high", "none")
    } else {
      message("ShatterSeek: no confidence, criteria or tier columns in ", path)
    }
  } else {
    conf[is.na(conf) | conf %in% c("", "no", "false", "na")] = "none"
  }
  dt[, confidence := conf]

  dt[, sample := get_col("sample")]
  cn_src = unique(na.omit(get_col("cn_source")))
  sv_src = unique(na.omit(get_col("sv_source")))
  combo = if (length(cn_src) == 1 && length(sv_src) == 1) {
    if (tolower(cn_src) == tolower(sv_src)) tolower(cn_src) else tolower(paste(sv_src, cn_src, sep = "_"))
  } else dir_combo
  dt[, combo := combo]

  setorder(dt[, .rank := .ss_chrom_rank(chrom)], .rank)
  dt[, .rank := NULL]
  setcolorder(dt, c("combo", "chrom"))
  list(combo = combo, dt = dt[])
}

# Rows with a call ("high" or any other non-"none" tier), across combos
shatterseek_calls = function(combos, labels = NULL) {
  empty = data.table(combo = character(0), combo_label = character(0), chrom = character(0),
                     confidence = character(0), criteria = character(0),
                     start = numeric(0), end = numeric(0), n_sv = numeric(0),
                     n_sv_incl_tra = numeric(0), n_tra = numeric(0), n_oscillating = numeric(0))
  if (length(combos) == 0) return(empty)
  rows = rbindlist(lapply(names(combos), function(cb) {
    d = combos[[cb]]
    if (is.null(d) || nrow(d) == 0 || !"confidence" %in% names(d)) return(NULL)
    d = d[!is.na(confidence) & confidence != "none"]
    if (nrow(d) == 0) return(NULL)
    d[, .(combo = cb, chrom, confidence, criteria, start, end, n_sv, n_sv_incl_tra, n_tra,
          n_oscillating)]
  }), fill = TRUE)
  if (nrow(rows) == 0) return(empty)
  lab = if (!is.null(labels)) unname(labels[rows$combo]) else shatterseek_combo_label(rows$combo)
  rows[, combo_label := ifelse(is.na(lab), shatterseek_combo_label(combo), lab)]
  rows[, .rank := .ss_chrom_rank(chrom)]
  setorder(rows, .rank, combo)
  rows[, .rank := NULL]
  setcolorder(rows, c("combo", "combo_label", "chrom", "confidence"))
  rows[]
}

# One-line verdict over every combo. Mentions disagreement explicitly when more than one
# combo ran and they do not call the same chromosomes.
shatterseek_verdict = function(x) {
  if (is.null(x) || !identical(x$status, "ok")) return("ShatterSeek was not run for this sample.")
  labels = x$labels
  high = x$high
  if (is.null(high) || nrow(high) == 0) {
    # One line for every combination: identical negatives are not repeated per combination
    n_tested = vapply(names(labels), function(cb) {
      d = x$combos[[cb]]
      if (is.null(d)) NA_integer_ else nrow(d)
    }, integer(1))
    tested = if (length(labels) == 1)
               sprintf("%s; %s chromosomes tested", unname(labels), n_tested)
             else if (length(unique(n_tested)) == 1)
               sprintf("%s; %s chromosomes tested each", paste(unname(labels), collapse = ", "),
                       n_tested[1])
             else paste(sprintf("%s: %s chromosomes tested", unname(labels), n_tested), collapse = "; ")
    lead = if (length(labels) == 1) "No chromothripsis called"
           else sprintf("No chromothripsis called by any of the %d caller combinations", length(labels))
    return(sprintf("%s (%s).", lead, tested))
  }
  by_chrom = high[, .(combos = list(unique(combo)),
                      conf   = paste(unique(confidence), collapse = "/")),
                  by = chrom]
  all_high = all(high$confidence == "high")
  multi = length(labels) > 1
  pieces = vapply(seq_len(nrow(by_chrom)), function(i) {
    cbs = by_chrom$combos[[i]]
    tag = if (!multi) "" else if (length(cbs) == length(labels)) " (all combinations)"
          else sprintf(" (%s only)", paste(unname(labels[cbs]), collapse = ", "))
    conf = if (all_high) "" else sprintf(" [%s]", by_chrom$conf[i])
    paste0(by_chrom$chrom[i], conf, tag)
  }, character(1))
  lead = if (all_high) "High-confidence chromothripsis on " else "Chromothripsis called on "
  out = paste0(lead, paste(pieces, collapse = ", "))
  if (!multi) out = paste0(out, " (", unname(labels[1]), ")")
  paste0(out, ".")
}

# TRUE when more than one combo ran and their called chromosome sets differ
shatterseek_disagree = function(x) {
  if (is.null(x) || !identical(x$status, "ok") || length(x$combos) < 2) return(FALSE)
  sets = lapply(names(x$combos), function(cb) sort(unique(x$high[combo == cb, chrom])))
  length(unique(vapply(sets, paste, character(1), collapse = ","))) > 1
}

register_section(list(
  id     = "shatterseek",
  title  = "Chromothripsis (ShatterSeek)",

  locate = function(sample_dir, sample_id) {
    hits = list.files(sample_dir, pattern = "\\.shatterseek\\.tsv$", recursive = TRUE,
                      full.names = TRUE)
    # Prefer this sample's files when the tree holds others (a cohort directory)
    own = hits[startsWith(basename(hits), paste0(sample_id, "."))]
    list(files = sort(if (length(own) > 0) own else hits))
  },

  parse = function(inputs, section_data, ctx) {
    files = inputs$files
    if (length(files) == 0) return(list(status = "not_run"))

    combos = list(); src = character(0)
    for (f in files) {
      p = parse_shatterseek_file(f)
      if (is.null(p)) next
      if (!is.null(combos[[p$combo]])) {
        message("ShatterSeek: second file for combination '", p$combo, "' ignored: ", f)
        next
      }
      combos[[p$combo]] = p$dt
      src[[p$combo]] = f
    }
    if (length(combos) == 0) return(list(status = "not_run"))

    # Display order: the pipeline's pairs first, anything else alphabetically after
    pref = c("severus_ascat", "severus_wakhan", "savana")
    ord = order(match(names(combos), pref, nomatch = length(pref) + 1L), names(combos))
    combos = combos[ord]
    src = src[names(combos)]
    labels = setNames(shatterseek_combo_label(names(combos)), names(combos))

    list(
      status = "ok",
      combos = combos,
      labels = labels,
      high   = shatterseek_calls(combos, labels),
      files  = src
    )
  }
))
