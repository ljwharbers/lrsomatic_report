# ecDNA / focal amplicon section module (CoRAL + AmpliconClassifier).
#
# Inputs, all under ecdna/ (pipeline subworkflow ECDNA):
#   coral/<id>_CNV_SEEDS.bed                       amplified seed intervals; 0 bytes = none above --coral_gain
#   coral/{reconstruct,cycles}/<id>_amplicon_summary.txt
#                                                  CoRAL's per-amplicon summary (cycles/ wins: with
#                                                  --coral_run_cycle that is what the classifier reads)
#   coral/plots/<id>_amplicon<N>_{graph,cycles}.png
#   amplicon_classifier/<id>_amplicon_classification_profiles.tsv
#                       <id>_gene_list.tsv, <id>_result_table.tsv   (AmpliconClassifier 2.0.0 formats)
#   (ecdna/<id>_coral_cn.bed is CoRAL's input, not a result: ignored.)
#
# Contract
#   id "ecdna"; parse(inputs, section_data, ctx) returns
#   list(
#     status     = "not_run"        no ecDNA output at all
#                | "none_detected"  seeds empty, or classifier ran and found no valid amplicon
#                | "detected"       classifier found >= 1 amplicon (not necessarily ecDNA: see n_ecdna)
#                | "unclassified"   CoRAL found amplified seeds / reconstructed amplicons, but no
#                                   AmpliconClassifier output exists (e.g. CHM13 run without
#                                   --aa_data_repo, or --skip_ampliconclassifier)
#     amplicons  = data.table(amplicon, class, ecdna, bfb, location, oncogenes, max_cn,
#                             features, size_bp, source)
#                    amplicon int; ecdna/bfb logical (NA when unclassified); location/oncogenes
#                    comma-separated character; source "AmpliconClassifier" | "CoRAL"
#     genes      = data.table(amplicon, feature, gene, gene_cn, oncogene, truncated)
#     plots      = character()   PNG paths, amplicon order, graph before cycles
#     n_ecdna    = integer, amplicons classified ecDNA+ (NA when unclassified)
#     n_bfb      = integer, amplicons classified BFB+   (NA when unclassified)
#     n_invalid  = integer, classifier rows dropped as "No amp/Invalid"
#     files      = named list of the source paths used
#   )
#   Rendered by templates/sections/_instability.qmd.

ECDNA_PROFILE_COLUMNS = list(
  sample    = c("sample_name", "sample"),
  amplicon  = c("amplicon_number", "amplicon", "AA amplicon number"),
  class     = c("amplicon_decomposition_class", "classification", "class"),
  ecdna     = c("ecDNA+", "ecDNA_positive", "ecdna"),
  bfb       = c("BFB+", "BFB_positive", "bfb"),
  n_species = c("ecDNA_amplicons", "ecDNA amplicons")
)

ECDNA_GENE_COLUMNS = list(
  amplicon  = c("amplicon_number", "amplicon", "AA amplicon number"),
  feature   = c("feature", "Feature ID", "feature_id"),
  gene      = c("gene", "gene_name", "symbol"),
  gene_cn   = c("gene_cn", "gene_copy_number", "copy_number"),
  truncated = c("truncated"),
  oncogene  = c("is_canonical_oncogene", "is_oncogene", "oncogene")
)

ECDNA_RESULT_COLUMNS = list(
  amplicon   = c("AA amplicon number", "amplicon_number", "amplicon"),
  feature    = c("Feature ID", "feature_id", "feature"),
  class      = c("Classification", "classification"),
  location   = c("Location", "location"),
  oncogenes  = c("Oncogenes", "oncogenes"),
  max_cn     = c("Feature maximum copy number", "feature_maximum_copy_number", "max_cn"),
  median_cn  = c("Feature median copy number", "feature_median_copy_number"),
  length     = c("Captured interval length", "captured_interval_length")
)

# "amplicon3" / "3" / "sample_amplicon3" -> 3L
.ec_amplicon_id = function(x) {
  suppressWarnings(as.integer(sub(".*?([0-9]+)$", "\\1", trimws(as.character(x)), perl = TRUE)))
}

# "Positive" / "TRUE" / "yes" -> TRUE; "None detected" / "Negative" / "FALSE" -> FALSE
.ec_flag = function(x) {
  v = tolower(trimws(as.character(x)))
  out = rep(NA, length(v))
  out[grepl("^(pos|true|t$|yes|y$|1$)", v)] = TRUE
  out[grepl("^(none|neg|false|f$|no|n$|0$)", v)] = FALSE
  out
}

# Python-list cells ("['chr8:1-2', 'chr8:5-9']", "['MYC']") or "|"/";"-joined -> "a, b"
.ec_list_cell = function(x) {
  vapply(as.character(x), function(s) {
    if (is.na(s)) return(NA_character_)
    s = gsub("[][\"']", "", s)
    parts = trimws(unlist(strsplit(s, "[,|;]")))
    parts = unique(parts[nzchar(parts) & !tolower(parts) %in% c("none", "nan", "na")])
    if (length(parts) == 0) NA_character_ else paste(parts, collapse = ", ")
  }, character(1), USE.NAMES = FALSE)
}

# Read a tab-separated file as all-character; NULL for absent/0-byte/unreadable files.
# A header-only file returns a zero-row table, which callers treat as "no rows".
.ec_read = function(path) {
  if (is.null(path) || !file.exists(path) || isTRUE(file.size(path) == 0)) return(NULL)
  tryCatch(fread(path, sep = "\t", header = TRUE, colClasses = "character",
                 na.strings = c("NA", "", "nan", "None")),
           error = function(e) {
             message("ecDNA: cannot read ", path, ": ", conditionMessage(e))
             NULL
           })
}

.ec_col = function(dt, cols, k) {
  if (is.na(cols[[k]])) rep(NA_character_, nrow(dt)) else as.character(dt[[cols[[k]]]])
}

parse_ac_profiles = function(path) {
  raw = .ec_read(path)
  if (is.null(raw) || nrow(raw) == 0) return(NULL)
  cols = resolve_columns(raw, ECDNA_PROFILE_COLUMNS)
  if (is.na(cols[["amplicon"]])) {
    message("ecDNA: no amplicon column in ", path, "; skipped")
    return(NULL)
  }
  data.table(
    amplicon  = .ec_amplicon_id(.ec_col(raw, cols, "amplicon")),
    class     = .ec_col(raw, cols, "class"),
    ecdna     = .ec_flag(.ec_col(raw, cols, "ecdna")),
    bfb       = .ec_flag(.ec_col(raw, cols, "bfb")),
    n_species = suppressWarnings(as.integer(.ec_col(raw, cols, "n_species")))
  )
}

parse_ac_gene_list = function(path) {
  raw = .ec_read(path)
  if (is.null(raw) || nrow(raw) == 0) return(NULL)
  cols = resolve_columns(raw, ECDNA_GENE_COLUMNS)
  if (is.na(cols[["gene"]])) {
    message("ecDNA: no gene column in ", path, "; skipped")
    return(NULL)
  }
  data.table(
    amplicon  = .ec_amplicon_id(.ec_col(raw, cols, "amplicon")),
    feature   = .ec_col(raw, cols, "feature"),
    gene      = .ec_col(raw, cols, "gene"),
    gene_cn   = suppressWarnings(as.numeric(.ec_col(raw, cols, "gene_cn"))),
    oncogene  = .ec_flag(.ec_col(raw, cols, "oncogene")),
    truncated = .ec_col(raw, cols, "truncated")
  )
}

parse_ac_result_table = function(path) {
  raw = .ec_read(path)
  if (is.null(raw) || nrow(raw) == 0) return(NULL)
  cols = resolve_columns(raw, ECDNA_RESULT_COLUMNS)
  if (is.na(cols[["amplicon"]])) {
    message("ecDNA: no amplicon column in ", path, "; skipped")
    return(NULL)
  }
  data.table(
    amplicon  = .ec_amplicon_id(.ec_col(raw, cols, "amplicon")),
    feature   = .ec_col(raw, cols, "feature"),
    class     = .ec_col(raw, cols, "class"),
    location  = .ec_list_cell(.ec_col(raw, cols, "location")),
    oncogenes = .ec_list_cell(.ec_col(raw, cols, "oncogenes")),
    max_cn    = suppressWarnings(as.numeric(.ec_col(raw, cols, "max_cn"))),
    length    = suppressWarnings(as.numeric(.ec_col(raw, cols, "length")))
  )
}

# CoRAL <id>_amplicon_summary.txt -> data.table(amplicon, location, size_bp, max_cn, heaviest_walk).
# Blocks start at "AmpliconID = N"; intervals read "Amplicon N > chr15:17,792,835-19,305,544";
# max_cn is the largest Copy_count of a decomposed cycle/path (a walk's copy count, not a
# segment copy number).
parse_coral_summary = function(path) {
  if (is.null(path) || !file.exists(path) || isTRUE(file.size(path) == 0)) return(NULL)
  lines = readLines(path, warn = FALSE)
  starts = grep("^\\s*AmpliconID\\s*=", lines)
  if (length(starts) == 0) return(NULL)
  ends = c(starts[-1] - 1L, length(lines))
  rbindlist(lapply(seq_along(starts), function(i) {
    blk = lines[starts[i]:ends[i]]
    id = .ec_amplicon_id(sub("^\\s*AmpliconID\\s*=\\s*", "", blk[1]))
    iv = regmatches(blk, regexpr("[^\\s>]+:[0-9,]+-[0-9,]+", blk, perl = TRUE))
    iv = gsub(",", "", iv)
    size = suppressWarnings(as.numeric(gsub("[^0-9]", "",
              sub(".*:", "", grep("Total Amplicon Size", blk, value = TRUE)[1]))))
    cc = suppressWarnings(as.numeric(sub(".*Copy_count=([0-9.eE+-]+).*", "\\1",
                                         grep("Copy_count=", blk, value = TRUE))))
    walk = grep("Heaviest graph walk solved was a", blk, value = TRUE)[1]
    data.table(
      amplicon      = id,
      location      = if (length(iv)) paste(unique(iv), collapse = ", ") else NA_character_,
      size_bp       = size,
      max_cn        = if (length(cc) && any(!is.na(cc))) max(cc, na.rm = TRUE) else NA_real_,
      heaviest_walk = if (is.na(walk)) NA_character_ else sub("\\.$", "", sub(".*solved was a ", "", walk))
    )
  }))
}

.ec_empty_amplicons = function() {
  data.table(amplicon = integer(0), class = character(0), ecdna = logical(0), bfb = logical(0),
             location = character(0), oncogenes = character(0), max_cn = numeric(0),
             features = character(0), size_bp = numeric(0), source = character(0))
}

.ec_empty_genes = function() {
  data.table(amplicon = integer(0), feature = character(0), gene = character(0),
             gene_cn = numeric(0), oncogene = logical(0), truncated = character(0))
}

# Assemble the per-amplicon table from whichever classifier files exist; CoRAL's summary
# fills location/size/max CN where the classifier has none.
build_ecdna_amplicons = function(profiles, results, genes, coral) {
  first_non_na = function(...) {
    for (v in list(...)) if (length(v) && !all(is.na(v))) return(v)
    list(...)[[1]]
  }
  per_amp_results = NULL
  if (!is.null(results) && nrow(results) > 0) {
    per_amp_results = results[!is.na(amplicon), .(
      rt_class  = paste(unique(na.omit(class)), collapse = ", "),
      rt_ecdna  = any(grepl("ecdna", tolower(class))),
      rt_bfb    = any(grepl("bfb", tolower(class))),
      location  = { v = unique(na.omit(unlist(strsplit(na.omit(location), ", ")))); if (length(v)) paste(v, collapse = ", ") else NA_character_ },
      oncogenes = { v = unique(na.omit(unlist(strsplit(na.omit(oncogenes), ", ")))); if (length(v)) paste(v, collapse = ", ") else NA_character_ },
      max_cn    = if (all(is.na(max_cn))) NA_real_ else max(max_cn, na.rm = TRUE),
      features  = paste(unique(na.omit(feature)), collapse = ", "),
      size_bp   = if (all(is.na(length))) NA_real_ else sum(length, na.rm = TRUE)
    ), by = amplicon]
  }

  if (!is.null(profiles) && nrow(profiles) > 0) {
    amps = profiles[!is.na(amplicon), .(amplicon, class, ecdna, bfb)]
  } else if (!is.null(per_amp_results)) {
    amps = per_amp_results[, .(amplicon, class = rt_class, ecdna = rt_ecdna, bfb = rt_bfb)]
  } else {
    return(NULL)
  }
  amps = unique(amps, by = "amplicon")
  # Flag columns absent or blank: read them off the class name instead
  amps[is.na(ecdna) & !is.na(class), ecdna := grepl("ecdna", tolower(class))]
  amps[is.na(bfb) & !is.na(class), bfb := grepl("bfb", tolower(class))]
  if (!is.null(per_amp_results)) {
    amps = merge(amps, per_amp_results[, .(amplicon, location, oncogenes, max_cn, features, size_bp)],
                 by = "amplicon", all.x = TRUE)
  } else {
    amps[, `:=`(location = NA_character_, oncogenes = NA_character_, max_cn = NA_real_,
                features = NA_character_, size_bp = NA_real_)]
  }

  # Gene list: oncogenes and max gene CN where the result table gave none
  if (!is.null(genes) && nrow(genes) > 0) {
    g = genes[!is.na(amplicon), .(
      g_onc  = { v = unique(gene[!is.na(oncogene) & oncogene]); if (length(v)) paste(v, collapse = ", ") else NA_character_ },
      g_cn   = if (all(is.na(gene_cn))) NA_real_ else max(gene_cn, na.rm = TRUE),
      g_feat = paste(unique(na.omit(feature)), collapse = ", ")
    ), by = amplicon]
    amps = merge(amps, g, by = "amplicon", all.x = TRUE)
    amps[is.na(oncogenes), oncogenes := g_onc]
    amps[is.na(max_cn), max_cn := g_cn]
    amps[is.na(features) | !nzchar(features), features := g_feat]
    amps[, c("g_onc", "g_cn", "g_feat") := NULL]
  }
  if (!is.null(coral) && nrow(coral) > 0) {
    amps = merge(amps, coral[, .(amplicon, c_loc = location, c_size = size_bp, c_cn = max_cn)],
                 by = "amplicon", all.x = TRUE)
    amps[is.na(location), location := c_loc]
    amps[is.na(size_bp), size_bp := c_size]
    amps[is.na(max_cn), max_cn := c_cn]
    amps[, c("c_loc", "c_size", "c_cn") := NULL]
  }
  amps[features == "", features := NA_character_]
  amps[, source := "AmpliconClassifier"]
  setorder(amps, amplicon)
  setcolorder(amps, names(.ec_empty_amplicons()))
  amps[]
}

# One-line verdict for the section and the overview
ecdna_verdict = function(x) {
  if (is.null(x) || identical(x$status, "not_run"))
    return("ecDNA analysis (CoRAL / AmpliconClassifier) was not run for this sample.")
  n = if (is.null(x$amplicons)) 0L else nrow(x$amplicons)
  switch(x$status,
    none_detected = if (isTRUE(x$seeds_empty))
        "No focal amplification above the CoRAL seed threshold: no amplicons to reconstruct."
      else "No focal amplicons classified by AmpliconClassifier.",
    unclassified = if (n > 0)
        sprintf("%d amplicon%s reconstructed by CoRAL; not classified (no AmpliconClassifier output), so ecDNA status is unknown.",
                n, if (n == 1) "" else "s")
      else "Amplified seed regions found, but no CoRAL reconstruction or AmpliconClassifier output: ecDNA status is unknown.",
    detected = {
      bits = c(if (isTRUE(x$n_ecdna > 0)) sprintf("%d ecDNA+", x$n_ecdna),
               if (isTRUE(x$n_bfb > 0)) sprintf("%d BFB+", x$n_bfb))
      if (length(bits) == 0)
        sprintf("%d amplicon%s found, none classified as ecDNA or BFB.", n, if (n == 1) "" else "s")
      else sprintf("%s of %d amplicon%s.", paste(bits, collapse = " and "), n, if (n == 1) "" else "s")
    },
    "ecDNA status unknown."
  )
}

# Footnote naming the folder(s) of the CoRAL plots past the embed cap, relative to --sample-dir
ecdna_more_plots_note = function(plots, cap, sample_dir) {
  rest = utils::tail(plots, max(0L, length(plots) - cap))
  if (length(rest) == 0) return(NULL)
  htmltools::tags$p(class = "table-footnote",
    sprintf("%d more CoRAL plot%s not embedded (--embed slim); see ", length(rest),
            if (length(rest) == 1) "" else "s"),
    htmltools::tags$code(paste(unique(sample_rel_path(dirname(rest), sample_dir)), collapse = ", ")), ".")
}

register_section(list(
  id     = "ecdna",
  title  = "ecDNA",

  locate = function(sample_dir, sample_id) {
    pat = paste0("(_CNV_SEEDS\\.bed|_amplicon_classification_profiles\\.tsv|_gene_list\\.tsv|",
                 "_result_table\\.tsv|_amplicon_summary\\.txt|_amplicon[0-9]+_(graph|cycles)\\.png)$")
    hits = list.files(sample_dir, pattern = pat, recursive = TRUE, full.names = TRUE)
    # _gene_list.tsv / _result_table.tsv are generic suffixes: only under ecdna/ or amplicon_classifier/
    rel = substring(hits, nchar(sample_dir) + 1L)
    generic = grepl("(_gene_list|_result_table)\\.tsv$", hits)
    hits = hits[!generic | grepl("(^|/)(ecdna|amplicon_classifier)(/|$)", rel)]
    own = hits[startsWith(basename(hits), paste0(sample_id, "_"))]
    if (length(own) > 0) hits = own

    pick = function(re, prefer = NULL) {
      h = sort(hits[grepl(re, hits)])
      if (length(h) == 0) return(NULL)
      if (!is.null(prefer)) {
        p = h[grepl(prefer, h)]
        if (length(p) > 0) return(p[1])
      }
      h[1]
    }
    plots = hits[grepl("_amplicon[0-9]+_(graph|cycles)\\.png$", hits)]
    if (length(plots) > 0) {
      n = .ec_amplicon_id(sub("_(graph|cycles)\\.png$", "", basename(plots)))
      kind = ifelse(grepl("_graph\\.png$", plots), 1L, 2L)
      plots = plots[order(n, kind)]
    }
    list(
      seeds    = pick("_CNV_SEEDS\\.bed$"),
      profiles = pick("_amplicon_classification_profiles\\.tsv$"),
      genes    = pick("_gene_list\\.tsv$"),
      results  = pick("_result_table\\.tsv$"),
      summary  = pick("_amplicon_summary\\.txt$", prefer = "/cycles/"),
      plots    = plots
    )
  },

  parse = function(inputs, section_data, ctx) {
    files = Filter(Negate(is.null), inputs[c("seeds", "profiles", "genes", "results", "summary")])
    plots = if (is.null(inputs$plots)) character(0) else inputs$plots
    base = list(amplicons = .ec_empty_amplicons(), genes = .ec_empty_genes(), plots = plots,
                n_ecdna = NA_integer_, n_bfb = NA_integer_, n_invalid = 0L, files = files,
                seeds_empty = FALSE)
    if (length(files) == 0 && length(plots) == 0) return(c(list(status = "not_run"), base))

    seeds_empty = !is.null(inputs$seeds) && isTRUE(file.size(inputs$seeds) == 0)
    base$seeds_empty = seeds_empty
    if (seeds_empty) return(c(list(status = "none_detected"), base))

    coral    = parse_coral_summary(inputs$summary)
    have_ac  = !is.null(inputs$profiles) || !is.null(inputs$results)
    profiles = parse_ac_profiles(inputs$profiles)
    results  = parse_ac_result_table(inputs$results)
    genes    = parse_ac_gene_list(inputs$genes)

    if (!have_ac) {
      # CoRAL ran but nothing classified its amplicons
      if (!is.null(coral) && nrow(coral) > 0) {
        base$amplicons = coral[, .(amplicon, class = NA_character_, ecdna = NA, bfb = NA,
                                   location, oncogenes = NA_character_, max_cn,
                                   features = NA_character_, size_bp, source = "CoRAL")]
        setorder(base$amplicons, amplicon)
      }
      return(c(list(status = "unclassified"), base))
    }

    amps = build_ecdna_amplicons(profiles, results, genes, coral)
    if (!is.null(amps) && nrow(amps) > 0) {
      invalid = !is.na(amps$class) & grepl("no amp|invalid", tolower(amps$class))
      base$n_invalid = sum(invalid)
      amps = amps[!invalid]
    }
    if (!is.null(genes) && nrow(genes) > 0) {
      base$genes = genes[is.na(amplicon) | amplicon %in% amps$amplicon]
    }
    if (is.null(amps) || nrow(amps) == 0) return(c(list(status = "none_detected"), base))

    base$amplicons = amps
    base$n_ecdna = sum(amps$ecdna %in% TRUE)
    base$n_bfb   = sum(amps$bfb %in% TRUE)
    c(list(status = "detected"), base)
  }
))
