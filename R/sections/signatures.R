# Mutational signatures section module: SigProfilerAssignment (COSMIC) fits + MatrixGenerator profiles
#
# Contract
#   id "signatures"; parse(inputs, section_data, ctx) returns NULL or
#   list(
#     cosmic_version,                       # "3.6", or NA when no COSMIC_v<ver> dir names it
#     lookup_version,                       # version of the bundled aetiology TSV actually used
#     contexts = list(
#       SBS = list(activities = data.table(signature, mutations, pct, aetiology, artefact, note,
#                                          in_lookup),
#                  stats      = list(total, cosine, l1_pct, l2_pct, kl, correlation),
#                  matrix     = data.table(channel, group, count),
#                  flags      = character(),       # subset of SIG_FLAG_LABELS names
#                  files      = list(activities, stats, matrix)),
#       DBS = <same shape>,
#       ID  = <same shape>),                 # a context is absent when its fit dir is absent
#     dominant = data.table(context, signature, mutations, pct, aetiology, artefact)
#   )
#   Rendered by templates/sections/_signatures.qmd.
#
# `stats` fields are NA when Assignment_Solution_Samples_Stats.txt is absent (older pipeline
# versions did not stage it); `total` then falls back to the matrix sum. `activities` keeps only
# signatures with > 0 mutations, sorted by pct desc; zero rows means "no signature assigned".

# ── Thresholds (one block; the footnote prints them) ─────────────────────────────────────────
SIG_MIN_COSINE       = 0.9
SIG_MIN_MUTATIONS    = c(SBS = 200, DBS = 50, ID = 100)
SIG_DOMINANT_MIN_PCT = 5

SIG_CONTEXTS     = c(SBS = "SBS96", DBS = "DBS78", ID = "ID83")
SIG_FLAG_LABELS  = c(few_mutations = "Few mutations", poor_fit = "Poor fit")

# SBS 6-class colours: the values of SNV_COLOURS in R/circos.R (not sourced here: needs circlize)
SBS6_COLOURS = c(
  "C>A" = "#2EBAED",
  "C>G" = "#1b1e22",
  "C>T" = "#b3402f",
  "T>A" = "#c7c2b8",
  "T>C" = "#ADCC54",
  "T>G" = "#F0D0CE"
)
# SigProfilerPlotting's DBS78 / ID83 class colours, so the profiles look like COSMIC's
DBS10_COLOURS = c(
  "AC>NN" = "#03BCEE", "AT>NN" = "#0366CB", "CC>NN" = "#A1CE63", "CG>NN" = "#036601",
  "CT>NN" = "#FE9898", "GC>NN" = "#E32926", "TA>NN" = "#FEB065", "TC>NN" = "#FE8002",
  "TG>NN" = "#CB98FE", "TT>NN" = "#4C0198"
)
ID16_COLOURS = c(
  "1bp Del C" = "#FBBD75", "1bp Del T" = "#FC8002", "1bp Ins C" = "#B0DB8E", "1bp Ins T" = "#3CA02F",
  "2bp Del repeat" = "#FDCAB5", "3bp Del repeat" = "#FC8A6A", "4bp Del repeat" = "#F14432",
  "5+bp Del repeat" = "#BC1A1A",
  "2bp Ins repeat" = "#D0E1F2", "3bp Ins repeat" = "#94C4DF", "4bp Ins repeat" = "#4A98C9",
  "5+bp Ins repeat" = "#1764AB",
  "2bp Del MH" = "#E1E1EF", "3bp Del MH" = "#B6B6D8", "4bp Del MH" = "#8683BD",
  "5+bp Del MH" = "#62409B"
)

# ── Standard COSMIC channel order ────────────────────────────────────────────────────────────
# MatrixGenerator writes SBS96 sorted alphabetically (A[C>A]A, A[C>A]C, ... by 5' base first),
# which is NOT the plotting order (substitution class first); always reorder through these.
sig_channels = function(context) {
  context = sig_context_key(context)
  b = c("A", "C", "G", "T")
  if (context == "SBS") {
    g = expand.grid(three = b, five = b, sub = names(SBS6_COLOURS), stringsAsFactors = FALSE)
    return(paste0(g$five, "[", g$sub, "]", g$three))
  }
  if (context == "DBS") {
    alts = list(
      AC = c("CA", "CG", "CT", "GA", "GG", "GT", "TA", "TG", "TT"),
      AT = c("CA", "CC", "CG", "GA", "GC", "TA"),
      CC = c("AA", "AG", "AT", "GA", "GG", "GT", "TA", "TG", "TT"),
      CG = c("AT", "GC", "GT", "TA", "TC", "TT"),
      CT = c("AA", "AC", "AG", "GA", "GC", "GG", "TA", "TC", "TG"),
      GC = c("AA", "AG", "AT", "CA", "CG", "TA"),
      TA = c("AT", "CG", "CT", "GC", "GG", "GT"),
      TC = c("AA", "AG", "AT", "CA", "CG", "CT", "GA", "GG", "GT"),
      TG = c("AA", "AC", "AT", "CA", "CC", "CT", "GA", "GC", "GT"),
      TT = c("AA", "AC", "AG", "CA", "CC", "CG", "GA", "GC", "GG")
    )
    return(unlist(lapply(names(alts), function(r) paste0(r, ">", alts[[r]])), use.names = FALSE))
  }
  # ID83
  one  = unlist(lapply(c("1:Del:C", "1:Del:T", "1:Ins:C", "1:Ins:T"),
                       function(p) paste0(p, ":", 0:5)))
  rep_ = unlist(lapply(c("Del", "Ins"), function(t)
    unlist(lapply(2:5, function(n) paste0(n, ":", t, ":R:", 0:5)))))
  mh   = unlist(lapply(2:5, function(n) paste0(n, ":Del:M:", seq_len(if (n == 5) 5 else n - 1))))
  c(one, rep_, mh)
}

# "SBS96"/"sbs"/"SBS" -> "SBS" (likewise DBS, ID)
sig_context_key = function(context) {
  k = toupper(sub("[0-9]+$", "", context))
  if (!k %in% names(SIG_CONTEXTS)) stop("Unknown signature context: ", context)
  k
}

# Channel -> plotting class: SBS 6 substitution classes, DBS 10 reference dinucleotides,
# ID 16 standard ID83 classes. NA for a channel that does not parse.
sig_channel_group = function(channel, context) {
  context = sig_context_key(context)
  channel = as.character(channel)
  if (context == "SBS") {
    g = sub("^.\\[(.>.)\\].$", "\\1", channel)
    g[g == channel] = NA_character_
    return(g)
  }
  if (context == "DBS") {
    ok = grepl("^[ACGT]{2}>[ACGT]{2}$", channel)
    return(ifelse(ok, paste0(substr(channel, 1, 2), ">NN"), NA_character_))
  }
  parts = strsplit(channel, ":", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) != 4) return(NA_character_)
    len = suppressWarnings(as.integer(p[1]))
    if (is.na(len) || !p[2] %in% c("Del", "Ins")) return(NA_character_)
    if (len == 1 && p[3] %in% c("C", "T")) return(paste("1bp", p[2], p[3]))
    lab = if (len >= 5) "5+bp" else paste0(len, "bp")
    if (len >= 2 && p[3] == "R") return(paste(lab, p[2], "repeat"))
    if (len >= 2 && p[3] == "M" && p[2] == "Del") return(paste(lab, "Del MH"))
    NA_character_
  }, character(1), USE.NAMES = FALSE)
}

sig_group_levels = function(context) {
  switch(sig_context_key(context),
         SBS = names(SBS6_COLOURS), DBS = names(DBS10_COLOURS), ID = names(ID16_COLOURS))
}

sig_group_colours = function(context) {
  switch(sig_context_key(context), SBS = SBS6_COLOURS, DBS = DBS10_COLOURS, ID = ID16_COLOURS)
}

# ── Parsers ──────────────────────────────────────────────────────────────────────────────────
# Assignment_Solution_Activities.txt: "Samples\tSBS1\tSBS2..." one row per sample, mostly zeros.
# Returns a long data.table(signature, mutations, pct) with zeros dropped, pct desc; NULL when
# the file is absent or unreadable; a zero-row table when nothing was assigned.
parse_spa_activities = function(f, sample_id = NULL) {
  if (is.null(f) || !file.exists(f)) return(NULL)
  dt = tryCatch(fread(f, sep = "\t", header = TRUE, colClasses = list(character = 1)),
                error = function(e) { message("Failed to read ", f, ": ", conditionMessage(e)); NULL })
  if (is.null(dt) || nrow(dt) == 0 || ncol(dt) < 2) return(NULL)
  # One sample per run in this pipeline; prefer the report's sample if several are present
  row = 1L
  if (!is.null(sample_id) && nrow(dt) > 1) {
    hit = which(dt[[1]] %in% c(sample_id, gsub("[^A-Za-z0-9_-]", "_", sample_id)))
    if (length(hit) > 0) row = hit[1]
  }
  sigs = names(dt)[-1]
  vals = suppressWarnings(as.numeric(unlist(dt[row, -1, with = FALSE], use.names = FALSE)))
  out  = data.table(signature = sigs, mutations = vals)
  out  = out[!is.na(mutations) & mutations > 0]
  tot  = sum(out$mutations)
  out[, pct := if (tot > 0) 100 * mutations / tot else numeric(0)]
  setorder(out, -pct)
  out[]
}

# Assignment_Solution_Samples_Stats.txt -> list(total, cosine, l1_pct, l2_pct, kl, correlation);
# NULL when absent. Percent columns carry a literal "%", stripped here.
parse_spa_stats = function(f, sample_id = NULL) {
  if (is.null(f) || !file.exists(f)) return(NULL)
  dt = tryCatch(fread(f, sep = "\t", header = TRUE, colClasses = "character"),
                error = function(e) { message("Failed to read ", f, ": ", conditionMessage(e)); NULL })
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  row = 1L
  if (!is.null(sample_id) && nrow(dt) > 1) {
    hit = which(dt[[1]] %in% c(sample_id, gsub("[^A-Za-z0-9_-]", "_", sample_id)))
    if (length(hit) > 0) row = hit[1]
  }
  cols = resolve_columns(dt, list(
    total = "Total Mutations", cosine = "Cosine Similarity", l1_pct = "L1_Norm_%",
    l2_pct = "L2_Norm_%", kl = "KL Divergence", correlation = "Correlation"))
  num = function(key) {
    col = cols[[key]]
    if (is.na(col)) return(NA_real_)
    suppressWarnings(as.numeric(sub("%\\s*$", "", trimws(dt[[col]][row]))))
  }
  list(total = num("total"), cosine = num("cosine"), l1_pct = num("l1_pct"),
       l2_pct = num("l2_pct"), kl = num("kl"), correlation = num("correlation"))
}

# MatrixGenerator <name>.{SBS96,DBS78,ID83}.all ("MutationType\t<sample>") ->
# data.table(channel, group, count) in standard COSMIC order, every channel present (absent
# ones count 0). Channels outside the standard set are dropped with a message.
parse_sig_matrix = function(f, context, sample_id = NULL) {
  if (is.null(f) || !file.exists(f)) return(NULL)
  dt = tryCatch(fread(f, sep = "\t", header = TRUE),
                error = function(e) { message("Failed to read ", f, ": ", conditionMessage(e)); NULL })
  if (is.null(dt) || ncol(dt) < 2) return(NULL)
  col = 2L
  if (!is.null(sample_id) && ncol(dt) > 2) {
    hit = match(c(sample_id, gsub("[^A-Za-z0-9_-]", "_", sample_id)), names(dt))
    hit = hit[!is.na(hit)]
    if (length(hit) > 0) col = hit[1]
  }
  obs = data.table(channel = as.character(dt[[1]]),
                   count = suppressWarnings(as.numeric(dt[[col]])))
  std = sig_channels(context)
  extra = setdiff(obs$channel, std)
  if (length(extra) > 0)
    message(sprintf("%s matrix: %d channel(s) outside the standard set dropped", context, length(extra)))
  out = data.table(channel = std)
  out[, count := obs$count[match(channel, obs$channel)]]
  out[is.na(count), count := 0]
  out[, group := sig_channel_group(channel, context)]
  setcolorder(out, c("channel", "group", "count"))
  out[]
}

# Bundled aetiology lookup assets/signatures/cosmic_v<version>_aetiology.tsv. Falls back to
# the highest bundled version (with a message) when this version has no file of its own;
# returns NULL when there is none at all. The chosen version is attr(, "version").
load_sig_aetiology = function(dir, version = NA) {
  files = list.files(dir, pattern = "^cosmic_v[0-9.]+_aetiology\\.tsv$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  vers = sub("^cosmic_v([0-9.]+)_aetiology\\.tsv$", "\\1", basename(files))
  pick = match(as.character(version), vers)
  if (is.na(pick)) {
    pick = order(numeric_version(vers), decreasing = TRUE)[1]
    if (!is.na(version))
      message(sprintf("No aetiology lookup for COSMIC v%s; using v%s", version, vers[pick]))
  }
  lk = fread(files[pick], sep = "\t", header = TRUE, colClasses = "character")
  lk[, artefact := toupper(artefact) %in% c("TRUE", "T", "1", "YES")]
  lk[is.na(note), note := ""]
  setattr(lk, "version", vers[pick])
  lk[]
}

# Attach aetiology/artefact/note to an activities table; a signature absent from the lookup
# reads "not in lookup", never an error
annotate_sig_activities = function(act, lookup) {
  if (is.null(act)) return(NULL)
  act = copy(act)
  i = if (is.null(lookup)) rep(NA_integer_, nrow(act)) else match(act$signature, lookup$signature)
  act[, in_lookup := !is.na(i)]
  act[, aetiology := if (is.null(lookup)) rep("not in lookup", .N)
                     else ifelse(is.na(i), "not in lookup", lookup$aetiology[i])]
  act[, artefact := if (is.null(lookup)) rep(FALSE, .N)
                    else ifelse(is.na(i), FALSE, lookup$artefact[i])]
  act[, note := if (is.null(lookup)) rep("", .N) else ifelse(is.na(i), "", lookup$note[i])]
  setcolorder(act, c("signature", "mutations", "pct", "aetiology", "artefact", "note", "in_lookup"))
  act[]
}

# Flags for one context: few_mutations (total below SIG_MIN_MUTATIONS) and poor_fit (cosine
# below SIG_MIN_COSINE). An unknown cosine (no stats file) raises no poor-fit flag.
sig_context_flags = function(context, total, cosine) {
  context = sig_context_key(context)
  f = character(0)
  if (!is.na(total) && total < SIG_MIN_MUTATIONS[[context]]) f = c(f, "few_mutations")
  if (!is.na(cosine) && cosine < SIG_MIN_COSINE) f = c(f, "poor_fit")
  f
}

# Signatures at >= SIG_DOMINANT_MIN_PCT in contexts without a flag; context order, pct desc
sig_dominant = function(contexts) {
  empty = data.table(context = character(), signature = character(), mutations = numeric(),
                     pct = numeric(), aetiology = character(), artefact = logical())
  rows = lapply(intersect(names(SIG_CONTEXTS), names(contexts)), function(k) {
    cx = contexts[[k]]
    if (length(cx$flags) > 0 || is.null(cx$activities) || nrow(cx$activities) == 0) return(NULL)
    a = cx$activities[pct >= SIG_DOMINANT_MIN_PCT]
    if (nrow(a) == 0) return(NULL)
    a[, .(context = k, signature, mutations, pct, aetiology, artefact)][order(-pct)]
  })
  out = rbindlist(c(list(empty), rows), use.names = TRUE)
  out[]
}

# ── Plot ─────────────────────────────────────────────────────────────────────────────────────
# Short strip labels for the ID83 facets (the full class names do not fit a 1-6 channel strip)
.sig_strip_label = function(group, context) {
  if (context != "ID") return(group)
  g = sub("^1bp ", "", group)
  g = sub("bp Del repeat$", " Del", g)
  g = sub("bp Ins repeat$", " Ins", g)
  g = sub("bp Del MH$", " MH", g)
  g
}

# Compact x label: SBS the trinucleotide (ref base in the middle), DBS the alt dinucleotide,
# ID the repeat / microhomology length (5 -> "5+")
.sig_x_label = function(channel, context) {
  if (context == "SBS") return(paste0(substr(channel, 1, 1), substr(channel, 3, 3), substr(channel, 7, 7)))
  if (context == "DBS") return(sub("^..>", "", channel))
  last = sub("^.*:", "", channel)
  ifelse(last == "5", "5+", last)
}

plot_sig_profile = function(matrix, context) {
  if (is.null(matrix) || nrow(matrix) == 0) return(NULL)
  k = sig_context_key(context)
  m = copy(matrix)
  lv = sig_group_levels(k)
  m[, group := factor(group, levels = lv)]
  m = m[!is.na(group)]
  m[, channel := factor(channel, levels = unique(channel))]
  strip_lv = .sig_strip_label(lv, k)
  m[, strip := factor(.sig_strip_label(as.character(group), k), levels = strip_lv)]
  ymax = max(c(m$count, 1))
  band = unique(m[, .(group, strip)])
  band[, `:=`(ymin = ymax * 1.03, ymax = ymax * 1.10)]
  labs_x = setNames(.sig_x_label(as.character(levels(m$channel)), k), levels(m$channel))

  ggplot2::ggplot(m, ggplot2::aes(x = channel, y = count, fill = group)) +
    ggplot2::geom_rect(data = band, inherit.aes = FALSE,
                       ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = group)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_grid(. ~ strip, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = sig_group_colours(k), guide = "none") +
    ggplot2::scale_x_discrete(labels = labs_x) +
    # Integer breaks only: a sparse DBS profile (max 1-2) would otherwise label 0.25 as "0"
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02)),
                                breaks = function(l) { b = pretty(l); b[b == round(b)] },
                                labels = function(x) formatC(x, format = "d", big.mark = ",")) +
    ggplot2::labs(x = NULL, y = "Mutations") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "#ece8df", linewidth = 0.3),
      panel.spacing.x    = grid::unit(if (k == "ID") 1.5 else 3, "pt"),
      strip.text         = ggplot2::element_text(size = if (k == "ID") 6.5 else 8,
                                                 face = "bold", colour = "#1b1e22",
                                                 angle = if (k == "ID") 90 else 0,
                                                 hjust = if (k == "ID") 0 else 0.5),
      strip.clip         = "off",
      axis.text.x        = ggplot2::element_text(size = if (k == "SBS") 4.2 else 5.5, angle = 90,
                                                 vjust = 0.5, hjust = 1, family = "mono",
                                                 colour = "#6a6f76"),
      axis.ticks.x       = ggplot2::element_blank(),
      plot.background    = ggplot2::element_rect(fill = "#ffffff", colour = NA)
    )
}

# ── Discovery ────────────────────────────────────────────────────────────────────────────────
# Assignment results live under .../COSMIC_v<ver>/<SBS96|DBS78|ID83>/Assignment_Solution/; the
# highest version present wins. Matrices are globbed by suffix: the stem is a sanitised sample
# name, not necessarily the report's id, so a stem match is only a tie-breaker.
locate_signatures = function(sample_dir, sample_id) {
  root = file.path(sample_dir, "signatures")
  if (!dir.exists(root)) root = sample_dir
  all = list.files(root, recursive = TRUE, full.names = TRUE,
                   pattern = "(_Activities\\.txt|_Samples_Stats\\.txt|\\.(SBS96|DBS78|ID83)\\.all)$")
  act = all[grepl("Assignment_Solution_Activities\\.txt$", all)]
  sts = all[grepl("Assignment_Solution_Samples_Stats\\.txt$", all)]

  ver_of = function(p) {
    v = regmatches(p, regexpr("COSMIC_v[0-9]+(\\.[0-9]+)*", p))
    if (length(v) == 0) NA_character_ else sub("^COSMIC_v", "", v)
  }
  act_ver = vapply(act, ver_of, character(1), USE.NAMES = FALSE)
  version = NA_character_
  if (any(!is.na(act_ver))) {
    vv = unique(act_ver[!is.na(act_ver)])
    # which.max() cannot coerce a numeric_version; order() can, and sorts 3.10 above 3.6
    version = vv[order(numeric_version(vv), decreasing = TRUE)[1]]
    act = act[!is.na(act_ver) & act_ver == version]
  }

  stem = gsub("[^A-Za-z0-9_-]", "_", sample_id)
  pick_matrix = function(ctx) {
    hits = all[grepl(paste0("\\.", ctx, "\\.all$"), all)]
    if (length(hits) == 0) return(NULL)
    pref = hits[basename(hits) == paste0(stem, ".", ctx, ".all")]
    if (length(pref) > 0) pref[1] else hits[1]
  }

  contexts = list()
  for (k in names(SIG_CONTEXTS)) {
    ctx = SIG_CONTEXTS[[k]]
    seg = paste0("/", ctx, "/Assignment_Solution/")
    a = act[grepl(seg, act, fixed = TRUE)]
    if (length(a) == 0) next
    s = sts[grepl(seg, sts, fixed = TRUE) & startsWith(sts, dirname(dirname(a[1])))]
    contexts[[k]] = list(activities = a[1], stats = if (length(s) > 0) s[1] else NULL,
                         matrix = pick_matrix(ctx))
  }
  list(version = version, contexts = contexts, sample_id = sample_id)
}

parse_signatures = function(inputs, lookup_dir) {
  if (is.null(inputs) || length(inputs$contexts) == 0) return(NULL)
  lookup = if (!is.null(lookup_dir)) load_sig_aetiology(lookup_dir, inputs$version) else NULL
  sid = inputs$sample_id

  contexts = list()
  for (k in names(inputs$contexts)) {
    f   = inputs$contexts[[k]]
    act = parse_spa_activities(f$activities, sid)
    if (is.null(act)) next
    mat = parse_sig_matrix(f$matrix, k, sid)
    st  = parse_spa_stats(f$stats, sid)
    if (is.null(st))
      st = list(total = NA_real_, cosine = NA_real_, l1_pct = NA_real_, l2_pct = NA_real_,
                kl = NA_real_, correlation = NA_real_)
    if (is.na(st$total))
      st$total = if (!is.null(mat)) sum(mat$count) else sum(act$mutations)
    contexts[[k]] = list(
      activities = annotate_sig_activities(act, lookup),
      stats      = st,
      matrix     = mat,
      flags      = sig_context_flags(k, st$total, st$cosine),
      files      = list(activities = f$activities, stats = f$stats, matrix = f$matrix)
    )
  }
  if (length(contexts) == 0) return(NULL)

  list(
    cosmic_version = inputs$version,
    lookup_version = if (is.null(lookup)) NA_character_ else attr(lookup, "version"),
    contexts       = contexts,
    dominant       = sig_dominant(contexts)
  )
}

register_section(list(
  id     = "signatures",
  title  = "Mutational signatures",
  locate = function(sample_dir, sample_id) locate_signatures(sample_dir, sample_id),
  parse  = function(inputs, section_data, ctx) {
    repo = if (!is.null(ctx$repo_dir)) ctx$repo_dir else "."
    parse_signatures(inputs, file.path(repo, "assets", "signatures"))
  }
))
