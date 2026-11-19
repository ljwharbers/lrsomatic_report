suppressPackageStartupMessages({
  library(data.table)
})

# Type-stable: fread() may read CHROM as integer and ifelse() on all-NA returns logical
ensure_chr_prefix = function(x) {
  x   = as.character(x)
  out = as.character(ifelse(startsWith(x, "chr"), x, paste0("chr", x)))
  out[is.na(x)] = NA_character_
  out
}

strip_chr_prefix = function(x) {
  sub("^chr", "", x)
}

# Parse VEP "Extra" key=value semicolon-delimited field into a named character vector
parse_extra_kv = function(extra_string) {
  if (is.na(extra_string) || extra_string == "" || extra_string == "-") return(character(0))
  pairs = strsplit(extra_string, ";", fixed = TRUE)[[1]]
  kv = strsplit(pairs, "=", fixed = TRUE)
  keys = vapply(kv, `[`, character(1), 1)
  vals = vapply(kv, function(x) if (length(x) >= 2) paste(x[-1], collapse = "=") else "", character(1))
  setNames(vals, keys)
}

# Vectorised: extract one key from VEP Extra column for each row
extract_extra_key = function(extra_vec, key) {
  extract_extra_keys(extra_vec, key)[[key]]
}

# Several keys in one pass over the Extra column: returns a named list of character vectors (NA where the key is absent). One split per row rather than one per row per key, which is what keeps the text path tractable with a dozen keys
extract_extra_keys = function(extra_vec, keys) {
  extra_vec = as.character(extra_vec)
  n = length(extra_vec)
  out = setNames(rep(list(rep(NA_character_, n)), length(keys)), keys)
  if (n == 0 || length(keys) == 0) return(out)

  pairs = strsplit(extra_vec, ";", fixed = TRUE)
  lens  = lengths(pairs)
  flat  = unlist(pairs, use.names = FALSE)
  row   = rep.int(seq_len(n), lens)
  eq    = regexpr("=", flat, fixed = TRUE)
  has   = eq > 0
  key   = substr(flat[has], 1L, eq[has] - 1L)
  val   = substr(flat[has], eq[has] + 1L, nchar(flat[has]))
  row   = row[has]
  for (k in keys) {
    hit = which(key == k)
    if (length(hit)) out[[k]][row[hit]] = val[hit]
  }
  out
}

# ---- Gene panels: plain lists (they round-trip through Quarto execute_params) of name, path, reference, has_coords, has_scopes, genes/genes_snv/genes_sv and, when has_coords, parallel chrom/start/end/interval_gene ----
#
# An optional `applies_to` column scopes each gene to one table ("snv" | "sv" | blank for
# both). `genes` stays the union; `genes_snv` and `genes_sv` are what the two tables match
# on, and the interval vectors hold the SV-scoped rows only — which is what scopes the
# whole coordinate path (panel_intervals, sv_panel_hits, bnd_panel_genes, the JS payload)
# without any of them knowing the column exists.

# Canonical reference names; a coordinate panel declares one and it is checked against the render
normalise_reference_name = function(x) {
  if (is.null(x) || length(x) != 1 || is.na(x) || !nzchar(trimws(x))) return(NA_character_)
  x = tolower(trimws(x))
  if (x %in% c("t2t", "chm13", "chm13v2", "chm13v2.0", "t2t-chm13")) return("t2t")
  if (x %in% c("hg38", "grch38", "hg38-noalt")) return("hg38")
  x
}

# Reference suffix of a builtin panel filename; list(name, reference), reference NA if unsuffixed
.split_panel_filename = function(path) {
  stem  = tools::file_path_sans_ext(basename(path))
  parts = strsplit(stem, ".", fixed = TRUE)[[1]]
  if (length(parts) >= 2) {
    ref = normalise_reference_name(parts[length(parts)])
    if (!is.na(ref) && ref %in% c("t2t", "hg38"))
      return(list(name = paste(parts[-length(parts)], collapse = "."), reference = ref))
  }
  list(name = stem, reference = NA_character_)
}

# Read the leading "#" comment block of a panel TSV (may declare "# reference: hg38")
.panel_header = function(path) {
  lines = readLines(path, warn = FALSE)
  n_comment = 0L
  declared  = NA_character_
  for (ln in lines) {
    if (!startsWith(ln, "#")) break
    n_comment = n_comment + 1L
    m = regmatches(ln, regexpr("^#\\s*reference\\s*:\\s*\\S+", ln))
    if (length(m) > 0) declared = sub("^#\\s*reference\\s*:\\s*", "", m)
  }
  list(n_comment = n_comment, reference = declared)
}

# Normalise the optional `applies_to` column to "both" | "snv" | "sv", one per row of the
# panel. Blank (or NA) is "both", so a column of mostly-empty cells reads the way it looks.
# An unrecognised value is an error rather than a silent default: the whole point of
# resolve_gene_panel()'s strictness is that a typo must not quietly change what is filtered.
.panel_scope = function(col, genes, path) {
  v = tolower(trimws(as.character(col)))
  v[is.na(v)] = ""
  out = rep(NA_character_, length(v))
  out[v %in% c("", "both", "all")]       = "both"
  out[v %in% c("snv", "small", "snv/indel", "smallvariant", "smallvariants")] = "snv"
  out[v %in% c("sv", "structural")]      = "sv"
  if (anyNA(out)) {
    bad = which(is.na(out))
    stop("Gene panel ", path, " has unrecognised applies_to value(s): ",
         paste(sprintf("'%s' (%s)", v[bad], genes[bad])[seq_len(min(5, length(bad)))],
               collapse = ", "),
         if (length(bad) > 5) paste0(" (and ", length(bad) - 5L, " more)") else "",
         ". Use 'snv', 'sv', or leave blank for both.")
  }
  out
}

# Normalise the optional `role` column to "oncogene" | "tsg" | "both" | "", one per row.
# Blank/NA is "" (no role declared). Like applies_to, an unrecognised value is an error.
.panel_role = function(col, genes, path) {
  v = tolower(trimws(as.character(col)))
  v[is.na(v)] = ""
  v = gsub("[[:space:]_-]+", " ", v)
  out = rep(NA_character_, length(v))
  out[v %in% c("", "na", "none")]                                        = ""
  out[v %in% c("og", "oncogene", "onc")]                                 = "oncogene"
  out[v %in% c("tsg", "tumour suppressor", "tumor suppressor",
               "tumour suppressor gene", "tumor suppressor gene")]       = "tsg"
  out[v %in% c("both", "og/tsg", "tsg/og", "og / tsg", "tsg / og",
               "oncogene/tsg", "tsg/oncogene")]                          = "both"
  if (anyNA(out)) {
    bad = which(is.na(out))
    stop("Gene panel ", path, " has unrecognised role value(s): ",
         paste(sprintf("'%s' (%s)", v[bad], genes[bad])[seq_len(min(5, length(bad)))],
               collapse = ", "),
         if (length(bad) > 5) paste0(" (and ", length(bad) - 5L, " more)") else "",
         ". Use 'OG'/'oncogene', 'TSG'/'tumour suppressor', 'both' ('OG/TSG'), or leave blank.")
  }
  out
}

# Load a gene panel TSV: `gene` column required; chrom/start/end all-or-nothing; a coordinate panel declaring another reference errors, one declaring none loads as ""
#
# Fields added in v2 (all additive — a panel without these columns loads as before):
#   has_roles, role_gene / role   parallel vectors over genes with a NON-EMPTY role, one entry
#                                 per gene (first non-empty role wins); empty when no `role`
#                                 column. role ∈ "oncogene" | "tsg" | "both". See panel_roles().
#   all_chrom / all_start / all_end / all_gene
#                                 every row with valid coordinates REGARDLESS of applies_to
#                                 (the interval vectors above stay SV-scoped); empty for a
#                                 symbol-only panel. See panel_gene_intervals(panel, "all").
# Empty vectors come back as NULL / list() after the Quarto execute_params YAML round trip,
# so consumers must unlist() every vector field.
load_gene_panel = function(path, reference = NULL) {
  if (!file.exists(path)) stop("Gene panel file not found: ", path)

  hdr = .panel_header(path)
  dt = tryCatch(
    fread(path, header = TRUE, sep = "\t", fill = TRUE, skip = hdr$n_comment),
    error = function(e) fread(path, header = FALSE, sep = "\t", fill = TRUE,
                              skip = hdr$n_comment)
  )
  if (nrow(dt) == 0 && ncol(dt) == 0) stop("Gene panel file is empty: ", path)

  # copy(): setnames() rewrites the names vector in place
  orig_names = copy(names(dt))
  setnames(dt, tolower(names(dt)))
  setnames(dt, old = c("chr", "chromosome"), new = c("chrom", "chrom"), skip_absent = TRUE)

  if ("gene" %in% names(dt)) {
    genes = as.character(dt[["gene"]])
  } else if (ncol(dt) == 1) {
    # Headerless one-column list: fread consumed the first symbol as the column name, put it back
    genes = c(orig_names[1], as.character(dt[[1]]))
  } else {
    genes = as.character(dt[[1]])
  }
  genes = trimws(genes)
  keep  = nzchar(genes) & !is.na(genes) & genes != "-"

  # Which table each gene filters. Absent column, or a blank cell, means both — so every
  # panel written before this column existed is unaffected.
  scope = if ("applies_to" %in% names(dt))
            .panel_scope(dt[["applies_to"]], genes, path) else rep("both", length(genes))
  snv_rows = keep & scope %in% c("both", "snv")
  sv_rows  = keep & scope %in% c("both", "sv")

  coord_cols = c("chrom", "start", "end")
  present    = intersect(coord_cols, names(dt))
  if (length(present) > 0 && length(present) < 3) {
    stop("Gene panel ", path, " carries coordinate column(s) ",
         paste(present, collapse = ", "), " but not all of ",
         paste(coord_cols, collapse = ", "),
         ". Supply all three, or none for symbol-only matching.")
  }
  # A panel whose every row is `snv`-scoped has nothing to match SVs on, coordinate
  # columns or not, so it is not in coordinate mode.
  has_coords = length(present) == 3 && any(sv_rows)

  # A `reference` column is an alternative to the "# reference:" comment line.
  declared = hdr$reference
  if (is.na(declared) && "reference" %in% names(dt)) {
    vals = unique(trimws(as.character(dt[["reference"]])))
    vals = vals[nzchar(vals) & !is.na(vals)]
    if (length(vals) > 1)
      stop("Gene panel ", path, " declares more than one reference: ",
           paste(vals, collapse = ", "))
    if (length(vals) == 1) declared = vals
  }
  declared_norm = normalise_reference_name(declared)
  want          = normalise_reference_name(reference)

  # Only coordinate panels are reference-specific; a symbol-only panel is agnostic.
  if (has_coords && !is.na(declared_norm) && !is.na(want) && declared_norm != want) {
    stop("Gene panel ", path, " declares reference '", declared_norm,
         "' but the report is being rendered against '", want,
         "'. Panel coordinates are only valid for the reference they were built on.")
  }

  out = list(
    name       = .split_panel_filename(path)$name,
    path       = path,
    reference  = if (is.na(declared_norm)) "" else declared_norm,
    has_coords = has_coords,
    has_scopes = "applies_to" %in% names(dt),
    # The union, kept so every reader predating the scope column still sees one gene set
    genes      = unique(genes[keep]),
    genes_snv  = unique(genes[snv_rows]),
    genes_sv   = unique(genes[sv_rows])
  )

  if (has_coords) {
    # nzchar() on the *raw* column, before ensure_chr_prefix(): that turns "" into "chr",
    # which is nzchar, so prefixing first would let a blank chromosome through
    chrom_raw = trimws(as.character(dt[["chrom"]]))
    start = suppressWarnings(as.integer(dt[["start"]]))
    end   = suppressWarnings(as.integer(dt[["end"]]))
    # Only rows that can match an SV need coordinates; an `snv` row is matched on its
    # symbol alone and may leave them blank.
    bad   = sv_rows & (is.na(chrom_raw) | !nzchar(chrom_raw) | is.na(start) | is.na(end))
    if (any(bad))
      stop("Gene panel ", path, " has missing or non-numeric coordinates for: ",
           paste(utils::head(genes[bad], 5), collapse = ", "),
           if (sum(bad) > 5) paste0(" (and ", sum(bad) - 5L, " more)") else "",
           ". Only rows with applies_to = 'snv' may leave them blank.")
    # Intervals are the SV-scoped rows only, which is what scopes every coordinate
    # consumer — panel_intervals(), sv_panel_hits(), bnd_panel_genes(), the JS payload
    out$chrom = ensure_chr_prefix(chrom_raw[sv_rows])
    out$start = start[sv_rows]
    out$end   = end[sv_rows]
    # `genes` is deduplicated; the interval vectors are not (one symbol can carry several loci)
    out$interval_gene = genes[sv_rows]
  }

  # Every row with usable coordinates, whatever its applies_to — for consumers that place
  # genes on the genome rather than match SVs (gene-level copy number, overview). Invalid
  # rows are skipped here; only SV-scoped rows are required to carry coordinates (above).
  out$all_chrom = character(0)
  out$all_start = integer(0)
  out$all_end   = integer(0)
  out$all_gene  = character(0)
  if (length(present) == 3) {
    chrom_raw = trimws(as.character(dt[["chrom"]]))
    start = suppressWarnings(as.integer(dt[["start"]]))
    end   = suppressWarnings(as.integer(dt[["end"]]))
    ok = keep & !is.na(chrom_raw) & nzchar(chrom_raw) & !is.na(start) & !is.na(end)
    out$all_chrom = ensure_chr_prefix(chrom_raw[ok])
    out$all_start = start[ok]
    out$all_end   = end[ok]
    out$all_gene  = genes[ok]
  }

  # Optional per-gene role (oncogene / tsg / both); only genes with a declared role are listed
  out$has_roles = "role" %in% names(dt)
  out$role_gene = character(0)
  out$role      = character(0)
  if (out$has_roles) {
    role = .panel_role(dt[["role"]], genes, path)
    has  = keep & nzchar(role)
    rg   = genes[has]
    first = !duplicated(rg)
    out$role_gene = rg[first]
    out$role      = role[has][first]
  }

  out
}

# Panel intervals as a keyed data.table, or NULL for a symbol-only panel.
panel_intervals = function(panel) {
  if (is.null(panel) || !isTRUE(panel$has_coords)) return(NULL)
  dt = data.table(gene  = as.character(panel$interval_gene),
                  chrom = as.character(panel$chrom),
                  start = as.integer(panel$start),
                  end   = as.integer(panel$end))
  setkey(dt, chrom, start, end)
  dt
}

# Panel gene intervals as an (unkeyed) data.table(chrom, start, end, gene); zero rows when the
# panel carries none. scope "sv" = the SV-scoped intervals panel_intervals() uses (only when
# has_coords); "all" = every row with valid coordinates regardless of applies_to. unlist()
# throughout: the panel may have come through the execute_params YAML round trip.
panel_gene_intervals = function(panel, scope = c("sv", "all")) {
  scope = match.arg(scope)
  empty = data.table(chrom = character(), start = integer(), end = integer(),
                     gene = character())
  if (is.null(panel)) return(empty)
  f = if (scope == "sv") {
    if (!isTRUE(as.logical(unlist(panel$has_coords)))) return(empty)
    c("chrom", "start", "end", "interval_gene")
  } else c("all_chrom", "all_start", "all_end", "all_gene")
  v = lapply(f, function(k) unlist(panel[[k]], use.names = FALSE))
  if (length(v[[1]]) == 0) return(empty)
  data.table(chrom = as.character(v[[1]]), start = as.integer(v[[2]]),
             end = as.integer(v[[3]]), gene = as.character(v[[4]]))
}

# Declared roles of a panel as data.table(gene, role), zero rows when it declares none
panel_roles = function(panel) {
  g = as.character(unlist(panel$role_gene, use.names = FALSE))
  r = as.character(unlist(panel$role, use.names = FALSE))
  if (length(g) == 0 || length(g) != length(r))
    return(data.table(gene = character(), role = character()))
  data.table(gene = g, role = r)
}

# Is a --gene-panel argument the "no filtering" sentinel?
is_no_gene_panel = function(panel_arg) {
  is.null(panel_arg) || length(panel_arg) != 1 || is.na(panel_arg) ||
    identical(tolower(trimws(panel_arg)), "none")
}

# Path of a builtin panel, preferring the `reference`-specific variant; NULL if none
builtin_panel_path = function(gene_lists_dir, name, reference = NULL) {
  ref = normalise_reference_name(reference)
  candidates = c(if (!is.na(ref)) file.path(gene_lists_dir, paste0(name, ".", ref, ".tsv")),
                 file.path(gene_lists_dir, paste0(name, ".tsv")))
  hit = candidates[file.exists(candidates)]
  if (length(hit) > 0) hit[1] else NULL
}

# Resolve a --gene-panel arg: "none", a builtin name, or a TSV path; anything else errors
resolve_gene_panel = function(panel_arg, gene_lists_dir, reference = NULL) {
  if (is_no_gene_panel(panel_arg)) return(NULL)
  builtin = builtin_panel_path(gene_lists_dir, panel_arg, reference)
  if (!is.null(builtin)) return(load_gene_panel(builtin, reference))
  if (file.exists(panel_arg)) return(load_gene_panel(panel_arg, reference))
  stop("Gene panel not found (tried builtin '", panel_arg, "' and as file path)")
}

# Load all builtin panels, resolving reference-specific files to one entry and skipping panels not shipped for this reference
load_all_gene_panels = function(gene_lists_dir, reference = NULL) {
  tsv_files = Sys.glob(file.path(gene_lists_dir, "*.tsv"))
  if (length(tsv_files) == 0) return(list())

  meta = lapply(tsv_files, .split_panel_filename)
  ref  = normalise_reference_name(reference)

  panels = list()
  for (nm in unique(vapply(meta, `[[`, character(1), "name"))) {
    idx  = which(vapply(meta, `[[`, character(1), "name") == nm)
    refs = vapply(meta[idx], function(m) m$reference, character(1))
    pick = if (!is.na(ref) && any(refs == ref, na.rm = TRUE)) idx[which(refs == ref)[1]]
           else if (any(is.na(refs))) idx[which(is.na(refs))[1]]
           else NA_integer_
    if (is.na(pick)) {
      message("Gene panel '", nm, "' ships only for reference(s) ",
              paste(unique(refs), collapse = ", "), " — not offered for ",
              if (is.na(ref)) "an unknown reference" else ref)
      next
    }
    p = tryCatch(load_gene_panel(tsv_files[pick], reference),
      error = function(e) {
        message("Skipping gene panel '", nm, "': ", conditionMessage(e)); NULL })
    if (!is.null(p)) { p$builtin = TRUE; panels[[nm]] = p }
  }
  panels
}

# Display label of a panel key: a builtin's name capitalised ("lymphoid" -> "Lymphoid"), a
# user list's file stem as given. Vectorised over `keys`; `panels` is the named panel list
panel_label = function(keys, panels = list()) {
  keys = as.character(keys)
  builtin = vapply(keys, function(k) nzchar(k) && isTRUE(as.logical(unlist(panels[[k]]$builtin))),
                   logical(1), USE.NAMES = FALSE)
  ifelse(builtin & nzchar(keys), paste0(toupper(substr(keys, 1, 1)), substr(keys, 2, nchar(keys))), keys)
}

# Non-colliding key for a user panel: the first collision gets "-custom", further ones are numbered
unique_panel_name = function(nm, taken) {
  if (!(nm %in% taken)) return(nm)
  cand = paste0(nm, "-custom")
  i = 1L
  while (cand %in% taken) {
    i = i + 1L
    cand = paste0(nm, "-custom", i)
  }
  cand
}

# Pull every occurrence of a repeatable flag out of argv (optparse has no action="append"); accepts `--flag value` and `--flag=value`
extract_repeated_option = function(args, flag) {
  args = as.character(args)
  vals = character(0)
  rest = character(0)
  i = 1L
  while (i <= length(args)) {
    a = args[i]
    if (identical(a, flag)) {
      if (i == length(args)) stop(flag, " requires a value")
      vals = c(vals, args[i + 1L])
      i = i + 2L
    } else if (startsWith(a, paste0(flag, "="))) {
      vals = c(vals, substring(a, nchar(flag) + 2L))
      i = i + 1L
    } else {
      rest = c(rest, a)
      i = i + 1L
    }
  }
  list(values = vals, rest = rest)
}

# ---- Small JS serialisation helpers: the R table and the JS indexing it are generated from the same object ----

js_quote = function(x) paste0('"', gsub('"', '\\\\"', as.character(x)), '"')

# {"col":0,"other":1} — column name to zero-based index, for a rownames=FALSE DT.
js_col_index_map = function(nms) {
  if (length(nms) == 0) return("{}")
  paste0("{", paste0(js_quote(nms), ":", seq_along(nms) - 1L, collapse = ","), "}")
}

# One JS number literal per element; element-wise so format() cannot pad or switch to scientific notation
js_num = function(x) {
  vapply(x, function(v) {
    if (is.na(v)) "null" else format(v, scientific = FALSE, trim = TRUE)
  }, character(1), USE.NAMES = FALSE)
}

.js_cell = function(v) if (is.numeric(v)) js_num(v) else ifelse(is.na(v), "null", js_quote(v))

# JS array literal from an atomic vector (numbers bare, NA as null); hand-written to avoid a jsonlite dependency
js_vec = function(x) {
  if (length(x) == 0) return("[]")
  paste0("[", paste(.js_cell(x), collapse = ","), "]")
}

# JS array of arrays, one per row of `dt`, columns in `cols` order (row-major to keep the payload small)
js_rows = function(dt, cols) {
  if (is.null(dt) || nrow(dt) == 0) return("[]")
  cells = lapply(cols, function(cl) .js_cell(dt[[cl]]))
  rows = do.call(paste, c(cells, sep = ","))
  paste0("[[", paste(rows, collapse = "],["), "]]")
}

# VEP impact severity order, as used by the styleEqual() palettes
IMPACT_LEVELS = c("HIGH", "MODERATE", "LOW", "MODIFIER")

# Distinct-value bounds for offering a dropdown; outside them a column keeps its text filter
FACET_MIN_VALUES = 2L
FACET_MAX_VALUES = 200L

# Distinct values with row counts for the tickbox column filters (assets/js/facet_filter.js): `seps` splits a column into tokens, `levels` fixes value order; returns a JS object literal {col: {sep, values: [[value, n], ...]}} with null as the "(none)" bucket. Columns absent or outside the FACET_*_VALUES bounds are omitted with a message
js_facet_defs = function(dt, cols, seps = character(0), levels = list()) {
  if (is.null(dt) || nrow(dt) == 0 || length(cols) == 0) return("{}")
  entries = character(0)

  for (nm in cols) {
    if (!nm %in% names(dt)) {
      message("Facet column '", nm, "' is not in the table - no value filter for it.")
      next
    }
    sp = if (nm %in% names(seps)) seps[[nm]] else NA_character_
    v  = as.character(dt[[nm]])

    if (!is.na(sp) && nzchar(sp)) {
      lst  = strsplit(v, sp, fixed = TRUE)
      lens = lengths(lst)
      # A cell that splits into nothing at all still contributes its row to the NA bucket.
      lst[lens == 0L] = NA_character_
      d = data.table(row = rep.int(seq_along(v), pmax(lens, 1L)),
                     tok = trimws(unlist(lst, use.names = FALSE)))
      d = unique(d, by = c("row", "tok"))   # "a,a" counts its row once
    } else {
      d = data.table(row = seq_along(v), tok = trimws(v))
    }
    d[is.na(tok) | !nzchar(tok), tok := NA_character_]

    cnt   = d[, .N, by = tok]
    n_val = nrow(cnt[!is.na(tok)])
    if (n_val < FACET_MIN_VALUES || n_val > FACET_MAX_VALUES) {
      message("Facet column '", nm, "': ", n_val,
              " distinct value(s) - keeping the plain text filter.")
      next
    }

    # Fixed levels first, then descending count, then alphabetical; the NA bucket sorts last
    lv   = intersect(as.character(levels[[nm]]), cnt$tok[!is.na(cnt$tok)])
    rank = ifelse(is.na(cnt$tok), length(lv) + 2L,
                  ifelse(cnt$tok %in% lv, match(cnt$tok, lv), length(lv) + 1L))
    cnt  = cnt[order(rank, -N, tok)]

    vals = sprintf("[%s,%s]",
                   ifelse(is.na(cnt$tok), "null", js_quote(cnt$tok)),
                   js_num(cnt$N))
    entries = c(entries,
                sprintf('%s:{"sep":%s,"values":[%s]}',
                        js_quote(nm),
                        if (is.na(sp) || !nzchar(sp)) "null" else js_quote(sp),
                        paste(vals, collapse = ",")))
  }

  paste0("{", paste(entries, collapse = ","), "}")
}

# Format a number for human-readable display
fmt_bp = function(x) {
  x = as.numeric(x)
  ifelse(abs(x) >= 1e6, paste0(round(x / 1e6, 1), " Mb"),
    ifelse(abs(x) >= 1e3, paste0(round(x / 1e3, 1), " kb"),
      paste0(x, " bp")))
}

# Embed a PNG as a framed figure (.report-figure in report.scss) with an optional caption
embed_png = function(path, max_width = "900px", caption = NULL) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  b64 = base64enc::base64encode(path)
  htmltools::tags$figure(
    class = "report-figure",
    htmltools::tags$img(
      src   = paste0("data:image/png;base64,", b64),
      style = paste0("max-width:", max_width, "; display:block; margin:auto;")
    ),
    if (!is.null(caption)) htmltools::tags$figcaption(class = "report-figure__caption", caption)
  )
}

# ggplot (or a zero-argument function drawing base graphics) -> svglite SVG -> base64
# <img src="data:image/svg+xml;base64,..."> tag; width/height in inches. Same route as the
# genome circos (svglite to a temp file, base64enc) so the report stays self-contained.
embed_svg_plot = function(p, width, height, alt = "", class = "report-figure") {
  if (is.null(p)) return(NULL)
  tmp = tempfile(fileext = ".svg")
  on.exit(unlink(tmp), add = TRUE)
  svglite::svglite(tmp, width = width, height = height)
  dev = grDevices::dev.cur()
  tryCatch({
    if (is.function(p)) p() else print(p)
  }, finally = grDevices::dev.off(dev))
  htmltools::tags$img(
    class = class,
    src   = paste0("data:image/svg+xml;base64,", base64enc::base64encode(tmp)),
    alt   = alt,
    style = "max-width:100%; height:auto; display:block; margin:auto;"
  )
}

# Map canonical column names to the columns actually present in `dt`.
# `spec` is a named character vector or list: canonical = c("candidate1", "candidate2", ...).
# Pass 1 takes, per canonical, the first candidate present by EXACT name; pass 2 resolves the
# rest by normalised name (lowercase, non-alphanumerics stripped), in candidate order, never
# reusing a column already claimed — so a table carrying both `pos_1` and `Pos_1` maps each
# to its own canonical when both are named exactly. Returns a named character vector
# canonical -> actual column name, NA where nothing matched.
resolve_columns = function(dt, spec) {
  nms   = if (is.character(dt)) dt else names(dt)
  canon = names(spec)
  if (is.null(canon) || any(!nzchar(canon))) stop("resolve_columns(): spec must be fully named")
  out   = setNames(rep(NA_character_, length(canon)), canon)
  taken = character(0)
  # Pass 1: exact names
  for (k in canon) {
    hit = intersect(as.character(unlist(spec[[k]])), setdiff(nms, taken))
    if (length(hit) > 0) { out[[k]] = hit[1]; taken = c(taken, hit[1]) }
  }
  # Pass 2: normalised names
  norm = function(x) gsub("[^a-z0-9]", "", tolower(x))
  for (k in canon[is.na(out)]) {
    free = setdiff(nms, taken)
    for (cand in as.character(unlist(spec[[k]]))) {
      m = free[norm(free) == norm(cand)]
      if (length(m) > 0) { out[[k]] = m[1]; taken = c(taken, m[1]); break }
    }
  }
  out
}

# Small-variant counts over unique chrom/pos/ref/alt: list(total, snv, indel). An SNV has a
# single-base ACGT ref and alt; anything else (incl. VEP's dash-form "-"/"A") is an indel.
# NA counts when vep_data is NULL (no VEP file), zeros for an empty table.
count_small_variants = function(vep_data) {
  if (is.null(vep_data))
    return(list(total = NA_integer_, snv = NA_integer_, indel = NA_integer_))
  if (nrow(vep_data) == 0) return(list(total = 0L, snv = 0L, indel = 0L))
  u = unique(as.data.table(vep_data)[, .(chrom, pos, ref, alt)])
  base = c("A", "C", "G", "T")
  is_snv = toupper(as.character(u$ref)) %in% base & toupper(as.character(u$alt)) %in% base
  list(total = nrow(u), snv = sum(is_snv), indel = nrow(u) - sum(is_snv))
}

# Coding TMB from a variant_table; consequence may be comma-joined; denominator_mb defaults to 30
compute_tmb = function(variant_table, denominator_mb = 30) {
  nonsyn_terms = c(
    "missense_variant", "frameshift_variant", "stop_gained", "stop_lost",
    "start_lost", "inframe_insertion", "inframe_deletion",
    "splice_acceptor_variant", "splice_donor_variant", "protein_altering_variant"
  )
  if (is.null(variant_table) || nrow(variant_table) == 0) {
    return(list(n_nonsyn = NA_integer_, tmb = NA_real_, denominator_mb = denominator_mb))
  }
  is_nonsyn = vapply(variant_table$consequence, function(csq) {
    if (is.na(csq) || csq == "") return(FALSE)
    any(trimws(unlist(strsplit(csq, ","))) %in% nonsyn_terms)
  }, logical(1))
  n_nonsyn = sum(is_nonsyn, na.rm = TRUE)
  list(
    n_nonsyn       = n_nonsyn,
    tmb            = round(n_nonsyn / denominator_mb, 2),
    denominator_mb = denominator_mb
  )
}
