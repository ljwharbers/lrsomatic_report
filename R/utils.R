suppressPackageStartupMessages({
  library(data.table)
})

ensure_chr_prefix = function(x) {
  ifelse(startsWith(x, "chr"), x, paste0("chr", x))
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
  vapply(extra_vec, function(x) {
    kv = parse_extra_kv(x)
    if (key %in% names(kv)) kv[[key]] else NA_character_
  }, character(1), USE.NAMES = FALSE)
}

# Load a gene panel TSV or plain text file; returns character vector of gene symbols
load_gene_panel = function(path) {
  if (!file.exists(path)) stop("Gene panel file not found: ", path)
  dt = tryCatch(
    fread(path, header = TRUE, sep = "\t", fill = TRUE),
    error = function(e) fread(path, header = FALSE, sep = "\t", fill = TRUE)
  )
  gene_col = if ("gene" %in% tolower(names(dt))) names(dt)[tolower(names(dt)) == "gene"][1] else names(dt)[1]
  unique(dt[[gene_col]])
}

# Filter a data frame to rows where the gene column matches the panel
filter_by_gene_panel = function(dt, panel_genes, gene_col = "gene") {
  dt[dt[[gene_col]] %in% panel_genes, ]
}

# Resolve a --gene-panel arg: either a builtin name ("lymphoid") or a file path
resolve_gene_panel = function(panel_arg, assets_dir) {
  builtin_path = file.path(assets_dir, "gene_lists", paste0(panel_arg, ".tsv"))
  if (file.exists(builtin_path)) return(load_gene_panel(builtin_path))
  if (file.exists(panel_arg)) return(load_gene_panel(panel_arg))
  stop("Gene panel not found (tried builtin '", panel_arg, "' and as file path)")
}

# Load all gene panels from assets/gene_lists/*.tsv
# Returns a named list (name = panel name, value = character vector of gene symbols)
load_all_gene_panels = function(assets_dir) {
  tsv_files = Sys.glob(file.path(assets_dir, "gene_lists", "*.tsv"))
  if (length(tsv_files) == 0) return(list())
  panels = lapply(tsv_files, load_gene_panel)
  names(panels) = tools::file_path_sans_ext(basename(tsv_files))
  panels
}

# Format a number for human-readable display
fmt_bp = function(x) {
  x = as.numeric(x)
  ifelse(abs(x) >= 1e6, paste0(round(x / 1e6, 1), " Mb"),
    ifelse(abs(x) >= 1e3, paste0(round(x / 1e3, 1), " kb"),
      paste0(x, " bp")))
}
