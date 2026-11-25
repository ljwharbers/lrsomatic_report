#!/usr/bin/env Rscript

# Repo root relative to this script: normalizePath() the script file (resolving a bin/ symlink) before dirname()
script_file = normalizePath(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
repo_dir    = normalizePath(file.path(dirname(script_file), ".."))

# --version answers before the library() block: base R only, so it still works when a
# dependency is broken, and stays cheap enough for a per-task Nextflow `eval` call.
if (any(commandArgs(trailingOnly = TRUE) %in% c("--version", "-V"))) {
  cat(readLines(file.path(repo_dir, "VERSION"), warn = FALSE)[1], "\n", sep = "")
  quit(status = 0)
}

suppressPackageStartupMessages({
  library(optparse)
  library(quarto)
  library(yaml)
})

# Source helpers (needed for detect_reference and locate_outputs below)
source(file.path(repo_dir, "R/utils.R"))
source(file.path(repo_dir, "R/references.R"))
source(file.path(repo_dir, "R/locate_outputs.R"))
source(file.path(repo_dir, "R/hmf.R"))

# ---- CLI argument parsing -----------------------------------------------
option_list = list(
  make_option("--sample-dir",  type = "character", default = NULL,
              help = "Path to the sample output directory (required)"),
  make_option("--sample-id",   type = "character", default = NULL,
              help = "Sample identifier, e.g. SAMPLE_ID (required)"),
  make_option("--reference",   type = "character", default = "auto",
              help = "Reference genome: t2t | hg38 | auto (default: auto)"),
  make_option("--sex",         type = "character", default = NULL,
              help = "Biological sex: male | female | XY | XX (required)"),
  make_option("--gene-panel",  type = "character", default = "none",
              help = paste("Gene panel applied on load: none | builtin name (lymphoid) | path to TSV",
                           "(default: none, i.e. unfiltered). Repeatable — pass it several times to",
                           "open with several panels applied at once; a variant or SV is kept if it",
                           "hits any of them.")),
  make_option("--gene-lists-dir", type = "character", default = NULL,
              help = paste("Directory of builtin gene panel TSVs, named <panel>.tsv or",
                           "<panel>.<hg38|t2t>.tsv (default: the bundled assets/gene_lists).",
                           "Replaces the bundled set rather than adding to it.")),
  make_option("--embed",       type = "character", default = "slim",
              help = paste("Small-variant rows embedded in the HTML table: slim | full",
                           "(default: slim — HIGH/MODERATE impact, ClinVar pathogenic and",
                           "panel-gene rows). The full table is always written next to the",
                           "report as <sample-id>_small_variants.tsv.gz.")),
  make_option("--embed-full-table", type = "character", default = "yes",
              help = paste("Also embed that full table in the HTML (gzipped, about the file's",
                           "size), behind a download button and a 'show all' action: yes | no",
                           "(default: yes; skipped above --embed-full-table-max-mb).")),
  make_option("--embed-full-table-max-mb", type = "character",
              default = as.character(FULL_TABLE_EMBED_MAX_MB),
              help = paste("Largest full table (MB, gzipped) --embed-full-table embeds; a larger",
                           "one stays only next to the report. The HTML grows by ~4/3 of it",
                           "(default: %default).")),
  make_option("--run-mode",    type = "character", default = NULL,
              help = paste("Run mode: paired | tumour_only. Overrides the mode inferred from",
                           "normal-side QC files (default: inferred).")),
  make_option("--purity-spread-warn", type = "character", default = NULL,
              help = paste("Warn prominently when the purity estimates of ASCAT, Wakhan and SAVANA",
                           "(top solution of each) differ by at least this much (default: 0.2).",
                           "A value above 1 never warns.")),
  make_option("--hmf-driver-panel", type = "character", default = NULL,
              help = paste("HMF DriverGenePanel.38.tsv (oncoanalyser resource bundle: GRCh38, or a",
                           "CHM13 copy lifted by the pipeline, whose Ensembl cache declares REFERENCE t2t).",
                           "Turns on the HMF driver annotation section and adds the panel as the",
                           "selectable gene list 'hmf_driver'. Needs --hmf-ensembl-data-dir.")),
  make_option("--hmf-ensembl-data-dir", type = "character", default = NULL,
              help = paste("HMF Ensembl data cache directory (common/ensembl_data, with",
                           "ensembl_gene_data.csv and ensembl_trans_exon_data.csv): canonical",
                           "transcripts and gene spans for the HMF section. A REFERENCE file in it",
                           "declares the build (absent: GRCh38).")),
  make_option("--hmf-known-fusions", type = "character", default = NULL,
              help = "HMF known_fusion_data.38.csv: known-fusion matching of the SV calls (optional)."),
  make_option("--hmf-hotspots", type = "character", default = NULL,
              help = "HMF KnownHotspots.somatic.38.vcf.gz: hotspot matching of small variants (optional)."),
  make_option("--hmf-actionable-bed", type = "character", default = NULL,
              help = paste("HMF ActionableCodingPanel.38.bed(.gz): flags small variants inside the",
                           "panel's coding regions (optional; absent from bundle v3).")),
  make_option("--version",     action = "store_true", default = FALSE,
              help = "Print the tool version and exit"),
  make_option("--output",      type = "character", default = NULL,
              help = "Output HTML path (default: <sample-id>_report.html in current dir)"),
  make_option("--title",       type = "character", default = NULL,
              help = "Report title (default: 'LRSomatic Report – <sample-id>')")
)

# --gene-panel is repeatable (optparse has no action="append"): strip it from argv before parse_args()
argv        = commandArgs(trailingOnly = TRUE)
gene_panel_args = extract_repeated_option(argv, "--gene-panel")
opt         = parse_args(OptionParser(option_list = option_list), args = gene_panel_args$rest)
gene_panels = if (length(gene_panel_args$values) == 0) "none" else gene_panel_args$values

# ---- Validate required arguments ----------------------------------------
abort = function(...) { cat("ERROR:", ..., "\n"); quit(status = 1) }

if (is.null(opt[["sample-dir"]])) abort("--sample-dir is required")
if (is.null(opt[["sex"]]))        abort("--sex is required")

sample_dir  = normalizePath(opt[["sample-dir"]], mustWork = TRUE)
sample_id   = if (!is.null(opt[["sample-id"]])) opt[["sample-id"]] else basename(sample_dir)
sex         = tolower(trimws(opt[["sex"]]))
sex         = switch(sex, xy = "male", xx = "female", sex)  # normalise XY/XX

embed       = tolower(trimws(opt[["embed"]]))
if (!embed %in% c("slim", "full"))
  abort(paste0("--embed must be 'slim' or 'full' (got '", opt[["embed"]], "')"))

embed_full_table = parse_yes_no(opt[["embed-full-table"]])
if (is.na(embed_full_table))
  abort(paste0("--embed-full-table must be 'yes' or 'no' (got '", opt[["embed-full-table"]], "')"))
embed_full_table_max_mb = parse_positive_mb(opt[["embed-full-table-max-mb"]])
if (is.na(embed_full_table_max_mb))
  abort(paste0("--embed-full-table-max-mb must be a number above 0 (got '",
               opt[["embed-full-table-max-mb"]], "')"))

purity_spread_warn = opt[["purity-spread-warn"]]
if (!is.null(purity_spread_warn)) {
  v = suppressWarnings(as.numeric(purity_spread_warn))
  if (length(v) != 1 || is.na(v) || v <= 0)
    abort(paste0("--purity-spread-warn must be a positive number (got '", purity_spread_warn, "')"))
  purity_spread_warn = v
}

run_mode   = opt[["run-mode"]]
if (!is.null(run_mode) && !tolower(trimws(run_mode)) %in% names(RUN_MODES))
  abort(paste0("--run-mode must be 'paired' or 'tumour_only' (got '", run_mode, "')"))

output     = if (!is.null(opt[["output"]])) opt[["output"]] else
             file.path(getwd(), paste0(sample_id, "_report.html"))
title      = if (!is.null(opt[["title"]])) opt[["title"]] else
             paste0("LRSomatic Report – ", sample_id)

# ---- HMF resources: checked before anything renders -----------------------
hmf_opts = c(driver_panel = "hmf-driver-panel", ensembl_data_dir = "hmf-ensembl-data-dir",
             known_fusions = "hmf-known-fusions", hotspots = "hmf-hotspots",
             actionable_bed = "hmf-actionable-bed")
hmf_files = lapply(hmf_opts, function(o) opt[[o]])
hmf_files = Filter(Negate(is.null), hmf_files)
if (length(hmf_files) > 0) {
  if (is.null(hmf_files$driver_panel))
    abort(paste0("--", hmf_opts[names(hmf_files)][1], " needs --hmf-driver-panel"))
  if (is.null(hmf_files$ensembl_data_dir))
    abort("--hmf-driver-panel needs --hmf-ensembl-data-dir (the HMF Ensembl data cache)")
  for (nm in names(hmf_files)) {
    p = hmf_files[[nm]]
    ok = if (nm == "ensembl_data_dir") dir.exists(p) else file.exists(p) && !dir.exists(p)
    if (!ok) abort(paste0("--", hmf_opts[[nm]], " not found: ", p))
    hmf_files[[nm]] = normalizePath(p)
  }
  sniff = function(path, n = 1) {
    con = if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con)); x = readLines(con, n = n, warn = FALSE)
    if (length(x) == 0) "" else x   # an empty file fails the checks below, not strsplit()
  }
  hdr = sniff(hmf_files$driver_panel)
  if (!all(c("gene", "likelihoodType") %in% strsplit(hdr, "\t")[[1]]))
    abort("--hmf-driver-panel does not look like an HMF DriverGenePanel (no gene / likelihoodType columns): ",
          hmf_files$driver_panel)
  for (f in c("ensembl_gene_data.csv", "ensembl_trans_exon_data.csv"))
    if (!file.exists(file.path(hmf_files$ensembl_data_dir, f)))
      abort("--hmf-ensembl-data-dir has no ", f, ": ", hmf_files$ensembl_data_dir)
  if (!is.null(hmf_files$known_fusions) &&
      !all(c("Type", "FiveGene", "ThreeGene") %in% strsplit(sniff(hmf_files$known_fusions), ",")[[1]]))
    abort("--hmf-known-fusions does not look like an HMF known_fusion_data.csv: ", hmf_files$known_fusions)
  if (!is.null(hmf_files$hotspots) && !startsWith(sniff(hmf_files$hotspots), "##fileformat=VCF"))
    abort("--hmf-hotspots is not a VCF: ", hmf_files$hotspots)
}

# ---- Locate per-tool outputs ---------------------------------------------
message("Locating outputs in: ", sample_dir)
outputs = apply_run_mode(locate_outputs(sample_dir, sample_id), run_mode)
message("Run mode: ", outputs$mode, if (identical(outputs$mode_source, "cli")) " (--run-mode)" else " (inferred from QC files)")
message("VEP somatic: ", ifelse(is.null(outputs$vep_somatic), "NOT FOUND", outputs$vep_somatic))
message("Somatic VAF VCF: ", ifelse(is.null(outputs$somatic_vaf_vcf), "NOT FOUND",
                                    paste(outputs$somatic_vaf_vcf, collapse = ", ")))
message("ASCAT segments: ", ifelse(is.null(outputs$ascat_segments), "NOT FOUND", outputs$ascat_segments))

# ---- Auto-detect reference (before the reference-specific gene panels) ----
reference = opt[["reference"]]
if (reference == "auto") {
  # Reuse already-resolved paths rather than a fixed vep/somatic/* glob
  vep_file = outputs$vep_somatic
  sv_file  = if (is.null(vep_file)) {
    hits = list.files(sample_dir, pattern = "severus_somatic\\.vcf\\.gz$", recursive = TRUE, full.names = TRUE)
    if (length(hits) > 0) hits[1] else NA_character_
  } else NA_character_
  probe    = if (!is.null(vep_file)) vep_file else if (!is.na(sv_file)) sv_file else NA_character_
  reference = if (!is.na(probe)) detect_reference(probe) else "t2t"
  message("Auto-detected reference: ", reference)
}
reference = tolower(reference)

# HMF resources must be on the report's reference: the published bundle is GRCh38 (and declares
# nothing); a CHM13 copy lifted by the pipeline declares t2t (hmf_declared_reference()). Plain
# GRCh38 files on another reference are left out with a warning; a declared mismatch is an error
if (length(hmf_files) > 0) {
  bad = hmf_check_reference(hmf_files, reference)
  if (!is.null(bad) && bad$kind == "mismatch") abort(paste0("HMF resources: ", bad$msg, "."))
  if (!is.null(bad)) {
    message("WARNING: ", bad$msg, ", so the HMF section is skipped. On CHM13, use the copy the ",
            "pipeline lifts over (its Ensembl data cache declares REFERENCE t2t).")
    hmf_files = list()
  }
}

# ---- Load all gene panels: every builtin ships in the report, --gene-panel only sets which are checked on load; "__all__" is the internal sentinel for none ----
gene_lists_dir = if (!is.null(opt[["gene-lists-dir"]])) {
  normalizePath(opt[["gene-lists-dir"]], mustWork = TRUE)
} else {
  file.path(repo_dir, "assets", "gene_lists")
}
all_panels = load_all_gene_panels(gene_lists_dir, reference)

# The HMF driver panel as a selectable gene list (unticked on load): spans from the HMF
# Ensembl cache, role from likelihoodType, read through the same loader as every list
if (!is.null(hmf_files$driver_panel)) {
  hmf_tsv = hmf_gene_list_tsv(hmf_files$driver_panel, hmf_files$ensembl_data_dir,
                              tempfile("hmf_driver", fileext = ".tsv"), reference)
  nm = unique_panel_name("hmf_driver", names(all_panels))
  all_panels[[nm]] = tryCatch(load_gene_panel(hmf_tsv, reference),
                              error = function(e) abort(conditionMessage(e)))
  all_panels[[nm]]$name = nm
  message("HMF driver annotation on: ", paste(names(hmf_files), collapse = ", "))
}

# "none" combined with a real panel is contradictory
if (any(vapply(gene_panels, is_no_gene_panel, logical(1))) && length(gene_panels) > 1) {
  abort("--gene-panel none cannot be combined with other panels; drop the 'none'.")
}
# Deduplicate by resolved path
gene_panels = unique(vapply(gene_panels, function(g)
  if (file.exists(g)) normalizePath(g) else g, character(1)))

default_panels = character(0)
for (gp in gene_panels) {
  if (is_no_gene_panel(gp)) next
  if (!is.null(builtin_panel_path(gene_lists_dir, gp, reference))) {
    # Load here so a builtin that fails against this reference aborts before the render
    invisible(tryCatch(resolve_gene_panel(gp, gene_lists_dir, reference),
                       error = function(e) abort(conditionMessage(e))))
    default_panels = c(default_panels, gp)
  } else if (file.exists(gp)) {
    # User-supplied TSV: register alongside the builtins so it is selectable in the report
    nm = unique_panel_name(tools::file_path_sans_ext(basename(gp)), names(all_panels))
    all_panels[[nm]] = tryCatch(load_gene_panel(gp, reference),
                                error = function(e) abort(conditionMessage(e)))
    default_panels = c(default_panels, nm)
  } else {
    abort(paste0("--gene-panel not found: tried builtin '", gp,
                 "' and as a file path. Use 'none' for no filtering."))
  }
}
default_panels = unique(default_panels)
# Sentinel rather than character(0): an empty vector round-trips through execute_params as NULL
if (length(default_panels) == 0) default_panels = "__all__"
message("Gene panels selected on load: ",
        if (identical(default_panels, "__all__")) "none" else paste(default_panels, collapse = ", "))

# ---- Render: copy templates/ and assets/ to a writable dir (the repo may be read-only and Quarto writes next to the .qmd) ----
work = file.path(getwd(), "._render")
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)
invisible(file.copy(file.path(repo_dir, "templates"), work, recursive = TRUE))
invisible(file.copy(file.path(repo_dir, "assets"), work, recursive = TRUE))
template = file.path(work, "templates", "per_sample.qmd")
if (!file.exists(template)) abort("Quarto template not found: ", template)

message("Rendering report to: ", output)
quarto::quarto_render(
  input          = template,
  output_file    = basename(output),
  output_format  = "html",
  execute_params = list(
    sample_id     = sample_id,
    sample_dir    = sample_dir,
    reference     = reference,
    sex           = sex,
    default_panels = default_panels,
    all_panels    = all_panels,
    title         = title,
    repo_dir      = repo_dir,
    outputs       = outputs,
    embed         = embed,
    embed_full_table = embed_full_table,
    embed_full_table_max_mb = embed_full_table_max_mb,
    purity_spread_warn = purity_spread_warn,
    # NULL (section off) or a named list of paths; plain strings survive the YAML round trip
    hmf           = if (length(hmf_files) > 0) hmf_files else NULL
  ),
  quiet = FALSE
)

# Move a file Quarto/knitr wrote next to the template to `dest`
move_out = function(rendered, dest) {
  if (!file.exists(rendered)) return(invisible(FALSE))
  if (normalizePath(rendered, mustWork = FALSE) == normalizePath(dest, mustWork = FALSE))
    return(invisible(TRUE))
  ok = file.copy(rendered, dest, overwrite = TRUE)
  if (ok) file.remove(rendered)
  invisible(ok)
}

# Move output if Quarto wrote it next to the template
move_out(file.path(dirname(template), basename(output)), output)

# The full small-variant table sidecar, written by the setup chunk next to the template
sidecar = paste0(sample_id, "_small_variants.tsv.gz")
if (move_out(file.path(dirname(template), sidecar), file.path(dirname(output), sidecar)))
  message("Full small-variant table written to: ", file.path(dirname(output), sidecar))
unlink(work, recursive = TRUE)

if (file.exists(output)) {
  message("Report written to: ", output)
} else {
  abort("Rendering completed but output file not found at: ", output)
}
