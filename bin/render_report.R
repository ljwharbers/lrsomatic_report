#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(quarto)
  library(yaml)
})

# Locate the repository root relative to this script
script_dir = normalizePath(dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())])))
repo_dir   = normalizePath(file.path(script_dir, ".."))

# Source helpers (needed for detect_reference and locate_outputs below)
source(file.path(repo_dir, "R/utils.R"))
source(file.path(repo_dir, "R/references.R"))
source(file.path(repo_dir, "R/locate_outputs.R"))

# ---- CLI argument parsing -----------------------------------------------
option_list = list(
  make_option("--sample-dir",  type = "character", default = NULL,
              help = "Path to the sample output directory (required)"),
  make_option("--sample-id",   type = "character", default = NULL,
              help = "Sample identifier, e.g. DLBCL3_pooled (required)"),
  make_option("--reference",   type = "character", default = "auto",
              help = "Reference genome: t2t | hg38 | auto (default: auto)"),
  make_option("--sex",         type = "character", default = NULL,
              help = "Biological sex: male | female | XY | XX (required)"),
  make_option("--gene-panel",  type = "character", default = "lymphoid",
              help = "Gene panel: builtin name (lymphoid) or path to TSV (default: lymphoid)"),
  make_option("--output",      type = "character", default = NULL,
              help = "Output HTML path (default: <sample-id>_report.html in current dir)"),
  make_option("--title",       type = "character", default = NULL,
              help = "Report title (default: 'LRSomatic Report – <sample-id>')")
)

opt = parse_args(OptionParser(option_list = option_list))

# ---- Validate required arguments ----------------------------------------
abort = function(...) { cat("ERROR:", ..., "\n"); quit(status = 1) }

if (is.null(opt[["sample-dir"]])) abort("--sample-dir is required")
if (is.null(opt[["sex"]]))        abort("--sex is required")

sample_dir = normalizePath(opt[["sample-dir"]], mustWork = TRUE)
sample_id  = if (!is.null(opt[["sample-id"]])) opt[["sample-id"]] else basename(sample_dir)
sex        = tolower(trimws(opt[["sex"]]))
sex        = switch(sex, xy = "male", xx = "female", sex)  # normalise XY/XX

gene_panel = opt[["gene-panel"]]
output     = if (!is.null(opt[["output"]])) opt[["output"]] else
             file.path(getwd(), paste0(sample_id, "_report.html"))
title      = if (!is.null(opt[["title"]])) opt[["title"]] else
             paste0("LRSomatic Report – ", sample_id)

# ---- Load all available gene panels ----------------------------------------
all_panels = load_all_gene_panels(file.path(repo_dir, "assets"))
default_panel = if (file.exists(file.path(repo_dir, "assets", "gene_lists",
                                            paste0(gene_panel, ".tsv")))) {
  gene_panel
} else {
  if (length(all_panels) > 0) names(all_panels)[1] else "custom"
}

# ---- Auto-detect reference -----------------------------------------------
reference = opt[["reference"]]
if (reference == "auto") {
  # Try VEP somatic VCF first (most likely to have ##contig lines)
  vep_file = Sys.glob(file.path(sample_dir, "vep/somatic/*_SOMATIC_VEP.vcf.gz"))[1]
  sv_file  = Sys.glob(file.path(sample_dir, "variants/severus/somatic_SVs/*.vcf.gz"))[1]
  probe    = if (!is.na(vep_file)) vep_file else if (!is.na(sv_file)) sv_file else NA_character_
  reference = if (!is.na(probe)) detect_reference(probe) else "t2t"
  message("Auto-detected reference: ", reference)
}
reference = tolower(reference)

# ---- Locate per-tool outputs ---------------------------------------------
message("Locating outputs in: ", sample_dir)
outputs = locate_outputs(sample_dir, sample_id)
message("Run mode: ", outputs$mode)
message("VEP somatic: ", ifelse(is.null(outputs$vep_somatic), "NOT FOUND", outputs$vep_somatic))
message("ASCAT segments: ", ifelse(is.null(outputs$ascat_segments), "NOT FOUND", outputs$ascat_segments))

# ---- Render the Quarto template -----------------------------------------
template = file.path(repo_dir, "templates/per_sample.qmd")
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
    gene_panel    = gene_panel,
    default_panel = default_panel,
    all_panels    = all_panels,
    title         = title,
    repo_dir      = repo_dir,
    outputs       = outputs
  ),
  quiet = FALSE
)

# Move output if Quarto wrote it next to the template
rendered = file.path(dirname(template), basename(output))
if (file.exists(rendered)) {
  dest = normalizePath(output, mustWork = FALSE)
  src  = normalizePath(rendered, mustWork = FALSE)
  if (src != dest) {
    ok = file.copy(rendered, output, overwrite = TRUE)
    if (ok) file.remove(rendered)
  }
}

if (file.exists(output)) {
  message("Report written to: ", output)
} else {
  abort("Rendering completed but output file not found at: ", output)
}
