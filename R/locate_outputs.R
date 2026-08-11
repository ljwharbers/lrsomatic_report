# Discover the per-tool output files for a sample run. Discovery is recursive
# under sample_dir: the pipeline may dump inputs flat rather than in a fixed
# directory tree, so files are matched by their distinctive filename suffix.
# Returns a named list; any missing optional file is NULL.

locate_outputs = function(sample_dir, sample_id, mode, somatic_vcf) {
  d = sample_dir  # shorthand

  # First recursive hit under `root` matching a filename pattern
  find1 = function(pattern, root = d) {
    hits = list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (length(hits) > 0) hits[1] else NULL
  }

  # Same, but excluding anything under a normal/ subtree (tumor-side QC)
  find1_tumor = function(pattern) {
    hits = list.files(d, pattern = pattern, recursive = TRUE, full.names = TRUE)
    hits = hits[!grepl("/normal/", hits)]
    if (length(hits) > 0) hits[1] else NULL
  }

  # The mirror of find1_tumor: normal-side files, wherever the pipeline puts them.
  # This has moved (a top-level normal/ historically, qc/normal/ today), so match on
  # the path component rather than rooting the search at a fixed directory.
  find1_normal = function(pattern) {
    hits = list.files(d, pattern = pattern, recursive = TRUE, full.names = TRUE)
    hits = hits[grepl("/normal/", hits)]
    if (length(hits) > 0) hits[1] else NULL
  }

  # --- small variants -------------------------------------------------------
  vep_somatic = find1("_SOMATIC_VEP\\.vcf\\.gz$")

  # The ClairS VCF used for VAF is ambiguous by filename alone (e.g. a generic
  # "somatic.vcf.gz" shared across callers) so it is passed in explicitly rather than
  # discovered; which slot it fills depends on the run mode. It may be several paths,
  # since matched-mode ClairS splits into snvs.vcf.gz + indel.vcf.gz.
  clairs_somatic   = if (mode == "matched")     somatic_vcf else NULL
  clairsto_somatic = if (mode == "tumour-only") somatic_vcf else NULL

  # DeepSomatic, by contrast, is discoverable: its output is named "<sample>.vcf.gz",
  # which only the containing directory disambiguates — so match on the full path.
  deepsomatic_vcf = {
    hits = list.files(d, pattern = "\\.vcf\\.gz$", recursive = TRUE, full.names = TRUE)
    hits = hits[grepl("/deepsomatic/", hits, fixed = TRUE)]
    if (length(hits) > 0) hits[1] else NULL
  }

  # --- structural variants ---------------------------------------------------
  # Severus paths are now located by R/sections/sv.R (section-module contract).

  # --- ASCAT ----------------------------------------------------------------
  ascat_segments_raw = find1("\\.segments_raw\\.txt$")
  ascat_purityploidy = find1("\\.purityploidy\\.txt$")
  ascat_plots = list(
    profile    = find1("\\.tumour\\.ASCATprofile\\.png$"),
    rawprofile = find1("\\.tumour\\.rawprofile\\.png$"),
    sunrise    = find1("\\.tumour\\.sunrise\\.png$"),
    aspcf      = find1("\\.tumour\\.ASPCF\\.png$"),
    before_gc  = find1("\\.before_correction\\..*\\.tumour\\.tumour\\.png$"),
    after_gc   = find1("\\.after_correction_gc.*\\.tumour\\.tumour\\.png$"),
    tumour_sep = find1("^tumorSep.*\\.tumour\\.png$")
  )

  # --- QC (tumor side) --------------------------------------------------------
  mosdepth_summary = find1_tumor("\\.mosdepth\\.summary\\.txt$")
  mosdepth_dist    = find1_tumor("\\.mosdepth\\.global\\.dist\\.txt$")
  cramino_aln      = find1_tumor("_cramino\\.txt$")
  flagstat         = find1_tumor("\\.flagstat$")
  samtools_stats   = find1_tumor("\\.stats$")

  # --- Normal-side QC (matched mode only) -----------------------------------
  normal_mosdepth_summary = find1_normal("\\.mosdepth\\.summary\\.txt$")
  normal_mosdepth_dist    = find1_normal("\\.mosdepth\\.global\\.dist\\.txt$")
  normal_cramino          = find1_normal("_cramino\\.txt$")
  normal_flagstat         = find1_normal("\\.flagstat$")
  normal_samtools_stats   = find1_normal("\\.stats$")

  # Driven by what was actually found rather than by directory layout or run mode:
  # the QC section renders a tumour/normal comparison only if there is normal data.
  has_normal = !is.null(normal_mosdepth_summary) || !is.null(normal_cramino)

  # --- Wakhan (optional) -----------------------------------------------------
  wakhan_dir = file.path(d, "wakhan")
  has_wakhan = dir.exists(wakhan_dir)
  wakhan_solutions = if (has_wakhan) {
    f = file.path(wakhan_dir, "solutions_ranks.tsv")
    if (file.exists(f)) f else NULL
  } else NULL
  wakhan_heatmap = if (has_wakhan) {
    hits = Sys.glob(file.path(wakhan_dir, "*heatmap_ploidy_purity.html"))
    if (length(hits) > 0) hits[1] else NULL
  } else NULL

  list(
    mode             = mode,
    vep_somatic      = vep_somatic,
    clairsto_somatic = clairsto_somatic,
    clairs_somatic   = clairs_somatic,
    deepsomatic_vcf  = deepsomatic_vcf,
    ascat_segments   = ascat_segments_raw,
    ascat_purityploidy = ascat_purityploidy,
    mosdepth_summary = mosdepth_summary,
    mosdepth_dist    = mosdepth_dist,
    cramino          = cramino_aln,
    flagstat         = flagstat,
    samtools_stats   = samtools_stats,
    has_normal       = has_normal,
    ascat_plots      = ascat_plots,
    normal_mosdepth_summary = normal_mosdepth_summary,
    normal_mosdepth_dist    = normal_mosdepth_dist,
    normal_cramino          = normal_cramino,
    normal_flagstat         = normal_flagstat,
    normal_samtools_stats   = normal_samtools_stats,
    has_wakhan       = has_wakhan,
    wakhan_dir       = wakhan_dir,
    wakhan_solutions = wakhan_solutions,
    wakhan_heatmap   = wakhan_heatmap
  )
}
