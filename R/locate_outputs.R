# Locate per-tool output files under sample_dir by filename suffix (recursive); missing optional files are NULL

locate_outputs = function(sample_dir, sample_id) {
  d = sample_dir  # shorthand

  # First recursive hit under `root` matching a filename pattern
  find1 = function(pattern, root = d) {
    hits = list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (length(hits) > 0) hits[1] else NULL
  }

  # QC files: one recursive listing, relative paths with a leading "/" so the "/normal/" side test and
  # the path preferences below see only the part under sample_dir, never sample_dir itself
  qc_rel = NULL
  qc_files = function() {
    if (is.null(qc_rel)) qc_rel <<- paste0("/", list.files(d, recursive = TRUE))
    qc_rel
  }

  # Pick one QC file for `side` whose BASENAME matches `pattern`. Several can match (cramino runs on the
  # aligned BAM and on each uBAM, all writing <id>_<type>_cramino.txt), so hits are ranked rather than
  # taken alphabetically: under /qc/ first, then by `prefer` (regexes, in order), then away from `avoid`.
  # `within`, when set, is a hard path requirement for loose basename patterns. The pipeline-staged
  # layout links files flat into qc/<side>/ with no tool subdirectories, so preferences only re-rank.
  pick_qc = function(pattern, side = c("tumor", "normal"), prefer = character(0), avoid = NULL, within = NULL) {
    side = match.arg(side)
    hits = qc_files()
    hits = hits[grepl(pattern, basename(hits))]
    is_normal = grepl("/normal/", hits, fixed = TRUE)
    hits = if (side == "normal") hits[is_normal] else hits[!is_normal]
    if (!is.null(within)) hits = hits[grepl(within, hits)]
    if (length(hits) == 0) return(NULL)
    keys = c(list(!grepl("/qc/", hits, fixed = TRUE)),
             lapply(prefer, function(p) !grepl(p, hits)),
             if (!is.null(avoid)) list(grepl(avoid, hits)),
             list(hits))
    file.path(d, substring(hits[do.call(order, keys)[1]], 2))
  }

  # mosdepth: "<prefix>.mosdepth.summary.txt" (older pipeline) or "<id>_<type>.summary.txt" (current,
  # prefix carries the type). The bare ".summary.txt" form is generic, so it must sit under qc/ or mosdepth/.
  pick_mosdepth_summary = function(side) {
    hit = pick_qc("\\.mosdepth\\.summary\\.txt$", side)
    if (is.null(hit)) hit = pick_qc("\\.summary\\.txt$", side, prefer = "/mosdepth/", within = "/(qc|mosdepth)/")
    hit
  }
  # ".global.dist.txt" is mosdepth's alone, with or without the ".mosdepth" infix
  pick_mosdepth_dist = function(side) pick_qc("\\.global\\.dist\\.txt$", side, prefer = "\\.mosdepth\\.global")
  # Aligned-BAM cramino is what the report documents; cramino_ubam* holds the raw-read twin
  pick_cramino = function(side) pick_qc("_cramino\\.txt$", side, prefer = "/cramino_aln/", avoid = "/cramino_ubam")

  # --- small variants -------------------------------------------------------
  vep_somatic = find1("_SOMATIC_VEP\\.vcf\\.gz$")

  # VAF/depth/phasing come from the VCF VEP annotated: prefer the phased somatic VCF, fall back to raw ClairS(-TO) output (may be several paths)
  somatic_vaf_vcf = {
    phased = file.path(d, "variants", "phased", "somatic_smallvariants.vcf.gz")
    if (file.exists(phased)) phased else {
      vcfs = list.files(d, pattern = "\\.vcf\\.gz$", recursive = TRUE, full.names = TRUE)
      named = vcfs[grepl("/clairs(to)?/somatic\\.vcf\\.gz$", vcfs)]
      any_c = vcfs[grepl("/clairs(to)?/", vcfs) & !grepl("/germline\\.vcf\\.gz$", vcfs)]
      if (length(named) > 0) named[1] else if (length(any_c) > 0) any_c else NULL
    }
  }

  # ClairS-TO Verdict germline tagging status (pipeline contract: variants/clairsto/<id>.verdict_status.tsv)
  verdict_status = {
    hits = list.files(d, pattern = "\\.verdict_status\\.tsv$", recursive = TRUE, full.names = TRUE)
    own = hits[basename(hits) == paste0(sample_id, ".verdict_status.tsv")]
    if (length(own) > 0) own[1] else if (length(hits) > 0) hits[1] else NULL
  }

  # --- structural variants: located by R/sections/sv.R ---

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
  mosdepth_summary = pick_mosdepth_summary("tumor")
  mosdepth_dist    = pick_mosdepth_dist("tumor")
  cramino_aln      = pick_cramino("tumor")
  flagstat         = pick_qc("\\.flagstat$", "tumor", prefer = "/samtools/")
  samtools_stats   = pick_qc("\\.stats$", "tumor", prefer = "/samtools/")

  # --- Normal-side QC (matched mode only) -----------------------------------
  normal_mosdepth_summary = pick_mosdepth_summary("normal")
  normal_mosdepth_dist    = pick_mosdepth_dist("normal")
  normal_cramino          = pick_cramino("normal")
  normal_flagstat         = pick_qc("\\.flagstat$", "normal", prefer = "/samtools/")
  normal_samtools_stats   = pick_qc("\\.stats$", "normal", prefer = "/samtools/")

  # Driven by what was found: the QC comparison renders only if there is normal data
  # (any normal-side QC file counts: a run can skip mosdepth and cramino but keep samtools stats)
  has_normal = !is.null(normal_mosdepth_summary) || !is.null(normal_cramino) ||
    !is.null(normal_flagstat) || !is.null(normal_samtools_stats)

  # Run mode derived from the same evidence; --run-mode overrides it (apply_run_mode())
  mode = if (has_normal) "matched" else "tumour-only"

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
    mode_source      = "inferred",
    vep_somatic      = vep_somatic,
    somatic_vaf_vcf  = somatic_vaf_vcf,
    verdict_status   = verdict_status,
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

# --run-mode values (the pipeline's contract) -> the report's internal mode names
RUN_MODES = c(paired = "matched", tumour_only = "tumour-only")

# Override the inferred run mode with --run-mode; NULL/"" keeps the inference. has_normal is
# left alone: it says whether normal-side QC files exist, which is what the QC section needs
apply_run_mode = function(outputs, run_mode = NULL) {
  if (is.null(run_mode) || length(run_mode) == 0 || is.na(run_mode) || !nzchar(trimws(run_mode)))
    return(outputs)
  rm = tolower(trimws(run_mode))
  if (!rm %in% names(RUN_MODES))
    stop("--run-mode must be 'paired' or 'tumour_only' (got '", run_mode, "')")
  outputs$mode = unname(RUN_MODES[[rm]])
  outputs$mode_source = "cli"
  outputs
}
