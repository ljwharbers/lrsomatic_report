# Discover the per-tool output files for a sample run and infer the run mode.
# Returns a named list; any missing optional file is NULL.

locate_outputs = function(sample_dir, sample_id) {
  d = sample_dir  # shorthand

  # --- helper ---------------------------------------------------------------
  pick = function(...) {
    candidates = c(...)
    for (p in candidates) {
      if (!is.null(p) && !is.na(p) && file.exists(p)) return(p)
    }
    NULL
  }

  # Glob the first file matching a pattern (for cases where exact name unknown)
  glob1 = function(pattern) {
    hits = Sys.glob(pattern)
    if (length(hits) > 0) hits[1] else NULL
  }

  # --- run mode -------------------------------------------------------------
  has_clairs = dir.exists(file.path(d, "variants/clairs")) &&
    length(Sys.glob(file.path(d, "variants/clairs/*.vcf.gz"))) > 0
  has_clairsto = file.exists(file.path(d, "variants/clairsto/somatic.vcf.gz"))
  has_deepsomatic = length(Sys.glob(file.path(d, "variants/deepsomatic/*.vcf.gz"))) > 0
  has_normal = dir.exists(file.path(d, "qc/normal"))

  mode = if (has_clairs) "matched" else "tumour-only"

  # --- small variants -------------------------------------------------------
  vep_somatic = pick(
    file.path(d, "vep/somatic", paste0(sample_id, "_SOMATIC_VEP.vcf.gz")),
    glob1(file.path(d, "vep/somatic/*_SOMATIC_VEP.vcf.gz"))
  )

  clairsto_somatic = if (has_clairsto) file.path(d, "variants/clairsto/somatic.vcf.gz") else NULL
  clairs_somatic   = if (has_clairs) pick(
    file.path(d, "variants/clairs/somatic.vcf.gz"),
    glob1(file.path(d, "variants/clairs/*somatic*.vcf.gz"))
  ) else NULL
  deepsomatic_vcf  = if (has_deepsomatic) glob1(file.path(d, "variants/deepsomatic/*.vcf.gz")) else NULL

  # --- structural variants --------------------------------------------------
  severus_somatic_vcf = pick(
    file.path(d, "variants/severus/somatic_SVs/severus_somatic.vcf.gz"),
    file.path(d, "variants/severus/severus_somatic.vcf.gz")
  )
  severus_gene_tsv = pick(
    file.path(d, "variants/severus/somatic_SVs/filtered_SV2/SV_filtered_with_gene_annotations.tsv"),
    file.path(d, "variants/severus/filtered_SV2/SV_filtered_with_gene_annotations.tsv")
  )

  # --- ASCAT ----------------------------------------------------------------
  ascat_segments_raw = pick(
    file.path(d, "ascat", paste0(sample_id, ".segments_raw.txt")),
    glob1(file.path(d, "ascat/*.segments_raw.txt"))
  )
  ascat_purityploidy = pick(
    file.path(d, "ascat", paste0(sample_id, ".purityploidy.txt")),
    glob1(file.path(d, "ascat/*.purityploidy.txt"))
  )

  # --- QC -------------------------------------------------------------------
  mosdepth_summary = pick(
    file.path(d, "qc/tumor/mosdepth", paste0(sample_id, ".mosdepth.summary.txt")),
    glob1(file.path(d, "qc/tumor/mosdepth/*.mosdepth.summary.txt"))
  )
  mosdepth_dist = pick(
    file.path(d, "qc/tumor/mosdepth", paste0(sample_id, ".mosdepth.global.dist.txt")),
    glob1(file.path(d, "qc/tumor/mosdepth/*.mosdepth.global.dist.txt"))
  )
  cramino_aln = pick(
    file.path(d, "qc/tumor/cramino_aln", paste0(sample_id, "_cramino.txt")),
    glob1(file.path(d, "qc/tumor/cramino_aln/*_cramino.txt"))
  )
  flagstat = pick(
    file.path(d, "qc/tumor/samtools", paste0(sample_id, ".flagstat")),
    glob1(file.path(d, "qc/tumor/samtools/*.flagstat"))
  )
  samtools_stats = pick(
    file.path(d, "qc/tumor/samtools", paste0(sample_id, ".stats")),
    glob1(file.path(d, "qc/tumor/samtools/*.stats"))
  )

  # --- Wakhan (optional v2) -------------------------------------------------
  has_wakhan = dir.exists(file.path(d, "wakhan"))

  list(
    mode             = mode,
    vep_somatic      = vep_somatic,
    clairsto_somatic = clairsto_somatic,
    clairs_somatic   = clairs_somatic,
    deepsomatic_vcf  = deepsomatic_vcf,
    severus_vcf      = severus_somatic_vcf,
    severus_gene_tsv = severus_gene_tsv,
    ascat_segments   = ascat_segments_raw,
    ascat_purityploidy = ascat_purityploidy,
    mosdepth_summary = mosdepth_summary,
    mosdepth_dist    = mosdepth_dist,
    cramino          = cramino_aln,
    flagstat         = flagstat,
    samtools_stats   = samtools_stats,
    has_wakhan       = has_wakhan
  )
}
