# Structural variants section. Reference implementation of the section-module
# contract (see R/sections.R and CLAUDE.md). Keyed by caller so a second SV
# caller can be added later without touching the plumbing below.

register_section(list(
  id = "sv",
  title = "Structural variants",

  locate = function(sample_dir, sample_id) {
    d = sample_dir

    pick = function(...) {
      candidates = c(...)
      for (p in candidates) {
        if (!is.null(p) && !is.na(p) && file.exists(p)) return(p)
      }
      NULL
    }
    glob1 = function(pattern) {
      hits = Sys.glob(pattern)
      if (length(hits) > 0) hits[1] else NULL
    }

    severus_vcf = pick(
      file.path(d, "variants/severus/somatic_SVs/severus_somatic.vcf.gz"),
      file.path(d, "variants/severus/severus_somatic.vcf.gz")
    )
    severus_gene_tsv = pick(
      file.path(d, "variants/severus/somatic_SVs/filtered_SV2/SV_filtered_with_gene_annotations.tsv"),
      file.path(d, "variants/severus/filtered_SV2/SV_filtered_with_gene_annotations.tsv")
    )

    list(callers = list(
      severus = list(vcf = severus_vcf, gene_tsv = severus_gene_tsv)
    ))
  },

  parse = function(inputs, section_data) {
    tabs = list()
    circ = list(nontrans = data.table(), translocations = data.table())
    any_gene_tsv = FALSE

    for (nm in names(inputs$callers)) {
      caller_inputs = inputs$callers[[nm]]
      v = parse_severus_vcf(caller_inputs$vcf)
      g = parse_severus_gene_tsv(caller_inputs$gene_tsv)
      if (!is.null(g)) any_gene_tsv = TRUE
      t = build_sv_table(g, gene_panel = NULL)
      if (!is.null(t) && nrow(t) > 0) {
        t[, caller := nm]
        tabs[[nm]] = t
      }
      # Circos tracks are drawn from raw breakpoints, not the gene table;
      # with a single caller today, last-write-wins is a no-op.
      circ$nontrans = v$nontrans
      circ$translocations = v$translocations
    }

    tbl = if (length(tabs) > 0) rbindlist(tabs, fill = TRUE) else data.table()

    list(
      table         = tbl,
      circos        = circ,
      gene_tsv_found = any_gene_tsv
    )
  }
))
