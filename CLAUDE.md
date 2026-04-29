# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project does

Standalone R/Quarto reporting tool for the [LRSomatic](https://github.com/nf-core/lrsomatic) Nextflow pipeline. Takes a single-sample output directory and generates a self-contained HTML report (circos plot, variant tables, QC summary).

## Running the report

```bash
Rscript bin/render_report.R \
  --sample-dir /path/to/DLBCL3_pooled \
  --sample-id  DLBCL3_pooled \
  --sex        male \
  --reference  auto
```

Output defaults to `<sample-id>_report.html` in the current directory. `--reference auto` reads `##contig` lines from the VEP somatic VCF to detect `t2t` vs `hg38`.

## R conventions

- Use `=` instead of `<-` for assignments.
- Data manipulation uses `data.table` (primary) and `dplyr` for secondary/display work.
- All R files in `R/` are `source()`d by the Quarto template — they are not packages or modules, just plain scripts.

## Code architecture

**Entry point:** `bin/render_report.R` — parses CLI args, auto-detects reference, calls `locate_outputs()`, then drives `quarto::quarto_render()` on `templates/per_sample.qmd`.

**Template:** `templates/per_sample.qmd` — the Quarto document. All R source files are re-sourced here via `params$repo_dir`. Receives a `params$outputs` named list (from `locate_outputs()`) containing file paths for every tool's output.

**R/ directory — responsibility split:**
| File | Role |
|---|---|
| `utils.R` | Shared helpers: gene panel loading (`resolve_gene_panel`), VEP Extra field parser (`parse_extra_kv` / `extract_extra_key`), `fmt_bp` |
| `references.R` | Load cytobands/chrom lengths from `assets/`, `detect_reference()`, `chromosomes_for_sex()` |
| `locate_outputs.R` | `locate_outputs(sample_dir, sample_id)` — discovers all tool output files and infers run mode (`matched` vs `tumour-only`) |
| `parse_smallvariants.R` | `parse_vep_text()` for VEP default text format (NOT VCF), `parse_caller_vcf()` for raw caller VCFs, `build_variant_table()` to join VEP + per-caller VAFs |
| `parse_severus.R` | Severus VCF + gene annotation TSV parsers |
| `parse_ascat.R` | ASCAT segments and purity/ploidy parsers |
| `parse_qc.R` | mosdepth, cramino, samtools flagstat/stats parsers |
| `circos.R` | `draw_circos()` — generates circos SVG using `circlize`; SBS 6-class colours hard-coded as `SNV_COLOURS` |

**Assets:** `assets/references/{t2t,hg38}/` contains bundled cytobands and chromosome length TSVs — no network access needed at render time. Gene panels live in `assets/gene_lists/` as TSVs with a `gene` column.

## Key design details

- **VEP file format:** `parse_vep_text()` reads VEP *default text output* (tab-delimited, `##`-commented header, column header line starts with `#Uploaded_variation`). This is NOT a VCF. The `Extra` column holds semicolon-delimited `KEY=VALUE` pairs parsed by `parse_extra_kv`.
- **Missing files are graceful:** Every parser returns `NULL` if its input file is absent; the template shows a "not available" notice per section.
- **Run mode detection:** `locate_outputs()` sets `mode = "matched"` if `variants/clairs/` exists, else `"tumour-only"`. This controls which VAF columns appear.
- **Caller VAF join:** `build_variant_table()` joins per-caller VCFs to VEP rows on `chrom|pos|ref|alt`. Caller columns are named `vaf_clairsto`, `vaf_clairs`, `vaf_deepsomatic`.

## R package requirements

```r
install.packages(c("data.table", "dplyr", "tidyr", "DT", "htmltools",
                   "optparse", "quarto", "yaml", "ggplot2", "svglite"))
BiocManager::install(c("circlize", "ComplexHeatmap", "GenomicRanges"))
```

Tested with R 4.4.1 and Quarto 1.5.57.

## Tests

`tests/` directory exists but is empty. Framework planned: `testthat`.
