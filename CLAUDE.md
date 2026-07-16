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
| `parse_smallvariants.R` | `parse_vep()` dispatches to `parse_vep_text()` (VEP default text format) or `parse_vep_vcf()` (genuine VCF with a CSQ INFO field) based on sniffed file contents; `parse_caller_vcf()` for raw caller VCFs; `build_variant_table()` to join VEP + per-caller VAFs |
| `parse_severus.R` | Severus VCF + gene annotation TSV parsers |
| `parse_ascat.R` | ASCAT segments and purity/ploidy parsers |
| `parse_qc.R` | mosdepth, cramino, samtools flagstat/stats parsers |
| `circos.R` | `draw_circos()` — generates circos SVG using `circlize`; SBS 6-class colours hard-coded as `SNV_COLOURS` |
| `sections.R` | Section-module contract: `register_section()`, `section_notice()` — see below |
| `sections/*.R` | One file per migrated section (e.g. `sv.R`), each calling `register_section()` |

**Assets:** `assets/references/{t2t,hg38}/` contains bundled cytobands and chromosome length TSVs — no network access needed at render time. Gene panels live in `assets/gene_lists/` as TSVs with a `gene` column.

## Section-module contract

Report sections are being migrated (opportunistically, one at a time) to a self-contained
module pattern so a new tool output can be added without touching `locate_outputs.R` or the
setup chunk of `per_sample.qmd`. `R/sections/sv.R` is the reference implementation.

To add a new section:

1. Create `R/sections/<id>.R` calling `register_section(list(...))` with:
   - `id`, `title`
   - `locate(sample_dir, sample_id)` — owns this section's own path discovery (globs/`pick()`),
     returns a named list of inputs (e.g. a `callers` list if multiple tools produce this
     output type, following the pattern in `sv.R`).
   - `parse(inputs, section_data)` — returns one data object for this section, or an
     empty/NULL-ish result if there's nothing to show. `section_data` holds already-parsed
     sections (registration order = parse order), for sections that depend on another (e.g.
     circos reads `SECTION_DATA$sv$circos`).
2. Create `templates/sections/_<id>.qmd` — reads `SECTION_DATA[["<id>"]]`, renders it, and uses
   `section_notice(msg)` (from `R/sections.R`) for the "nothing to show" case instead of a raw
   `tags$div(...)`.
3. Add `{{< include sections/_<id>.qmd >}}` to `templates/per_sample.qmd` in display order.

`R/sections.R` holds the registry (`SECTIONS`, populated at source time) and `section_notice()`.
The setup chunk in `per_sample.qmd` sources every file in `R/sections/`, then runs
`SECTION_DATA[[s$id]] = s$parse(s$locate(sample_dir, sample_id), SECTION_DATA)` for each
registered section. Sections not yet migrated (SNV, ASCAT, QC, circos) keep parsing directly
from `outputs$<key>` in the setup chunk — migrate them the same way when they next need a
change.

## Key design details

- **VEP file format:** `outputs$vep_somatic` (`*_SOMATIC_VEP.vcf.gz`) ships in either of two incompatible formats depending on the VEP invocation — the filename doesn't tell you which. `parse_vep()` sniffs the header (`#Uploaded_variation` vs `#CHROM`) and dispatches accordingly:
  - `parse_vep_text()` — VEP *default text output* (tab-delimited, `##`-commented header, column header line starts with `#Uploaded_variation`, NOT a VCF). The `Extra` column holds semicolon-delimited `KEY=VALUE` pairs parsed by `parse_extra_kv`.
  - `parse_vep_vcf()` — genuine VCF (`--vcf` VEP output) with annotations in a pipe-delimited `CSQ` INFO field; the field order is read from the `##INFO=<ID=CSQ,...Format: ...>` header line rather than hard-coded.
  Both return the same column contract (`chrom, pos, ref, alt, symbol, gene_id, consequence, impact, hgvsp, existing, dbsnp, cosmic, sift, polyphen`); `derive_dbsnp_cosmic()` is shared between them.
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
