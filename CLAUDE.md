# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project does

Standalone R/Quarto reporting tool for the [LRSomatic](https://github.com/nf-core/lrsomatic) Nextflow pipeline. Takes a single-sample output directory and generates a self-contained HTML report (circos plot, variant tables, QC summary).

## Running the report

```bash
S=/path/to/CCS15-ONT
Rscript bin/render_report.R \
  --sample-dir  $S \
  --sample-id   CCS15-ONT \
  --sex         male \
  --reference   auto
```

Output defaults to `<sample-id>_report.html` in the current directory. `--reference auto` reads `##contig` lines from the VEP somatic VCF to detect `t2t` vs `hg38`. Matched and tumour-only runs take the same command — every input, and the run mode itself, is discovered from `--sample-dir`.

`--mode` and `--somatic-vcf` were **removed** (they used to be required). Run mode is derived from `has_normal`, and the VCF supplying VAF is discovered by `locate_outputs()`. Anything still passing them fails on an unknown option.

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
| `locate_outputs.R` | `locate_outputs(sample_dir, sample_id)` — recursively discovers all tool output files by filename suffix; derives `mode` and `has_normal` from what it finds |
| `parse_smallvariants.R` | `parse_vep()` dispatches to `parse_vep_text()` (VEP default text format) or `parse_vep_vcf()` (genuine VCF with a CSQ INFO field) based on sniffed file contents; `parse_caller_vcf()` reads VAF/DP/GT/PS from a VCF; `variant_key()` normalises VEP-vs-VCF indel representation; `build_variant_table()` joins the two |
| `parse_severus.R` | Severus VCF + gene annotation TSV parsers |
| `parse_ascat.R` | ASCAT segments and purity/ploidy parsers |
| `parse_qc.R` | mosdepth, cramino, samtools flagstat/stats parsers |
| `circos.R` | `draw_circos()` — generates circos SVG using `circlize`; SBS 6-class colours hard-coded as `SNV_COLOURS` |
| `sections.R` | Section-module contract: `register_section()`, `section_notice()` — see below |
| `sections/*.R` | One file per migrated section (`sv.R`, `whatshap.R`), each calling `register_section()` |

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

One deliberate exception: the SNV section was reworked without being migrated. Its
`vep_somatic` path is needed by `bin/render_report.R` for `--reference auto` *before* Quarto
sources any section module, so moving that discovery into a section module would drag
reference auto-detection into the change. Migrate it only together with a plan for that.

## Key design details

- **VEP file format:** `outputs$vep_somatic` (`*_SOMATIC_VEP.vcf.gz`) ships in either of two incompatible formats depending on the VEP invocation — the filename doesn't tell you which. `parse_vep()` sniffs the header (`#Uploaded_variation` vs `#CHROM`) and dispatches accordingly:
  - `parse_vep_text()` — VEP *default text output* (tab-delimited, `##`-commented header, column header line starts with `#Uploaded_variation`, NOT a VCF). The `Extra` column holds semicolon-delimited `KEY=VALUE` pairs parsed by `parse_extra_kv`.
  - `parse_vep_vcf()` — genuine VCF (`--vcf` VEP output) with annotations in a pipe-delimited `CSQ` INFO field; the field order is read from the `##INFO=<ID=CSQ,...Format: ...>` header line rather than hard-coded.
  Both return the same column contract (`chrom, pos, ref, alt, symbol, gene_id, consequence, impact, hgvsp, existing, dbsnp, cosmic, sift, polyphen, caller`); `derive_dbsnp_cosmic()` is shared between them.
- **A CSQ-format "somatic" VEP VCF may be a merged multi-caller VCF.** This applies to the `parse_vep_vcf()` path only — the pipeline's default `vep_args` produce text output (somatic-only, 27k–167k rows on real samples), and a `--vcf` run is what can carry germline calls (DeepVariant, Clair3) into `*_SOMATIC_VEP.vcf.gz` alongside the somatic ones — on a real sample that's ~99% of the records — distinguished only by an `INFO/CALLER` tag. `parse_vep_vcf()` therefore keeps only `FILTER` ∈ {`PASS`, `.`} records whose caller is in `SOMATIC_CALLERS` (`clairs`, `clairs-to`, `clairsto`, `deepsomatic` — ClairS is tagged `clairs` in matched mode and `clairs-to` in tumour-only). The caller filter is applied **only when a `##INFO=<ID=CALLER` header line is present**, so older single-caller VEP VCFs aren't filtered to nothing. This filter is also what keeps the parse tractable: unfiltered it is ~5.4M variants and ~7 minutes.
- **Missing files are graceful:** Every parser returns `NULL` if its input file is absent; the template shows a "not available" notice per section.
- **Run mode is derived, not declared:** `locate_outputs()` sets `mode = if (has_normal) "matched" else "tumour-only"`. Its only consumer is the hero badge in `_header.qmd`.
- **Small variants use the VEP output as the single source of truth.** The VEP file defines the variant set; VAF/DP/GT/PS are joined from the VCF VEP annotated, which `locate_outputs()` resolves into `somatic_vaf_vcf` by a fallback chain: `variants/phased/somatic_smallvariants.vcf.gz` (what the pipeline actually feeds VEP, and the only rung carrying `FORMAT/PS`) → `variants/clairs{,to}/somatic.vcf.gz` → any non-germline VCF in those directories → `NULL`. Only ~20 of the sample dirs on disk have `variants/phased/`, hence the chain. `--somatic-vcf` and the per-caller `vaf_clairs`/`vaf_clairsto`/`vaf_deepsomatic` columns are gone; the table now has plain `vaf`, `dp`, `gt`, `ps`.
- **VEP and the VCF disagree about indels — `variant_key()` reconciles them.** VEP always reports an indel one base right of the VCF anchor, and depending on version writes either the raw VCF allele pair (`chr1_1871655_TG/T`) or its own trimmed dash form (`chr1_192937_AATA/-`). Both are normalised to *trimmed alleles at anchor+1* before joining. Getting this wrong is silent — it produces `NA` VAF rather than an error — and it was: before this key existed, **0% of deletions joined**. Verified at 100% across all three variant classes on 10 real samples covering both notations. This also requires `parse_vep_text()` to take `chrom`/`pos` from `variant_id`, **not** from `Location`: for a dash-form insertion, `Location`'s start is one base left of the position `variant_id` names, and mixing the two sources is what caused the original breakage.
- **Phasing statistics are germline.** `R/sections/whatshap.R` reads `qc/whatshap_stats/*_whatshap_stats.tsv`, which the pipeline generates from the phased *germline* VCF — every row's `file_name` is `germline_smallvariants.vcf.gz`. The header line starts with `#`, so `setnames(dt, sub("^#", "", names(dt)))` is required after `fread`, and `bp_per_block_sum` arrives as `integer64` and must be widened. `phased_fraction` is a 0–1 fraction. `qc/whatshap_stats/` is *not* `<type>`-scoped, unlike the rest of `qc/`, so a plain recursive match is correct.
- **No methylation section, by decision.** The only methylation output the pipeline publishes is `methylation/<type>/modkit_pileup/<id>.bed.gz` — 38–42 GB gzipped per sample, unfiltered, no tabix index — which cannot be read at render time. `modkit pileup` already runs `--bgzf`, so publishing a `tabix -p bed` index alongside it would unblock this (region queries measured ~0.16 s/Mb on an indexed pileup). The ~1 GB indexed pileups and `*_island_methylation_summary*.tsv` files present under some sample dirs are **ad-hoc analysis outputs, not pipeline outputs** — inconsistent filenames, ambiguous suffix variants, partial coverage — so they are deliberately not wired up. A per-gene methylation view would additionally need a gene-coordinate asset for both references; `assets/references/` has only cytobands and chromosome lengths.
- **Normal-side QC is found by path, not by a fixed root:** `find1_normal()` keeps recursive hits containing `/normal/` (mirroring `find1_tumor()`, which drops them), because the pipeline has moved this from a top-level `normal/` to `qc/normal/`. `has_normal` is derived from what was actually found rather than from `dir.exists()`.
- **Gene-panel filtering is opt-in, and only client-side.** `--gene-panel` defaults to `none` — LRSomatic is a general somatic pipeline, so hiding non-panel variants on load is wrong for most runs. Three-way resolution in `resolve_gene_panel()`: the `none` sentinel returns `NULL`, a builtin name or an existing TSV path loads, and **anything else is an error** so a typo can't silently produce an unfiltered report. `filter_by_gene_panel(dt, NULL)` returns `dt` untouched, while an empty character vector still means a genuinely empty panel — the distinction matters, so don't collapse `NULL` and `character(0)`. Tables are always *built* unfiltered and filtered in the browser; `--gene-panel` only picks which option is selected on load, and a user-supplied TSV is registered alongside the builtins so it can be selected (previously the fallback silently reselected `lymphoid` for any non-builtin value). The `"__all__"` no-filter sentinel is shared by `bin/render_report.R`, `templates/sections/_gene_filter.qmd`, and the search hook in `per_sample.qmd` — change it in all three or none. The "Panel variants"/"Panel SVs" cards read `N/A` when no panel is applied, rather than a `0` that reads as "no panel genes hit".

## R package requirements

```r
install.packages(c("data.table", "dplyr", "tidyr", "DT", "htmltools",
                   "optparse", "quarto", "yaml", "ggplot2", "svglite"))
BiocManager::install(c("circlize", "ComplexHeatmap", "GenomicRanges"))
```

Tested with R 4.4.1 and Quarto 1.5.57.

## Tests

`testthat`, in `tests/testthat/`: `test-utils.R`, `test-severus.R`,
`test-smallvariants.R`. `setup.R` sources the `R/` files under test.

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Run it from the repository root — `setup.R` derives `repo_root` from
`dirname(dirname(getwd()))`, so it relies on testthat setting the working
directory to `tests/testthat`. There is no CI; this and a full render are the
only checks.
