# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project does

Standalone R/Quarto reporting tool for the [LRSomatic](https://github.com/nf-core/lrsomatic) Nextflow pipeline. Takes a single-sample output directory and generates a self-contained HTML report (circos plot, variant tables, QC summary).

## Running the report

```bash
S=/path/to/sample-dir
Rscript bin/render_report.R \
  --sample-dir  $S \
  --sample-id   SAMPLE_ID \
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
| `utils.R` | Shared helpers: gene panel loading (`load_gene_panel` / `resolve_gene_panel` / `panel_intervals`), VEP Extra field parser (`parse_extra_kv` / `extract_extra_key`), `fmt_bp`, `js_col_index_map` |
| `references.R` | Load cytobands/chrom lengths from `assets/`, `detect_reference()`, `chromosomes_for_sex()` |
| `locate_outputs.R` | `locate_outputs(sample_dir, sample_id)` — recursively discovers all tool output files by filename suffix; derives `mode` and `has_normal` from what it finds |
| `parse_smallvariants.R` | `parse_vep()` dispatches to `parse_vep_text()` (VEP default text format) or `parse_vep_vcf()` (genuine VCF with a CSQ INFO field) based on sniffed file contents; `parse_caller_vcf()` reads VAF/DP/GT/PS from a VCF; `variant_key()` normalises VEP-vs-VCF indel representation; `build_variant_table()` joins the two |
| `parse_severus.R` | Severus VCF parsing (`parse_severus_somatic_records` — mate-collapsed, one row per rearrangement), circos tracks (`severus_circos_tracks`), the VEP per-breakend join (`build_sv_table_from_vep`), coordinate panel matching (`sv_panel_hits`), gene-annotation TSV fallback |
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
- **That single `vaf` column is ambiguous after a consensus merge, and the report says so rather than fixing it.** `somatic_vaf_vcf` resolves to *one* file, and the LRSomatic module stages the pipeline's final somatic VCF at exactly that path (which is also what drives run-mode detection). When several somatic callers were combined by consensus, that VCF's FORMAT fields come from whichever caller won each variant — so `vaf`/`dp`/`gt`/`ps` need not correspond to the `callers` value in the same row. The variant **set** is unaffected: it comes from the VEP file, where `INFO/CALLER` is read per record. v1.2.0 makes this visible instead of joining per caller: `vaf_provenance()` (`R/parse_smallvariants.R`) reads just the VCF header for its path and any `##source=` line, and `templates/sections/_smallvariants.qmd` prints a `.table-footnote` naming the file, the join coverage (`n_vaf` of `nrow(variant_table)` — a near-zero count here is how a broken `variant_key()` surfaces), and, when `n_multi > 0`, the consensus caveat plus highlighting on the multi-caller `callers` cells. **`n_multi` is always 0 on the VEP text path**, which carries no per-variant caller; the caveat therefore fires only on the CSQ path, which is exactly where the ambiguity lives. Joining per caller instead would need `modules/local/lrsomaticreport/main.nf` to stage each caller's VCF — a module input-tuple change, so it was deliberately not done here.
- **VEP and the VCF disagree about indels — `variant_key()` reconciles them.** VEP always reports an indel one base right of the VCF anchor, and depending on version writes either the raw VCF allele pair (`chr1_1871655_TG/T`) or its own trimmed dash form (`chr1_192937_AATA/-`). Both are normalised to *trimmed alleles at anchor+1* before joining. Getting this wrong is silent — it produces `NA` VAF rather than an error — and it was: before this key existed, **0% of deletions joined**. Verified at 100% across all three variant classes on 10 real samples covering both notations. This also requires `parse_vep_text()` to take `chrom`/`pos` from `variant_id`, **not** from `Location`: for a dash-form insertion, `Location`'s start is one base left of the position `variant_id` names, and mixing the two sources is what caused the original breakage.
- **One row per rearrangement, not per breakend record.** Severus writes both sides of a rearrangement as separate `_1`/`_2` records linked by `INFO/MATE_ID` (verified 1:1 in every sample checked), and tags single breakends `sBND`, which fails a bare `== "BND"` test. `parse_severus_somatic_records()` collapses mate pairs into one row carrying `chrom_a/pos_a` **and** `chrom_b/pos_b`, derives `svclass` (`translocation` vs `intra-chr breakend` vs `single breakend` vs the plain SV type — 24 of 35 collapsed BNDs in COLO829-ONT are intra-chromosomal, which the old code filed under "translocations"), and is the single VCF read behind the table, the SV count and the circos links, so those three cannot disagree. The partner locus comes from the ALT bracket notation — matched **reference-agnostically**, because the old `chr[^:]+` pattern would drop every BND on a bare-contig reference — and is filled in even for a record whose mate was filtered out. Groups that aren't exactly a pair stay singletons rather than being dropped. Three crashes were latent along this path and are now guarded: `END`/`SVLEN` are created unconditionally (a BND/INS-only VCF has neither), `.info_val()` is length-preserving (bare `regmatches()` drops non-matching elements — an error waiting for the first optional INFO key), and `ensure_chr_prefix()` coerces its input (`fread` types a bare-contig `CHROM` as integer, and `ifelse` on an all-NA test returns logical NA).
- **The SV table's breakend genes are display context; they are not the filter key.** `build_sv_table_from_vep()` joins the VEP SV VCF **on record ID**, with a locus key only as a fallback for annotation files whose IDs don't overlap at all (IDs matched 100% in every sample checked) — the old locus join leaked annotations between distinct records starting at the same base. Both breakends are annotated, giving `gene_a`/`gene_b` in place of the single `gene_hits`; `consequence`/`impact` come from the higher-impact side. The choice of key is made once for the whole table, not per side: per-side would ID-match side A and silently locus-match side B. `parse_severus_somatic_records()` applies no FILTER filter while `parse_vep_vcf()` keeps only `PASS`/`.`, so some rows legitimately have no symbols. A missing VEP SV VCF empties those columns but keeps every SV — panel matching is positional, so the section still filters correctly, and the gene-annotated TSV fallback (`build_sv_table()`) is mapped onto the same schema for the same reason.
- **Phasing statistics are germline.** `R/sections/whatshap.R` reads `qc/whatshap_stats/*_whatshap_stats.tsv`, which the pipeline generates from the phased *germline* VCF — every row's `file_name` is `germline_smallvariants.vcf.gz`. The header line starts with `#`, so `setnames(dt, sub("^#", "", names(dt)))` is required after `fread`, and `bp_per_block_sum` arrives as `integer64` and must be widened. `phased_fraction` is a 0–1 fraction. `qc/whatshap_stats/` is *not* `<type>`-scoped, unlike the rest of `qc/`, so a plain recursive match is correct.
- **No methylation section, by decision.** The only methylation output the pipeline publishes is `methylation/<type>/modkit_pileup/<id>.bed.gz` — 38–42 GB gzipped per sample, unfiltered, no tabix index — which cannot be read at render time. `modkit pileup` already runs `--bgzf`, so publishing a `tabix -p bed` index alongside it would unblock this (region queries measured ~0.16 s/Mb on an indexed pileup). The ~1 GB indexed pileups and `*_island_methylation_summary*.tsv` files present under some sample dirs are **ad-hoc analysis outputs, not pipeline outputs** — inconsistent filenames, ambiguous suffix variants, partial coverage — so they are deliberately not wired up. A per-gene methylation view would additionally need a gene-coordinate asset for both references; `assets/references/` has only cytobands and chromosome lengths.
- **Normal-side QC is found by path, not by a fixed root:** `find1_normal()` keeps recursive hits containing `/normal/` (mirroring `find1_tumor()`, which drops them), because the pipeline has moved this from a top-level `normal/` to `qc/normal/`. `has_normal` is derived from what was actually found rather than from `dir.exists()`.
- **Gene-panel filtering is opt-in, and only client-side.** `--gene-panel` defaults to `none` — LRSomatic is a general somatic pipeline, so hiding non-panel variants on load is wrong for most runs. Three-way resolution in `resolve_gene_panel()`: the `none` sentinel returns `NULL`, a builtin name or an existing TSV path loads, and **anything else is an error** so a typo can't silently produce an unfiltered report. Tables are always *built* unfiltered and filtered in the browser; `--gene-panel` only picks which option is selected on load, and a user-supplied TSV is registered alongside the builtins so it can be selected (previously the fallback silently reselected `lymphoid` for any non-builtin value). The `"__all__"` no-filter sentinel is shared by `bin/render_report.R`, `templates/sections/_gene_filter.qmd`, and the search hook in `per_sample.qmd` — change it in all three or none. The "Panel variants"/"Panel SVs" cards read `N/A` when no panel is applied, rather than a `0` that reads as "no panel genes hit".
- **A panel is a list, not a character vector, and coordinates are what SVs match on.** `load_gene_panel()` returns `list(name, path, reference, has_coords, genes, chrom, start, end, interval_gene)` — a plain list because it round-trips through Quarto's `execute_params` as YAML. Small variants always match on `genes` (symbols); SVs match on **position** when `has_coords`, via `sv_panel_hits()`: a panel interval within `SV_PANEL_WINDOW_BND` (1 Mb) of either breakend of a BND, or `SV_PANEL_WINDOW_OTHER` (100 kb) of the span of any other type. Both windows are defined once, in `R/parse_severus.R`, and emitted to JS from the `panel-js-data` chunk — the R "Panel SVs" card and the browser filter run the same test on the same numbers, which is what stopped them drifting before. The coordinate mode exists because **whether VEP annotates a breakend at all is a property of the sample's VEP invocation** (measured 1.6%–90% of BND rows across samples), so symbol matching hid the very `MYC`/`BCL2`/`BCL6` translocations the `lymphoid` panel is for. A symbol-only panel still works — direct hits on `gene_a`/`gene_b`, no windows — and the SV section footnote names which mode is live, so "symbol-only" can't be misread as "windows found nothing".
- **Panel coordinates are reference-specific, and a mismatch is an error.** A coordinate-carrying TSV must declare its reference (`# reference: hg38` comment line, or a `reference` column); a declaration that disagrees with the rendered reference `stop()`s, and one that declares nothing loads with `reference = ""` and reads "unverified" in the footnote. Symbol-only panels are reference-agnostic. Builtins ship per reference (`lymphoid.hg38.tsv`, `lymphoid.t2t.tsv`) and `load_all_gene_panels(assets, reference)` presents them as one selectable `lymphoid`; a panel that ships only for *other* references is dropped with a message rather than offered. This is why `bin/render_report.R` resolves the reference **before** loading panels — reordering those two blocks back would silently pick the wrong file. Coordinate columns are all-or-nothing: a TSV with `chrom` but no `end` is an error, not a downgrade to symbol mode.

## R package requirements

`recipe/meta.yaml` is the source of truth — it is what the Bioconda package and the
pipeline container are built from, and it must stay in sync with the pipeline's
`modules/local/lrsomaticreport/environment.yml`. Do not add a dependency here without
adding it there.

```r
install.packages(c("data.table", "dplyr", "DT", "htmltools", "optparse",
                   "quarto", "yaml", "ggplot2", "svglite", "knitr",
                   "R.utils", "base64enc"))
BiocManager::install("circlize")
```

`R.utils` is never called directly; `data.table::fread()` needs it for gzipped input.
`svglite`, `base64enc` and `knitr` are used namespaced (`R/circos.R`, `R/utils.R`,
`templates/sections/_qc.qmd`) rather than via `library()`, so grepping for `library(`
under-reports the real dependency set.

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
