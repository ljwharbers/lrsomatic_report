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
| `utils.R` | Shared helpers: gene panel loading (`load_gene_panel` / `resolve_gene_panel` / `panel_intervals`), VEP Extra field parser (`parse_extra_kv` / `extract_extra_key`), `fmt_bp`, JS serialisation (`js_col_index_map` / `js_vec` / `js_rows`) |
| `references.R` | Load cytobands/chrom lengths from `assets/`, `detect_reference()`, `chromosomes_for_sex()` |
| `locate_outputs.R` | `locate_outputs(sample_dir, sample_id)` — recursively discovers all tool output files by filename suffix; derives `mode` and `has_normal` from what it finds |
| `parse_smallvariants.R` | `parse_vep()` dispatches to `parse_vep_text()` (VEP default text format) or `parse_vep_vcf()` (genuine VCF with a CSQ INFO field) based on sniffed file contents; `parse_caller_vcf()` reads VAF/DP/GT/PS from a VCF; `variant_key()` normalises VEP-vs-VCF indel representation; `build_variant_table()` joins the two |
| `parse_severus.R` | Severus VCF parsing (`parse_severus_somatic_records` — mate-collapsed, one row per rearrangement), circos tracks (`severus_circos_tracks`), the VEP per-breakend join (`build_sv_table_from_vep`), coordinate panel matching (`sv_panel_hits`), gene-annotation TSV fallback |
| `parse_ascat.R` | ASCAT segments and purity/ploidy parsers |
| `parse_qc.R` | mosdepth, cramino, samtools flagstat/stats parsers |
| `circos.R` | `draw_circos()` — generates the genome-wide circos SVG using `circlize`; SBS 6-class colours hard-coded as `SNV_COLOURS` |
| `circos_bnd.R` | The breakend circos in the SV section, R half only: `bnd_links()` picks the drawable rearrangements, `bnd_panel_genes()` the panel genes to label, `bnd_circos_data()` serialises both into `window.BND_DATA`. The plot itself is `assets/js/bnd_circos.js` |
| `sections.R` | Section-module contract: `register_section()`, `section_notice()`, `table_details()` — see below |
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
3. Add `{{< include sections/_<id>.qmd >}}` to `templates/per_sample.qmd` in display order —
   or, for a subsection of an existing section, to that section's own `.qmd` — `_qc.qmd` includes
   `_whatshap.qmd` that way. Nested include paths are still resolved relative to the **root**
   document (`templates/per_sample.qmd`), not to the including file, so a nested include keeps
   the same `sections/_<id>.qmd` prefix a top-level one has.

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
  Both return the same column contract (`chrom, pos, ref, alt, symbol, gene_id, consequence, impact, hgvsp, existing, dbsnp, cosmic, sift, polyphen, caller, coord_space`); `derive_dbsnp_cosmic()` is shared between them.
- **A CSQ-format "somatic" VEP VCF may be a merged multi-caller VCF.** This applies to the `parse_vep_vcf()` path only — the pipeline's default `vep_args` produce text output (somatic-only, 27k–167k rows on real samples), and a `--vcf` run is what can carry germline calls (DeepVariant, Clair3) into `*_SOMATIC_VEP.vcf.gz` alongside the somatic ones — on a real sample that's ~99% of the records — distinguished only by an `INFO/CALLER` tag. `parse_vep_vcf()` therefore keeps only `FILTER` ∈ {`PASS`, `.`} records whose caller is in `SOMATIC_CALLERS` (`clairs`, `clairs-to`, `clairsto`, `deepsomatic` — ClairS is tagged `clairs` in matched mode and `clairs-to` in tumour-only). The caller filter is applied **only when a `##INFO=<ID=CALLER` header line is present**, so older single-caller VEP VCFs aren't filtered to nothing. This filter is also what keeps the parse tractable: unfiltered it is ~5.4M variants and ~7 minutes.
- **Missing files are graceful:** Every parser returns `NULL` if its input file is absent; the template shows a "not available" notice per section.
- **Run mode is derived, not declared:** `locate_outputs()` sets `mode = if (has_normal) "matched" else "tumour-only"`. Its only consumer is the hero badge in `_header.qmd`.
- **Small variants use the VEP output as the single source of truth.** The VEP file defines the variant set; VAF/DP/GT/PS are joined from the VCF VEP annotated, which `locate_outputs()` resolves into `somatic_vaf_vcf` by a fallback chain: `variants/phased/somatic_smallvariants.vcf.gz` (what the pipeline actually feeds VEP, and the only rung carrying `FORMAT/PS`) → `variants/clairs{,to}/somatic.vcf.gz` → any non-germline VCF in those directories → `NULL`. Only ~20 of the sample dirs on disk have `variants/phased/`, hence the chain. `--somatic-vcf` and the per-caller `vaf_clairs`/`vaf_clairsto`/`vaf_deepsomatic` columns are gone; the table now has plain `vaf`, `dp`, `gt`, `ps`.
- **That single `vaf` column is ambiguous after a consensus merge, and the report says so rather than fixing it.** `somatic_vaf_vcf` resolves to *one* file, and the LRSomatic module stages the pipeline's final somatic VCF at exactly that path (which is also what drives run-mode detection). When several somatic callers were combined by consensus, that VCF's FORMAT fields come from whichever caller won each variant — so `vaf`/`dp`/`gt`/`ps` need not correspond to the `callers` value in the same row. The variant **set** is unaffected: it comes from the VEP file, where `INFO/CALLER` is read per record. v1.2.0 makes this visible instead of joining per caller: `vaf_provenance()` (`R/parse_smallvariants.R`) reads just the VCF header for its path and any `##source=` line, and `templates/sections/_smallvariants.qmd` prints a `.table-footnote` naming the file, the join coverage (`n_vaf` of `nrow(variant_table)` — a near-zero count here is how a broken `variant_key()` surfaces), and, when `n_multi > 0`, the consensus caveat plus highlighting on the multi-caller `callers` cells. **`n_multi` is always 0 on the VEP text path**, which carries no per-variant caller; the caveat therefore fires only on the CSQ path, which is exactly where the ambiguity lives. Joining per caller instead would need `modules/local/lrsomaticreport/main.nf` to stage each caller's VCF — a module input-tuple change, so it was deliberately not done here.
- **VEP and the VCF disagree about indels — `variant_key()` reconciles them.** VEP always reports an indel one base right of the VCF anchor, and depending on version writes either the raw VCF allele pair (`chr1_1871655_TG/T`) or its own trimmed dash form (`chr1_192937_AATA/-`). Both are normalised to *trimmed alleles at anchor+1* before joining. Getting this wrong is silent — it produces `NA` VAF rather than an error — and it was: before this key existed, **0% of deletions joined**. Verified at 100% across all three variant classes on 10 real samples covering both notations. This also requires `parse_vep_text()` to take `chrom`/`pos` from `variant_id`, **not** from `Location`: for a dash-form insertion, `Location`'s start is one base left of the position `variant_id` names, and mixing the two sources is what caused the original breakage.
- **Which coordinate space a VEP table is in is declared by its parser, not assumed by the joiner.** The two parsers genuinely differ: `parse_vep_text()` reads coordinates from `variant_id` (VEP space — indels already shifted, possibly dash-form), while `parse_vep_vcf()` reads the VCF's own `POS`/`REF`/`ALT` (VCF space — the *original* input coordinate, which VEP does not rewrite in `--vcf` mode). Each sets a `coord_space` column and `build_variant_table()` keys `variant_key()` from it. Hard-coding `space = "vep"` there shifted only one side of the key on the CSQ path, so **no indel could ever join** and its `vaf`/`dp`/`gt`/`ps` came back silently `NA` — measured on COLO829-ONT as 717 of 49,003 rows, exactly the 395 deletions + 322 insertions. Now 100% on all three classes. `parse_vep_vcf()` deliberately does **not** normalise into VEP space: that would change the `pos`/`ref`/`alt` a reader *sees* for every CSQ-path indel. Only the key changes.
- **The allele-fraction FORMAT tag is not always `AF`.** LRSomatic renames it by prioritised caller — `STANDARDIZE_AF` in `subworkflows/local/small_variant_consensus.nf` maps `AF → VAF` for `deepvariant`/`deepsomatic` and `VAF → AF` for `clair`, and a `consensus` merge skips the renaming entirely — so both names reach `parse_caller_vcf()`. It hard-coded `AF`, which is why a DeepSomatic run rendered an **entirely empty VAF column** while `dp`/`gt` beside it were fine (same merge, same key — the join was never the problem). `allele_fraction()` now tries `AF`, then `VAF`, then derives from `AD`. The pipeline default is `prioritize_caller_somatic = 'clair'`, i.e. `AF`, which is why this went unnoticed.
- **The sample column is chosen by name.** `parse_caller_vcf()` read column 10 positionally; on a matched two-sample VCF where NORMAL sorts first that reports the **normal's** VAF/DP for every somatic variant, with no error and nothing in the report to reveal it. `pick_sample_column()` prefers the sample named for the report, then the first not named like a normal, and warns when it has to guess. `vaf_provenance()` names the chosen sample in the footnote — but only when the VCF held more than one, since that is the only case where the choice can be wrong.
- **One row per rearrangement, not per breakend record.** Severus writes both sides of a rearrangement as separate `_1`/`_2` records linked by `INFO/MATE_ID` (verified 1:1 in every sample checked), and tags single breakends `sBND`, which fails a bare `== "BND"` test. `parse_severus_somatic_records()` collapses mate pairs into one row carrying `chrom_a/pos_a` **and** `chrom_b/pos_b`, derives `svclass` (`translocation` vs `intra-chr breakend` vs `single breakend` vs the plain SV type — 24 of 35 collapsed BNDs in COLO829-ONT are intra-chromosomal, which the old code filed under "translocations"), and is the single VCF read behind the table, the SV count and the circos links, so those three cannot disagree. The partner locus comes from the ALT bracket notation — matched **reference-agnostically**, because the old `chr[^:]+` pattern would drop every BND on a bare-contig reference — and is filled in even for a record whose mate was filtered out. Groups that aren't exactly a pair stay singletons rather than being dropped. Three crashes were latent along this path and are now guarded: `END`/`SVLEN` are created unconditionally (a BND/INS-only VCF has neither), `.info_val()` is length-preserving (bare `regmatches()` drops non-matching elements — an error waiting for the first optional INFO key), and `ensure_chr_prefix()` coerces its input (`fread` types a bare-contig `CHROM` as integer, and `ifelse` on an all-NA test returns logical NA).
- **The SV table's breakend genes are display context; they are not the filter key.** `build_sv_table_from_vep()` joins the VEP SV VCF **on record ID**, with a locus key only as a fallback for annotation files whose IDs don't overlap at all (IDs matched 100% in every sample checked) — the old locus join leaked annotations between distinct records starting at the same base. Both breakends are annotated, giving `gene_a`/`gene_b` in place of the single `gene_hits`; `consequence`/`impact` come from the higher-impact side. The choice of key is made once for the whole table, not per side: per-side would ID-match side A and silently locus-match side B. `parse_severus_somatic_records()` applies no FILTER filter while `parse_vep_vcf()` keeps only `PASS`/`.`, so some rows legitimately have no symbols. A missing VEP SV VCF empties those columns but keeps every SV — panel matching is positional, so the section still filters correctly, and the gene-annotated TSV fallback (`build_sv_table()`) is mapped onto the same schema for the same reason.
- **Phasing statistics are germline.** `R/sections/whatshap.R` reads `qc/whatshap_stats/*_whatshap_stats.tsv`, which the pipeline generates from the phased *germline* VCF — every row's `file_name` is `germline_smallvariants.vcf.gz`. The header line starts with `#`, so `setnames(dt, sub("^#", "", names(dt)))` is required after `fread`, and `bp_per_block_sum` arrives as `integer64` and must be widened. `phased_fraction` is a 0–1 fraction. `qc/whatshap_stats/` is *not* `<type>`-scoped, unlike the rest of `qc/`, so a plain recursive match is correct. The section renders as a collapsible **Phasing** block inside **QC details** rather than as a top-level section: `_whatshap.qmd` carries its own `::: {.callout-note collapse="true"}` fence and an `###` heading, and `_qc.qmd` includes it as its last block, so it matches its three sibling QC blocks. It is still a registered section module, so `SECTION_DATA$whatshap` — which the header's "Phased variants" card reads — is unaffected by where it is drawn.
- **No methylation section, by decision.** The only methylation output the pipeline publishes is `methylation/<type>/modkit_pileup/<id>.bed.gz` — 38–42 GB gzipped per sample, unfiltered, no tabix index — which cannot be read at render time. `modkit pileup` already runs `--bgzf`, so publishing a `tabix -p bed` index alongside it would unblock this (region queries measured ~0.16 s/Mb on an indexed pileup). The ~1 GB indexed pileups and `*_island_methylation_summary*.tsv` files present under some sample dirs are **ad-hoc analysis outputs, not pipeline outputs** — inconsistent filenames, ambiguous suffix variants, partial coverage — so they are deliberately not wired up. A per-gene methylation view would additionally need a gene-coordinate asset for both references; `assets/references/` has only cytobands and chromosome lengths.
- **Normal-side QC is found by path, not by a fixed root:** `find1_normal()` keeps recursive hits containing `/normal/` (mirroring `find1_tumor()`, which drops them), because the pipeline has moved this from a top-level `normal/` to `qc/normal/`. `has_normal` is derived from what was actually found rather than from `dir.exists()`.
- **Gene-panel filtering is opt-in, and only client-side.** `--gene-panel` defaults to `none` — LRSomatic is a general somatic pipeline, so hiding non-panel variants on load is wrong for most runs. Three-way resolution in `resolve_gene_panel()`: the `none` sentinel returns `NULL`, a builtin name or an existing TSV path loads, and **anything else is an error** so a typo can't silently produce an unfiltered report. Tables are always *built* unfiltered and filtered in the browser; `--gene-panel` only picks which option is selected on load, and a user-supplied TSV is registered alongside the builtins so it can be selected (previously the fallback silently reselected `lymphoid` for any non-builtin value). The `"__all__"` no-filter sentinel is shared by `bin/render_report.R`, `templates/sections/_gene_filter.qmd`, and the search hook in `per_sample.qmd` — change it in all three or none. The "Panel variants"/"Panel SVs" cards read `N/A` when no panel is applied, rather than a `0` that reads as "no panel genes hit".
- **A panel is a list, not a character vector, and coordinates are what SVs match on.** `load_gene_panel()` returns `list(name, path, reference, has_coords, genes, chrom, start, end, interval_gene)` — a plain list because it round-trips through Quarto's `execute_params` as YAML. Small variants always match on `genes` (symbols); SVs match on **position** when `has_coords`, via `sv_panel_hits()`: a panel interval within `SV_PANEL_WINDOW_BND` (1 Mb) of either breakend of a BND, or `SV_PANEL_WINDOW_OTHER` (100 kb) of the span of any other type. Both windows are defined once, in `R/parse_severus.R`, and emitted to JS from the `panel-js-data` chunk — the R "Panel SVs" card and the browser filter run the same test on the same numbers, which is what stopped them drifting before. The coordinate mode exists because **whether VEP annotates a breakend at all is a property of the sample's VEP invocation** (measured 1.6%–90% of BND rows across samples), so symbol matching hid the very `MYC`/`BCL2`/`BCL6` translocations the `lymphoid` panel is for. A symbol-only panel still works — direct hits on `gene_a`/`gene_b`, no windows — and the SV section footnote names which mode is live, so "symbol-only" can't be misread as "windows found nothing".
- **A `panel_hit` entry says how it matched, because the windows are wide enough that it has to.** Each reads `GENE (side, how)`: `side` is `A`/`B`/`span`, and `how` is `direct` when the locus overlaps the gene interval or `fmt_bp()` of the gap when it merely fell inside the window. The label originally carried only the side, and on COLO829-ONT that rendered a breakend 188.5 kb past `RB1` identically to two deletions that swallow `RB1`, `GPC5` and `PTEN` whole — proximity reading as disruption. The distance is measured from the **unpadded** locus, so `sv_panel_hits()` carries `qlo`/`qhi` through `foverlaps()` alongside the padded `start`/`end` it matches on; measuring from the window edge would report zero for everything. The symbol path is always `direct` — a VEP symbol sits on the breakend, and there is nothing to measure. As with the windows themselves, the label is built **twice** and the table's copy is the client-side one (`svPanelHits()` in the `panel-js-data` chunk fills the cells; R's copy only feeds the "Panel SVs" card), so both must change together. JS mirrors `fmt_bp()` as `Math.round(x * 10) / 10`, not `toFixed(1)` — R prints the rounded number, so `toFixed` would render `188.5 kb` as `188.5 kb` but `306 kb` as `306.0 kb`. Verified identical on all 154 SV rows of COLO829-ONT.
- **Panel coordinates are reference-specific, and a mismatch is an error.** A coordinate-carrying TSV must declare its reference (`# reference: hg38` comment line, or a `reference` column); a declaration that disagrees with the rendered reference `stop()`s, and one that declares nothing loads with `reference = ""` and reads "unverified" in the footnote. Symbol-only panels are reference-agnostic. Builtins ship per reference (`lymphoid.hg38.tsv`, `lymphoid.t2t.tsv`) and `load_all_gene_panels(assets, reference)` presents them as one selectable `lymphoid`; a panel that ships only for *other* references is dropped with a message rather than offered. This is why `bin/render_report.R` resolves the reference **before** loading panels — reordering those two blocks back would silently pick the wrong file. Coordinate columns are all-or-nothing: a TSV with `chrom` but no `end` is an error, not a downgrade to symbol mode.
- **The breakend circos is drawn in the browser, and that is what lets it re-lay-out.** `templates/sections/_sv.qmd` shows a second, BND-only circos over just the chromosomes a breakend touches, and clicking a row in the SV table highlights that row's arc. R selects the records (`bnd_links()`, `bnd_panel_genes()`) and serialises them (`bnd_circos_data()` → `window.BND_DATA`); `assets/js/bnd_circos.js` builds the SVG and **rebuilds it from scratch on every filter change**. Three things about that split:
  - **It has to be client-side, because the layout is a function of the filter.** Which chromosomes get a sector — and therefore every other sector's angular width, since they share the circle — depends on which links survive. A server-drawn SVG bakes those angles in, and circlize's arcs are flattened `<polyline>` point lists carrying no genomic coordinates, so no script can recompute them. The previous version could only *dim* arcs in place, leaving every chromosome on the plot however the panel was set. Pre-rendering one SVG per panel was the alternative and was rejected: it cannot serve `Custom…` or the per-column search boxes, which produce row subsets no pre-render anticipates.
  - **This retired the sentinel-colour tagging, and its whole failure mode with it.** The old design inlined svglite output, which emits no ids or classes, so every arc, label and gene body was *drawn* in a unique colour and the markup rewritten afterwards to attach `class`/`data-*`; a mistag downgraded the plot to a static figure, because highlighting the wrong arc is worse than highlighting none. Nodes built by `createElementNS()` carry their own attributes by construction, so there is nothing left to mis-tag. Gone with it: the `.bnd-svglite` CSS-leak rename, the `viewBox`-case repair (Quarto lowercases attributes on *markup* it processes, not on script-created nodes), and the `interactive = FALSE` fallback.
  - **The payload is row-major arrays, hand-serialised.** `js_vec()`/`js_rows()` in `R/utils.R` emit `[[...],[...]]` rather than arrays of objects, and cytobands are trimmed to the plotted sectors — 6 kB for a 4-arc plot against 70–160 kB for the SVG it replaced. Hand-rolled because `jsonlite` would be a new dependency in both `recipe/meta.yaml` and the pipeline's `environment.yml`. `js_num()` formats element-wise on purpose: vectorised `format()` pads to a common format (`c(1, 2.5)` → `"1.0"`), and scientific notation on a base-pair coordinate reads back as a different number.

  Two rules about *what* it draws. Arcs are `svclass ∈ {translocation, intra-chr breakend}` with both loci mapped — intra-chromosomal breakends are included so that every BND row a user can click has an arc, and a single breakend has no partner locus so it necessarily has none. **The gene track is the `panel_hit` column's genes**: R ships the superset (every coordinate-carrying panel's genes within `SV_PANEL_WINDOW_BND` of a drawn breakend, via the same `foverlaps()` test as `sv_panel_hits()`), and the client draws only the genes named by currently-visible rows. Two consequences that look like bugs but are not: under the default `--gene-panel none` no row has a panel hit, so **the gene track is empty on load**; and a symbol-only or custom panel produces `panel_hit` labels but no gene bodies, because it carries no coordinates. Gene bodies have a minimum angular width — MYC is 7 kb against a 145 Mb sector — so they are position markers, not spans to scale; their labels are de-overlapped by a two-pass greedy sweep, with a connector line showing the displacement. Selection state lives in the section's own script rather than DT's (`selection = "none"`): `Scroller` + `deferRender` destroy `<tr>` nodes on scroll, and it is re-applied from a `Set` of SV ids on every `draw.dt`.

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
`test-smallvariants.R`, `test-circos-bnd.R`. `setup.R` sources the `R/` files under
test; it sources `R/circos_bnd.R` only when `circlize` is installed, and the tests that
need it `skip_if_not_installed()`, so the rest of the suite still runs without it.

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Run it from the repository root — `setup.R` derives `repo_root` from
`dirname(dirname(getwd()))`, so it relies on testthat setting the working
directory to `tests/testthat`. There is no CI; this and a full render are the
only checks.
