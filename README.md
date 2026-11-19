# lrsomatic_report

Standalone reporting tool for the [LRSomatic](https://github.com/nf-core/lrsomatic) Nextflow
pipeline. It turns one sample's pipeline output directory into a single self-contained HTML
report, plus the full small-variant table as a TSV beside it.

> **Research use only.** This report is not a validated diagnostic product and must not be
> used on its own for clinical decisions. Every call it shows comes from the pipeline's tools,
> filtered and summarised by the rules documented below; confirm anything actionable with a
> validated method. The printed report carries the same statement in its footer.

## What the report contains

In page order:

- **Overview** — the sample band (run mode, reference, sex), the gene-list chips, and:
  - **Key findings**: alterations in the genes of the ticked gene list(s) that carry evidence
    of being drivers — small variants, structural variants and copy-number changes — each with
    a tier (see [Key findings](#key-findings)). **Key findings need a gene list**: with none
    ticked, the table is replaced by a prompt listing how many tier 1–2 findings each embedded
    list has.
  - **Genome-wide biomarkers**: purity / ploidy, whole-genome duplication, LOH fraction, coding
    TMB, small-variant counts, SV counts per caller, chromothripsis, ecDNA / focal amplicons
    and dominant mutational signatures. A biomarker whose tool did not run reads "not run".
  - **Sequencing** cards: mean coverage, read N50, error rate, phased germline variants.
- **Circos overview** — ideogram, somatic **SNV density** (per 1 Mb window), non-breakend SVs
  and ASCAT copy number, with translocation links, and the tier 1–2 key-finding genes of the
  ticked gene list(s) labelled around the outside.
- **Copy number profile** — ASCAT plots; a static haplotype-resolved **Wakhan** profile of the
  rank-1 solution with the solutions table; and **gene-level copy number** for every gene of the
  embedded gene lists (see [Copy-number calls](#copy-number-calls)).
- **Genome instability** — chromothripsis (ShatterSeek, per SV/CN caller combination) and
  ecDNA / focal amplicons (CoRAL, AmpliconClassifier). A negative result is one line: no
  chromothripsis in any combination is a single sentence naming them all, and calls are one
  table with a column for the combination; amplicons of which none is ecDNA or BFB are listed
  in a collapsed block under the verdict. Each call names the PCAWG tier(s) it meets (HC1,
  HC2, HC3), and a minimum-tier selector (all tiers, HC1–HC2, HC1 only), shared with the
  overview's chromothripsis card, narrows the calls shown; it opens on all tiers.
- **Mutational signatures** — SigProfilerAssignment COSMIC fits for SBS96, DBS78 and ID83, with
  the mutational profile and a proposed aetiology per signature (see
  [Mutational signatures](#mutational-signatures)).
- **Small variants** — the VEP-annotated somatic table (VAF, depth, phasing, plugin scores,
  ClinVar), with tickbox column filters and the gene-list filter. By default only the relevant
  rows are embedded (see [Report size](#report-size-and---embed)); the full table is always
  written next to the report as `<sample-id>_small_variants.tsv.gz`.
- **Structural variants** — one row per rearrangement from **Severus** and, where it ran,
  **SAVANA** (listed side by side, never merged), with VEP breakend genes, **Padfoot**
  annotations, a breakend circos that follows the table filter, and **ReConPlot** figures.
- **QC details** — mosdepth coverage, samtools flagstat, cramino read statistics and WhatsHap
  phasing statistics (germline).

## Quick start

```bash
S=/path/to/sample-dir

Rscript bin/render_report.R \
  --sample-dir  $S \
  --sample-id   SAMPLE_ID \
  --sex         male \
  --reference   auto \
  --gene-panel  lymphoid        # key findings need at least one gene list
```

This writes `SAMPLE_ID_report.html` and `SAMPLE_ID_small_variants.tsv.gz` to the current
directory. Matched and tumour-only runs take the same command: every input file is discovered
from the sample directory, and so is the run mode unless `--run-mode` states it.

## All options

```
--sample-dir      Path to the sample output directory (required)
--sample-id       Sample identifier (default: directory name)
--reference       t2t | hg38 | auto  (default: auto — read from the VEP somatic VCF's
                  ##contig lines, else the Severus VCF's, else t2t)
--sex             male | female | XY | XX  (required). Sets the expected chrX/chrY copy number
--gene-panel      none | builtin panel name (lymphoid, sarcoma) | path to a custom TSV
                  (default: none — nothing ticked on load). Repeatable: pass it several
                  times to open with several lists ticked (union). See Gene panels.
--gene-lists-dir  Directory to take the builtin panels from, replacing the bundled
                  assets/gene_lists (default: the bundled set). See Gene panels.
--embed           slim | full  (default: slim). Which small-variant rows and which
                  secondary figures are embedded in the HTML; see Report size.
--run-mode        paired | tumour_only  (default: inferred — matched when any normal-side QC
                  file is found). Sets the header badge and the tumour-only caveat; the
                  pipeline passes it, since a pair whose QC was skipped has no normal-side file
--output          Output HTML path  (default: <sample-id>_report.html in current dir).
                  The small-variant TSV is written to the same directory.
--title           Report title
--version, -V     Print the version and exit
```

## Key findings

A **finding** is an event in a gene of a ticked gene list that carries evidence of being a
driver. There is no built-in cancer-gene catalogue, by decision: a gene is only considered when
it is on a list. Findings are computed for every embedded list at render time, and the table
shows the rows of the lists ticked in the report. A gene's `role` column in the list
(oncogene / tumour suppressor, see [Gene roles](#gene-role-role)) adjusts some tiers.

The table has **one row per gene**, listing each of its alterations (for example a homozygous
deletion from copy number and the structural variant that causes it); the tier shown is the
gene's strongest, and genes are ordered by it. A rearrangement with breakpoints in two listed
genes is one row for the pair (`EWSR1–ATF1`; an en dash, because which partner is 5′ is not
established), matched by either gene's list or a Custom… entry for either gene.

**Tumour-only samples.** Without a matched normal, small-variant findings can include germline
variants. The report says so above the table and shows VAF for judgement (a heterozygous
germline variant sits near 0.5, or near 1 under LOH). When the pipeline reports that ClairS-TO's
germline tagging (Verdict) was skipped — it skips above purity 0.6 — a second warning names the
reason and the purity used, here and in the small-variant section: the calls then likely contain
many germline variants.

### Small variants

Considered when the gene is in the list's small-variant scope:

| Tier | Rule |
|---|---|
| 1 | ClinVar pathogenic or likely pathogenic ("Conflicting classifications" never counts) |
| 1 | VEP impact HIGH (truncating: stop gained, frameshift, splice donor/acceptor, start lost, transcript ablation; or another HIGH consequence) |
| 2 | VEP impact MODERATE with a damaging prediction: AlphaMissense *likely pathogenic*; where AlphaMissense gives no call, EVE pathogenic; failing both, SIFT *deleterious* together with PolyPhen *probably damaging* |

Role adjustments: a truncating variant in an **oncogene** drops to tier 2 (unless ClinVar
calls it pathogenic); a qualifying variant in a **tumour suppressor** whose gene also shows
loss of heterozygosity is tier 1 and reads "possible biallelic (LOH)".

### Structural variants

Both callers are evaluated. Calls of the two callers whose breakpoints both agree within
**1 kb** (and within half the event's size, so neighbouring small deletions stay apart) are one
event: it is listed once, naming both callers, and typed by whichever caller typed it — a
Severus "intra-chr breakend" matching a SAVANA deletion reads as a deletion. Within a gene,
events are listed per event class.

| Tier | Rule |
|---|---|
| 1 | A breakpoint lies inside the gene |
| 1 | A deletion or duplication of **at most 3 Mb** covers the gene |
| 1 | Padfoot rates the event `HIGH*` and names the gene |
| 2 | One of the first two, but Padfoot annotated the call and rated it below HIGH |
| 2 | A deletion/duplication/insertion with a breakpoint in the gene and neither a VEP nor a Padfoot annotation (exon effect unknown) |
| 3 | A deletion/duplication/insertion with a breakpoint in the gene that VEP annotated without a HIGH (exon) effect and Padfoot did not rate HIGH — mostly small intronic deletions. Translocations and other junctions are exempt: an intronic breakpoint is how a fusion looks |
| 3 | A non-translocation event of at most 3 Mb at an immunoglobulin / T-cell receptor locus (`IGH`, `IGK`, `IGL`, `TRA/D`, `TRB`, `TRG`, …): physiological V(D)J recombination |
| 3 | "Near": a breakend within 1 Mb of the gene (100 kb for other event types) without touching it; at most 3 such rows per gene, and none for an event already listed with stronger evidence for the gene |

A gene that merely lies inside a larger event (for example a 40 Mb inversion or deletion) is
not an SV finding; it is left to the copy-number rows.

### Copy number

From the gene-level calls of the copy-number section:

| Tier | Rule |
|---|---|
| 1 | Amplification — not listed at all for a tumour suppressor |
| 1 | Homozygous deletion — tier 2 for an oncogene |
| 3 | Either call in a **pseudoautosomal region** (reads from the chrX and chrY copies map to one of them, so a PAR gene can read CN 0 with nothing deleted) or at an **IG/TR locus** (physiologically rearranged and deleted in lymphocytes) |

Loss of heterozygosity is not a finding on its own; it is noted on the gene ("LOH in gene").
Tier 3 alterations are hidden until "Show lower-confidence findings" is clicked.

**Custom… lists.** A list typed into the report matches by gene symbol among the findings of the
embedded lists and among up to 500 small variants outside them that are ClinVar (likely)
pathogenic or HIGH impact (the strongest first). SV and copy-number findings exist only for
genes on an embedded list.

## Copy-number calls

Segments come from Wakhan's rank-1 `integer_profile.bed` (haplotype-resolved), else ASCAT's
`segments_raw.txt`. The purity and ploidy used are those of the same source, so the thresholds
and the segments share one ploidy scale; every source found (Wakhan cell purity, ASCAT, SAVANA)
is listed, because they can disagree — the purity biomarker is flagged when they differ by
0.15 or more.

Calls follow the COSMIC Cancer Gene Census CNA definition, on the total copy number of the
segments overlapping each gene:

| Ploidy | Amplification | Homozygous deletion |
|---|---|---|
| ≤ 2.7 | total CN ≥ 5 | total CN = 0 |
| > 2.7 | total CN ≥ 9 | total CN < ploidy − 2.7 |

Otherwise a gene is a **gain** or **loss** against the expected copy number, round(ploidy) —
and on chrX/chrY in a male, max(1, round(ploidy / 2)). A gene spanning several segments is
called amplified on its highest segment and homozygously deleted on its lowest, and marked
*partial* when its segments disagree. **LOH** means no minor allele on a non-zero segment
where two alleles are expected. A gene that no segment overlaps has no call. If ASCAT
inferred a different sex from the one passed with `--sex`, the section warns that the chrX/chrY
calls may be wrong.

Whole-genome duplication is called when the major allele is ≥ 2 over more than half of the
autosomal length; the LOH fraction is the autosomal fraction with no minor allele. ASCAT's own
WGD / LOH / GI metrics are shown alongside when present.

Gene positions come from the gene list's coordinates. A gene without coordinates is called
only where Wakhan's own per-gene file lists it (rows with a haplotype state above 100 are
dropped as implausible).

## Mutational signatures

COSMIC signatures (v3.6 when the pipeline used it) as fitted by SigProfilerAssignment. Per
context:

- **Few mutations** is flagged below 200 SBS, 50 DBS or 100 ID mutations; **poor fit** when the
  cosine similarity of the reconstruction is below 0.90.
- Signatures contributing **at least 5%** in an *unflagged* context are carried to the
  overview's "Dominant signatures". A flagged context is still shown in its section, but never
  reaches the overview.
- ID83 signatures are fitted against COSMIC's GRCh37 set whatever the alignment reference, and
  long-read indel calls in homopolymers can inflate replication-slippage and mismatch-repair
  indel signatures. On the long-read samples checked, the ID fit fell just below the cosine
  threshold.
- Aetiologies come from a bundled, deliberately conservative lookup
  (`assets/signatures/cosmic_v3.6_aetiology.tsv`): where COSMIC's proposed aetiology is not
  certain it reads *Unknown*, and a signature absent from the lookup reads *not in lookup*.
  Possible sequencing artefacts are marked.

## Report size and `--embed`

With the default `--embed slim`, reports measured **4.5–6.4 MB** on five real samples (v1.x
reports that embedded Wakhan's interactive plots ran to 30–90 MB); `--embed full` was 15 MB on
one of them. Slim mode embeds:

- small-variant rows that are VEP impact HIGH or MODERATE, ClinVar pathogenic / likely
  pathogenic, or in an embedded gene list's small-variant genes — the table says how many of
  the total are shown, and the full table is always in `<sample-id>_small_variants.tsv.gz`;
- the rank-1 Wakhan profile only;
- no genome-wide ReConPlot figure (at page width it is unreadable), and per-chromosome
  figures only for chromosomes with a high-impact call: a ShatterSeek chromothripsis call or
  an SV with a Padfoot cancer-gene rating, `HIGH_ONCO` or `HIGH_TSG` (always shown), or a
  tier 1–2 SV or copy-number key finding of a gene list (shown while that list is ticked).
  Padfoot's plain `HIGH` does not count: it is common in large fragile-site genes. Tickboxes above the
  figures switch each kind of reason off, and a chromothripsis reason follows the ShatterSeek
  tier selector;
- every ReConPlot **focus** figure (`reconplot/<combination>/focus/`, the regions and genes the
  pipeline was asked to draw), as a "Your regions" block that opens the ReConPlot section; the
  automatic per-chromosome figures then follow in a collapsed block;
- the first six CoRAL amplicon plots;
- no ASCAT diagnostic plots (GC correction, tumour separation).

`--embed full` embeds all of the above. Every count in the report (TMB, small-variant counts,
key findings) is computed from the full table in either mode.

**Wakhan's interactive Plotly plots are no longer embedded.** They remain in the pipeline's
output, under `wakhan/solution_rank_<N>/integer_profile.html` and
`wakhan/*_heatmap_ploidy_purity.html`; the report names them.

> **Changed in 2.0.0** (tested in the pipeline first through the pre-releases `2.0.0rc1` and
> `2.0.0rc2`; the release itself changes only `VERSION`, `recipe/meta.yaml` and this note):
> - The genome circos **labels the key findings**: the tier 1–2 genes of the ticked gene
>   list(s), outside the ideogram, bold for tier 1. They follow the chips, like the overview;
>   with no list ticked there are none.
> - **ReConPlot**: slim mode drops the genome-wide figure and embeds per-chromosome figures for
>   chromosomes with a high-impact call (see [Report size](#report-size-and---embed)).
> - A **light / dark mode toggle** in the header. The choice is remembered by the browser; a
>   printout is always light.
> - New **clinical overview** at the top: key findings by gene list, genome-wide biomarkers
>   and a print layout that puts the overview on the first page(s). It replaces the summary
>   cards. The old **"SNVs" card counted SNVs and indels together**; its successor is the
>   "Small variants" biomarker, which splits the two.
> - **Key findings need a gene list.** With none ticked the overview lists no findings. The
>   default is still `--gene-panel none`, so pass a list (or tick one in the report).
> - **SAVANA** is read as a second SV caller, next to Severus. SV counts are given **per
>   caller and never summed** (two callers report the same rearrangements); the genome-wide
>   circos draws one caller — Severus when it has calls, else SAVANA. With both, one
>   rearrangement can appear twice in the SV table; filter on `caller`.
> - SAVANA's same-chromosome breakend pairs are typed DEL / DUP / INV from `BP_NOTATION` and
>   treated like Severus' own spans (size, 100 kb panel window, no breakend arc); the hidden
>   `svtype_raw` column keeps `BND`.
> - **Padfoot** SV annotations (`padfoot_*` columns) and **ReConPlot** figures are shown.
> - New sections: **Mutational signatures**, **Genome instability** (ShatterSeek, ecDNA) and
>   **gene-level copy number**.
> - **Wakhan's interactive plots are gone from the report**, replaced by a static
>   haplotype-resolved profile; the originals are in the pipeline's `wakhan/` output. Reports
>   shrink from tens of MB to about 5 MB.
> - The genome circos shows **SNV density per 1 Mb** instead of one point per SNV coloured by
>   substitution class; the substitution classes are in the signatures section.
> - `--embed slim|full` (default `slim`) and the `<sample-id>_small_variants.tsv.gz` sidecar.
>   A slim report's small-variant table holds only the relevant rows (see above).
> - **samtools flagstat "mapped" now counts every mapped record**, secondary and
>   supplementary included — what samtools means by it. It previously showed the *primary*
>   mapped count under that name; that figure is now its own row, "Primary mapped reads".
> - cramino 2.0 is read by exact key: "Read N50" is the full-read-length N50, and the aligned-
>   portion N50 is shown beside it.
> - The `lymphoid` list has 234 genes (was 235): `FAM46C` dropped as a duplicate of `TENT5C`,
>   `PRKBC`/`RCK` renamed `PRKCB`/`DDX6`, and three hg38 rows that carried CHM13 coordinates
>   fixed. Every row now has coordinates (for gene-level copy number) and the small-variant
>   rows carry a `role`. See [`assets/gene_lists/README.md`](assets/gene_lists/README.md).
> - Gene lists may carry a `role` column (oncogene / tumour suppressor), used by the key-finding
>   tiers.

> **Changed in v1.1.0:**
> - `--mode` and `--somatic-vcf` were removed. Run mode is derived from whether normal-side
>   QC is present, and the VCF supplying VAF is now discovered (see below), so neither needs
>   to be declared. Scripts passing them will fail on an unknown option.
> - `--gene-panel` now defaults to `none` instead of `lymphoid`: reports open unfiltered
>   unless a panel is asked for. Pass `--gene-panel lymphoid` to restore the old default.

> **Changed in v1.2.1:**
> - The SV table is **one row per rearrangement**, not one per breakend record: Severus's
>   `_1`/`_2` mate records are collapsed, and each row carries both loci
>   (`chrom_a`/`pos_a`, `chrom_b`/`pos_b`) plus a `svclass` separating interchromosomal
>   translocations from intra-chromosomal breakends. The SV count and the circos links
>   halve accordingly — they were double-counting. The single `gene_hits` column is
>   replaced by per-breakend `gene_a`/`gene_b`, and a `panel_hit` column names which panel
>   gene matched, on which side, and how: each entry reads `GENE (side, how)` — for
>   example `RB1 (span, direct)` when the span overlaps the gene, or `RB1 (B, 188.5 kb)`
>   when breakend B only fell inside the search window. The windows are wide (1 Mb around a
>   BND), so that second token is what separates a disrupted gene from a nearby one.
> - SVs are matched against a gene panel by **coordinate**, not gene symbol, whenever the
>   panel carries `chrom`/`start`/`end` — see [Gene panels](#gene-panels). The bundled
>   `lymphoid.tsv` is replaced by `lymphoid.hg38.tsv` and `lymphoid.t2t.tsv`; a custom
>   symbol-only TSV still works and still matches on symbols.

> **Changed in v1.3.0:**
> - The categorical columns of both variant tables now filter by **tickbox dropdown**
>   instead of a free-text box: `consequence`, `impact` and `callers` on the small-variant
>   table, `svclass`, `svtype`, `impact`, `consequence` and `caller` on the SV table. Each
>   dropdown lists the values actually present in that sample with a row count, so it also
>   answers "what is even in this column?". Every other column keeps its text box, and the
>   table's own search box still does substring across all columns.
> - Ticking several values in one column is **OR**; ticking values in two columns is
>   **AND**. A ticked term matches a cell holding several — ticking `missense_variant` also
>   shows a row whose consequence is `missense_variant,splice_region_variant`. Counts are
>   over all rows and do not change as you filter, and `(none)` selects the rows with no
>   value in that column.
> - A column with fewer than two distinct values keeps its plain text box rather than
>   offering an empty dropdown. That is expected for `callers` when VEP produced its default
>   text output (which carries no per-variant caller).
> - Ticks survive changing the gene panel, are reflected in the "N shown" count, and are
>   honoured by the copy/CSV buttons. **Clear value filters** in the panel bar resets them —
>   it appears only while something is ticked, since a filter set on a header that has been
>   scrolled past is otherwise easy to lose track of.
> - Fixed: the "N small variants · N SVs shown" line went stale when a per-column filter was
>   used, having only followed the gene-panel selector.
> - `--gene-panel` is now **repeatable**, and the report's panel selector is a row of
>   checkboxes rather than a dropdown: tick any number and a row is kept if it hits **any**
>   of them. No box ticked is the unfiltered state, so the "All genes" entry is gone. With
>   two or more ticked, each `panel_hit` entry ends with the panel it matched in square
>   brackets — with one, the labels read exactly as before. `--gene-panel none` combined
>   with a real panel is now an error rather than a case where one quietly wins.
> - The SV table's footnote about how the panel matches now follows the ticked boxes; it
>   previously described the load-time panel and went stale the moment a reader switched.

> **Changed in v1.3.1:**
> - A calmer look. The page texture, gradients, card animations and shadows are gone;
>   every card, table and plot frame shares one flat surface. The circos drops the
>   per-chromosome Mb axes and the alternating track bands, and its legend sits beside the
>   plate. ASCAT plots are framed and captioned, and the QC tables match the variant tables.
>   Nothing about the data, filters or exports changed.
> - The tickbox dropdown counts are now **live**: each value's count is measured over the
>   rows passing every *other* filter (gene panel, text boxes, the other dropdowns), so it
>   says how many rows ticking it would leave. Values with no rows left grey out in place;
>   the unfiltered count moves into the tooltip.
> - Fixed: opening a dropdown on a right-hand column of a wide table no longer snaps the
>   table back to its first column.

> **Changed in v1.4.0:**
> - The small-variant table now shows the VEP plugin annotations LRSomatic adds:
>   AlphaMissense (`am_class`, `am_score`), ClinVar (`clinvar`, `clinvar_id`, linked to
>   the ClinVar record), CADD, REVEL and EVE where the run had them. SIFT and PolyPhen
>   are read from the `PolyPhen_SIFT` plugin on T2T, whose VEP cache carries neither — those
>   two columns were empty on every T2T report — and each predictor's class and score are
>   separate columns, so the class has a tickbox filter and the score sorts numerically.
>   Pathogenic-class calls are tinted.
> - A column appears only when the VEP run declared its source: a plugin that was not on
>   gives no column rather than an empty one, and `cosmic`/`dbsnp` are gone on T2T, whose
>   cache has no variation data. A new **Annotation sources** footnote under the table
>   names what was present and what was not.

## Gene panels

A gene list ("panel") does two jobs: it filters the two variant tables, and it is the set of
genes the [key findings](#key-findings) and the gene-level copy number are computed for.

`--gene-panel` only chooses which lists are ticked when the report opens; the rendered HTML
always contains every builtin list, so a reader can tick and untick lists (or paste a custom
gene list) in the browser without re-rendering. Tick any number and a variant or SV is kept if
it hits **any** of them; with none ticked the tables are unfiltered and the overview lists no
key findings.

Built-in panels live in `assets/gene_lists/`. Each is a TSV with a `gene` column (HGNC symbols)
and, optionally, `chrom`/`start`/`end` — which changes how structural variants are matched:

| Panel columns | Small variants | Structural variants |
|---|---|---|
| `gene` only | symbol match | symbol match on the annotated breakend genes — no positional window |
| `gene, chrom, start, end` | symbol match | within **1 Mb of either breakend** of a BND, or **100 kb of the span** of any other type |

Coordinate matching is the reliable mode: whether VEP annotates a breakend with a gene symbol
at all depends on the sample's VEP invocation (1.6%–90% of breakend rows across the samples
measured), so a symbol-only panel can hide exactly the translocations it exists to find. The
note under the SV table says which mode is in force. Coordinates are also what gene-level copy
number is called on; a symbol-only list gets copy-number calls only where Wakhan's per-gene
file has the gene.

Because coordinates are only valid for one genome, a coordinate panel must declare its
reference (a leading `# reference: hg38` line, or a `reference` column) and a mismatch with the
rendered reference is a hard error. Builtins ship one file per reference and are offered as a
single entry, resolved against the detected one.

| Panel | Description |
|---|---|
| `lymphoid` | 234 non-Hodgkin lymphoma genes, merged from an SV list (including the IG/TR rearrangement partners) and a coding-mutation list; each gene scoped to the table it belongs in (`applies_to`), with oncogene / tumour-suppressor roles on the coding-mutation genes. `lymphoid.hg38.tsv`, `lymphoid.t2t.tsv` |
| `sarcoma` | 140 soft-tissue and bone sarcoma genes — tumour suppressors, amplification targets and recurrent fusion partners. No roles. `sarcoma.hg38.tsv`, `sarcoma.t2t.tsv` |

```bash
--gene-panel lymphoid                # open with the builtin lymphoid panel ticked
--gene-panel /path/to/my_genes.tsv   # must have a 'gene' column or be a single-column file

# Repeatable — open with both ticked, a row kept if it hits either:
--gene-panel lymphoid --gene-panel /path/to/my_genes.tsv
```

With more than one panel ticked, each `panel_hit` entry ends with the panel it matched in
square brackets (`MYC (A, direct) [lymphoid]`); with one, there is no suffix. Two TSVs sharing
a basename both stay selectable — the second is registered as `<name>-custom`.

A `--gene-panel` value that is neither `none`, a builtin name, nor an existing file is an error —
a typo will not silently produce an unfiltered report. `none` combined with a real panel is an
error too, rather than a case where one of the two quietly wins.

### Optional columns: `applies_to` and `role`

- `applies_to` scopes a row to one table: blank / `both` / `all` (both tables), `snv` / `small`
  (small variants only), `sv` / `structural` (SVs only). Coordinates are required only on rows
  that can match an SV.
- <a id="gene-role-role"></a>`role` says whether a gene is an oncogene (`OG`, `oncogene`), a
  tumour suppressor (`TSG`, `tumour suppressor`, `tumor suppressor`) or both (`both`,
  `OG/TSG`, `TSG/OG`); blank means *not declared*. It only adjusts key-finding tiers (see
  [Key findings](#key-findings)); a gene that two lists give different roles is treated as
  "both", which adjusts nothing.

Any other value in either column is a hard error. See
[`assets/gene_lists/README.md`](assets/gene_lists/README.md) for the full file format and how the
bundled lists were built.

### Supplying your own set of builtins

`--gene-lists-dir` points the builtin lookup at another directory, **replacing** the bundled
set rather than adding to it:

```bash
--gene-lists-dir /path/to/panels --gene-panel mypanel
```

Files there follow the same convention — `<name>.tsv`, or a `<name>.hg38.tsv` /
`<name>.t2t.tsv` pair resolved against the rendered reference — and every one of them is
embedded in the report as a tickbox, exactly as the bundled panels are. This is what lets a
pipeline curate its own panel set without a release of this tool: a panel passed this way is a
first-class builtin, whereas the same file passed as `--gene-panel /path/to/mypanel.tsv` is a
user panel, named after the file and only embedded because it was named.

## Expected input layout

The `--sample-dir` must be the root of a single-sample LRSomatic output (`<outdir>/<id>/`).
Files are discovered **recursively**, mostly by their distinctive filename suffix, so the layout
below is what the pipeline publishes but not a requirement — except where a directory is named
as part of the pattern. When several files match, one named for the sample is preferred.

| Section | Files read (relative to the sample directory) |
|---|---|
| Small variants | `*_SOMATIC_VEP.vcf.gz` (variant set and annotation); VAF / depth / genotype / phase set from `variants/phased/somatic_smallvariants.vcf.gz`, else `variants/clairs{,to}/somatic.vcf.gz`, else any non-germline VCF in those directories; ClairS-TO germline tagging status from `variants/clairsto/<id>.verdict_status.tsv` (`status`, `purity`, `ploidy`, `purity_source`, `message`; status `applied`, `skipped` or `not_run`) |
| Structural variants | Severus: `severus_somatic.vcf.gz`, `*_SV_VEP.vcf.gz` (fallback `SV_filtered_with_gene_annotations.tsv`). SAVANA: `*.pon_filtered.classified.somatic.vcf[.gz]`, else `*.classified.somatic.vcf[.gz]`; `*_VEP_SAVANA.vcf.gz`. Padfoot: `padfoot/<combination>/annotated_svs.tsv` (`severus_wakhan` for Severus, `savana` for SAVANA) |
| ReConPlot | `reconplot/<combination>/genome_wide/*_genome_wide.png`, `reconplot/<combination>/per_chromosome/*_chr*.png`, `reconplot/<combination>/focus/*.png` |
| Copy number | ASCAT: `*.segments_raw.txt`, `*.purityploidy.txt`, `*.metrics.txt` (beside them), `*.tumour.{ASCATprofile,ASPCF,sunrise}.png` and diagnostic PNGs. Wakhan: `wakhan/solutions_ranks.tsv`, `wakhan/solution_rank_<N>/integer_profile.bed` and `…/genes/genes_copynumber_states.bed` (also `solution_<N>/` and the `repository_name` directory). SAVANA: `*_fitted_purity_ploidy.tsv` |
| Genome instability | `shatterseek/<combination>/*.shatterseek.tsv`; `ecdna/coral/*_CNV_SEEDS.bed`, `ecdna/coral/{cycles,reconstruct}/*_amplicon_summary.txt`, `ecdna/coral/plots/*_amplicon<N>_{graph,cycles}.png`, `ecdna/amplicon_classifier/*_amplicon_classification_profiles.tsv`, `*_gene_list.tsv`, `*_result_table.tsv` |
| Mutational signatures | `signatures/**/COSMIC_v<ver>/{SBS96,DBS78,ID83}/Assignment_Solution/**/Assignment_Solution_Activities.txt` and `Assignment_Solution_Samples_Stats.txt`; matrices `*.SBS96.all`, `*.DBS78.all`, `*.ID83.all` |
| QC (tumour) | mosdepth `*.mosdepth.summary.txt` (or `*.summary.txt` under `qc/` or `mosdepth/`) and `*.global.dist.txt`; cramino `*_cramino.txt` (the aligned-BAM one); samtools `*.flagstat`, `*.stats`; `qc/whatshap_stats/*_whatshap_stats.tsv` |
| QC (normal) | the same QC files under any `normal/` directory (e.g. `qc/normal/`) — their presence is what makes the run *matched*, unless `--run-mode` says otherwise |

Missing inputs are handled gracefully: the corresponding section or biomarker says it was not
run or not found.

**Small variants come from the VEP annotation only.** `*_SOMATIC_VEP.vcf.gz` defines the
variant set; VAF, depth, genotype and phase set are joined from the VCF that VEP annotated.
Runs predating `variants/phased/` yield VAF and depth but no phase set. If none is found the
table still renders, without those columns.

A footnote under the variant table names the file those columns actually came from and how
many variants they cover. **One VCF supplies them for the whole table**, so if the run
combined several somatic callers by consensus, the VAF, depth, genotype and phase set come
from whichever caller won each merge and need not match the `callers` column beside them.
Variants reported by more than one caller are highlighted in that column, and the footnote
says so. (Only the joined columns are ambiguous — the variant *set* is taken per record from
the VEP file.)

VEP writes indels at a different position and sometimes in a different allele notation than
the VCF it was given, so the join is made on a normalised key — see `variant_key()` in
`R/parse_smallvariants.R`.

`*_SOMATIC_VEP.vcf.gz` ships in two formats. Usually it is VEP *default text output* despite
the `.vcf.gz` name. If VEP was run with `--vcf` it is a genuine VCF with a `CSQ` field, and
may be a *merged* multi-caller VCF carrying germline calls (DeepVariant, Clair3) alongside
somatic ones, tagged in `INFO/CALLER`; the report then keeps only `PASS` records from a
somatic caller (ClairS, ClairS-TO, DeepSomatic). Both formats are detected automatically.

### Not covered: methylation

There is no methylation section. The only methylation output the pipeline publishes is
`methylation/<type>/modkit_pileup/<id>.bed.gz` — measured at 38–42 GB gzipped per sample,
unfiltered and with no tabix index, which cannot be read at render time.

`modkit pileup` already runs with `--bgzf`, so emitting a `tabix -p bed` index next to the
`.bed.gz` would be enough to unblock this: region queries on an indexed pileup measured
~0.16 s per Mb, making a binned genome-wide profile or per-locus lookup practical.

## Supported references

| `--reference` | Cytobands source | chr1 length |
|---|---|---|
| `t2t` | CHM13v2.0 | 248,387,328 bp |
| `hg38` | GRCh38 (UCSC) | 248,956,422 bp |

Auto-detection reads `##contig` lines from the VEP somatic VCF.

## Installing

### Container (recommended)

Every release is published as a container, so nothing needs installing (`<version>` is the
content of `VERSION`, e.g. `2.0.0`):

```bash
# Docker / Podman
docker run --rm ghcr.io/ljwharbers/lrsomatic-report:<version> render_report.R --version

# Apptainer / Singularity — a native SIF, no docker:// conversion needed
apptainer pull lrsomatic-report.sif oras://ghcr.io/ljwharbers/lrsomatic-report-sif:<version>
apptainer exec lrsomatic-report.sif render_report.R --version
```

`render_report.R` is on `PATH` inside the image and the tool is installed at
`/opt/conda/share/lrsomatic_report` — the same layout `recipe/build.sh` produces, so a
future Bioconda package is a drop-in replacement.

Images are built for `linux/amd64` only, by `.github/workflows/container.yml` on each `v*`
tag. The build refuses to run unless the tag, `VERSION` and `render_report.R --version`
all agree, and unless `container/environment.yml` and `recipe/meta.yaml` name the same
packages.

### R package requirements (running from a checkout)

`recipe/meta.yaml` is the source of truth for runtime dependencies — it is what the
Bioconda package is built from, and `container/environment.yml` (the image's pinned
version of the same list) is checked against it at build time. The list below mirrors it;
if they ever disagree, the recipe is right. v2.0 added no dependency.

Install in your R environment if missing:

```r
install.packages(c("data.table", "dplyr", "DT", "htmltools", "optparse",
                   "quarto", "yaml", "ggplot2", "svglite", "knitr",
                   "R.utils", "base64enc"))
BiocManager::install("circlize")
```

Plus the `quarto` CLI itself. `R.utils` is not called directly — `data.table::fread()`
requires it to read the gzipped VCFs.

Tested with R 4.4.1 and Quarto 1.5.57.

## Repository structure

```
lrsomatic_report/
├── bin/render_report.R          CLI entrypoint
├── R/
│   ├── utils.R                  Shared helpers (gene panels, JS serialisation, plot embedding)
│   ├── references.R             Cytoband + chrom-length loading, reference auto-detection
│   ├── locate_outputs.R         Discover small-variant, ASCAT and QC files in a sample directory
│   ├── parse_smallvariants.R    VEP text + VCF parsers; variant table, slim embed, SNV density
│   ├── parse_severus.R          SV VCF parsing (Severus + SAVANA), SV table, panel matching
│   ├── parse_padfoot.R          Padfoot annotation parser and join
│   ├── parse_ascat.R            ASCAT, Wakhan and SAVANA copy-number / purity parsers
│   ├── parse_qc.R               mosdepth, cramino, flagstat, samtools stats parsers
│   ├── circos.R                 draw_circos() — the genome-wide circos SVG
│   ├── circos_bnd.R             Breakend circos: selects and serialises it for the browser
│   ├── sections.R               Section-module registry
│   ├── sections/                One module per section: sv, whatshap, cna, shatterseek,
│   │                            ecdna, reconplot, signatures
│   └── key_findings.R           Clinical overview: key findings and biomarkers
├── templates/                   per_sample.qmd + sections/_*.qmd (Quarto)
├── assets/
│   ├── references/{t2t,hg38}/  Cytobands + chrom lengths (bundled, no network needed)
│   ├── gene_lists/             lymphoid / sarcoma per reference + README (override with
│   │                           --gene-lists-dir)
│   ├── signatures/             COSMIC aetiology lookup
│   ├── styles/                 report.scss and its section / print partials
│   └── js/                     bnd_circos.js, circos_drivers.js, facet_filter.js,
│                               key_findings.js — inlined at render time
├── tests/                      testthat unit tests + tests/js (node, no dependencies)
├── recipe/                     Conda recipe (Bioconda); build.sh defines the install layout
├── container/                  environment.yml (pinned deps, checked against the recipe),
│                               smoke test and release checks
├── Dockerfile                  Container image, reproducing recipe/build.sh's layout
└── VERSION                     Single source of truth; --version reads it
```

## Roadmap

- Cohort report (oncoprint, recurrence tables across multiple samples)
- Methylation, once the pipeline publishes an indexed pileup (see above)
