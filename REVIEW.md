# Review guide for lrsomatic_report

Guidance for automated (Claude) and human PR reviews. Every check below comes from a
problem that reached this repository before; `CLAUDE.md` has the full design notes. The
automated reviewer only sees the diff and cannot run R, Quarto or node, so it must never
claim that tests or a render pass or fail.

## Tests and checks

- CI (`.github/workflows/container.yml`) only builds the image and smoke-renders one tiny
  fixture. **Neither test suite runs in CI**: `Rscript -e 'testthat::test_dir("tests/testthat")'`
  and `node tests/js/*.js` are run by hand, plus a render of a real sample. Flag new parsing,
  selection or tiering logic that has no test in `tests/testthat/`, and new pure JS without a
  node test in `tests/js/`.
- A change to the rendered HTML, the JS or the CSS cannot be checked by the suites at all. Ask
  for a real-sample render (and a browser look for anything interactive) in the summary.
- Discovery of a new tool output needs a test in the **pipeline-staged layout** too (the
  `test-*.R` files with a staged twin), see below.

## Silent wrong answers (the failure mode to look for)

Most past bugs produced a plausible report with wrong or empty values, not an error:

- Joins that can silently yield `NA`: small variants join VEP to the VCF through
  `variant_key()`, whose coordinate space comes from the parser's `coord_space` column. A
  change to either side, or a hard-coded space, can make every indel's VAF `NA`.
- Columns read by position instead of by name (the VCF sample column, VEP/Padfoot columns),
  and FORMAT tags assumed fixed (`AF` vs `VAF` vs `AD`).
- `regmatches()` dropping non-matching elements, `fread` typing a bare-contig `CHROM` as
  integer, `ifelse()` on an all-`NA` test returning logical — all broke length or type
  assumptions before.
- A filter or panel test that returns "nothing found" where it should say "not run" or
  "not applicable" (an `N/A` card, a "not run" biomarker).

## R conventions

- Assignment is `=`, not `<-`. `data.table` first, `dplyr` for display only.
- Files in `R/` are `source()`d by the Quarto template, not a package. A new section is a
  module in `R/sections/<id>.R` calling `register_section()` plus `templates/sections/_<id>.qmd`
  (contract in `R/sections.R` and `CLAUDE.md`); flag discovery added to `locate_outputs.R` or
  the setup chunk of `templates/per_sample.qmd` instead.
- Every parser returns `NULL` (or an empty table) when its input is absent, and the section
  shows `section_notice()`; flag code that errors on a missing file.

## Pipeline coupling

- This tool is vendored into the LRSomatic pipeline (`assets/lrsomatic_report/`). The pipeline
  does **not** hand over its published directory: it builds a `sample_dir/` of symlinked
  **files**, never directories, with some files flat at the root and QC flat in
  `qc/tumor/`/`qc/normal/`. `list.files(recursive = TRUE)` does not descend into a symlinked
  directory. Flag a discovery pattern that relies on a directory name the staging does not
  reproduce, or that only works on the published layout. A new tool output also needs the
  pipeline's staging channel updated; say so in the summary.
- Dependencies: `recipe/meta.yaml` is the source of truth and must stay identical to the
  pipeline's `modules/local/lrsomaticreport/environment.yml`. Flag any new package (including
  namespaced `pkg::fn` calls, which `library()` greps miss), and any use of one not listed
  there. `jsonlite` in particular is deliberately not a dependency (see `js_rows()`).
- Releases: `VERSION` must equal the recipe's version, and the recipe's `sha256` is only
  recorded in a follow-up commit after the tag exists — a release PR keeps the old digest.

## R ↔ JS mirrors and the browser layer

- Several rules exist twice, in R and in the inline JS, and must change together: the SV
  panel windows (`SV_PANEL_WINDOW_*` / `svPanelHits()`), the `panel_hit` label, the
  `"__all__"` no-filter sentinel (`bin/render_report.R`, `_gene_filter.qmd`, the search hook
  in `per_sample.qmd`), and the key-findings row test (`kfRowMatches()`), which other
  chip-following parts reuse rather than copy. Flag a change to one side only.
- The JS resolves table columns **by name** (`window.SV_COLS`, `idx0()`), so reordering columns
  is safe but **renaming** one (`chrom_a`, `pos_a`, …) silently breaks the panel filter and
  empties the breakend circos.
- Hand-written serialisation: `js_quote()` escapes only double quotes, so free text going into
  a `<script>` must be sanitised (backslashes, newlines, `</`); `js_num()` is element-wise on
  purpose — vectorised `format()` pads (`c(1, 2.5)` → `"1.0"`) or switches to scientific.
- Quarto lowercases attributes on markup it processes, so SVG `viewBox` written in R markup
  arrives as `viewbox` and is ignored; set it from JS.
- Panel listeners (`window.registerPanelListener`) are called more than once per change, so
  they must be idempotent.

## Report size and embedding

- The report is self-contained: no network access at render time, no external scripts or
  fonts beyond what is bundled.
- Slim mode (`--embed slim`, the default) is what keeps reports at a few MB. Flag new
  base64-embedded figures or table rows that ignore `embed`, and slim-mode omissions that do
  not tell the reader where the rest is.

## Clinical rules and docs

- The README documents every clinical rule (tiers, copy-number thresholds, signature flags,
  which figures slim mode embeds) in plain language: a constant changed in `R/` must be
  changed there too. Flag a rule change without the matching README edit.
- Gene-list coordinates are reference-specific; a panel that declares a different reference
  must stop with an error, never load silently. Check new logic against both `hg38` and `t2t`.

## Styling

- Colours are CSS custom properties on `:root`. Dark mode goes through the SCSS
  `dark-mode` mixin in `assets/styles/report.scss` where it exists; flag a raw
  `@media (prefers-color-scheme: dark)` block, which would ignore the header toggle.
- Plot plates (circos, PNG figures) stay light in dark mode; their ink is fixed, not a theme
  token.

## Do not comment on

- Rendered reports (`*.html`), `*_small_variants.tsv.gz` sidecars and binary fixtures under
  `tests/fixtures/`.
- Line wrapping, comment wording or formatting.
- `HANDOFF.md` and other untracked notes.
