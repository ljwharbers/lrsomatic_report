# Gene Panel Selector Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the static per-render gene panel with an interactive dropdown in the HTML report that switches between all built-in panels, "All genes", and a custom text input — with the default panel pre-selected on page load.

**Architecture:** All variant data is embedded in the HTML at render time (no panel filtering at the R layer); JavaScript `$.fn.dataTable.ext.search` handles row visibility client-side. Both the SNV and SV tables are driven by a single shared dropdown. The default panel (first built-in, or the `--gene-panel` CLI arg) is applied immediately on DataTables initialisation.

**Tech Stack:** R (data.table), Quarto / HTML, DataTables (DT R package), vanilla JS, testthat

---

## File map

| File | Change |
|---|---|
| `R/utils.R` | Add `load_all_gene_panels(assets_dir)` |
| `R/parse_smallvariants.R` | Make `gene_panel` arg in `build_variant_table()` optional (NULL → no filter) |
| `R/parse_severus.R` | Make `gene_panel` arg in `build_sv_table()` optional (NULL → no filter, no row explosion) |
| `bin/render_report.R` | Call `load_all_gene_panels()`, pass `all_panels` + `default_panel` to template params |
| `templates/per_sample.qmd` | Add `all_panels`/`default_panel` params; unfiltered table calls; JS panel selector |
| `tests/testthat/setup.R` | Source R modules for all tests |
| `tests/testthat/test-utils.R` | Tests for `load_all_gene_panels()` |
| `tests/testthat/test-smallvariants.R` | Tests for updated `build_variant_table()` |
| `tests/testthat/test-severus.R` | Tests for updated `build_sv_table()` |

---

## Task 1: Bootstrap testthat

**Files:**
- Create: `tests/testthat/setup.R`
- Create: `tests/testthat/test-utils.R` (placeholder)

- [ ] **Step 1: Create `tests/testthat/setup.R`**

```r
# Assumes tests run from repo root (default for testthat::test_dir())
library(testthat)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})
source("R/utils.R")
source("R/parse_smallvariants.R")
source("R/parse_severus.R")
```

- [ ] **Step 2: Create `tests/testthat/test-utils.R` with a trivial passing test**

```r
test_that("setup loads without error", {
  expect_true(is.function(load_gene_panel))
})
```

- [ ] **Step 3: Verify the suite runs**

```bash
Rscript -e "testthat::test_dir('tests/testthat')"
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1 ]`

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/setup.R tests/testthat/test-utils.R
git commit -m "test: bootstrap testthat scaffold"
```

---

## Task 2: Add `load_all_gene_panels()` to `utils.R`

**Files:**
- Modify: `R/utils.R`
- Modify: `tests/testthat/test-utils.R`

- [ ] **Step 1: Write the failing test**

Add to `tests/testthat/test-utils.R`:

```r
test_that("load_all_gene_panels returns named list of character vectors", {
  # Use the real assets dir
  panels = load_all_gene_panels(file.path(getwd(), "assets"))
  expect_type(panels, "list")
  expect_true(length(panels) >= 1)
  expect_true("lymphoid" %in% names(panels))
  expect_type(panels[["lymphoid"]], "character")
  expect_true(length(panels[["lymphoid"]]) > 0)
})

test_that("load_all_gene_panels panel names are lowercase filenames without extension", {
  panels = load_all_gene_panels(file.path(getwd(), "assets"))
  expect_true(all(names(panels) == tolower(names(panels))))
})
```

- [ ] **Step 2: Run to confirm failure**

```bash
Rscript -e "testthat::test_dir('tests/testthat')"
```

Expected: FAIL with "could not find function \"load_all_gene_panels\""

- [ ] **Step 3: Implement in `R/utils.R`**

Add after the existing `resolve_gene_panel` function:

```r
# Load all gene panels from assets/gene_lists/*.tsv
# Returns a named list (name = panel name, value = character vector of gene symbols)
load_all_gene_panels = function(assets_dir) {
  tsv_files = Sys.glob(file.path(assets_dir, "gene_lists", "*.tsv"))
  if (length(tsv_files) == 0) return(list())
  panels = lapply(tsv_files, load_gene_panel)
  names(panels) = tools::file_path_sans_ext(basename(tsv_files))
  panels
}
```

- [ ] **Step 4: Run to confirm pass**

```bash
Rscript -e "testthat::test_dir('tests/testthat')"
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 3 ]`

- [ ] **Step 5: Commit**

```bash
git add R/utils.R tests/testthat/test-utils.R
git commit -m "feat: add load_all_gene_panels() to utils"
```

---

## Task 3: Make `build_variant_table()` panel-optional

**Files:**
- Modify: `R/parse_smallvariants.R`
- Create: `tests/testthat/test-smallvariants.R`

The change: when `gene_panel` is `NULL`, skip the panel-filter block and return all variants.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-smallvariants.R`:

```r
# Minimal synthetic VEP data (mimics parse_vep_text output)
make_vep = function() {
  data.table(
    chrom        = c("chr1", "chr1", "chr2"),
    pos          = c(100L, 200L, 300L),
    ref          = c("A", "C", "G"),
    alt          = c("T", "G", "A"),
    symbol       = c("MYC", "TP53", "KRAS"),
    gene_id      = c("ENSG001", "ENSG002", "ENSG003"),
    consequence  = c("missense_variant", "stop_gained", "synonymous_variant"),
    impact       = c("MODERATE", "HIGH", "LOW"),
    hgvsp        = c("p.A1T", "p.Q2*", NA_character_),
    existing     = c(NA_character_, NA_character_, NA_character_),
    dbsnp        = c(NA_character_, NA_character_, NA_character_),
    cosmic       = c(NA_character_, NA_character_, NA_character_),
    sift         = c(NA_character_, NA_character_, NA_character_),
    polyphen     = c(NA_character_, NA_character_, NA_character_)
  )
}

test_that("build_variant_table with NULL panel returns all variants", {
  vep = make_vep()
  result = build_variant_table(vep, list(), gene_panel = NULL)
  expect_equal(nrow(result), 3L)
  expect_true(all(c("MYC", "TP53", "KRAS") %in% result$symbol))
})

test_that("build_variant_table with panel filters to panel genes", {
  vep = make_vep()
  result = build_variant_table(vep, list(), gene_panel = c("MYC", "TP53"))
  expect_equal(nrow(result), 2L)
  expect_false("KRAS" %in% result$symbol)
})

test_that("build_variant_table with empty panel returns empty table", {
  vep = make_vep()
  result = build_variant_table(vep, list(), gene_panel = character(0))
  expect_equal(nrow(result), 0L)
})
```

- [ ] **Step 2: Run to confirm failure**

```bash
Rscript -e "testthat::test_dir('tests/testthat')"
```

Expected: FAIL — the NULL case currently still applies the panel filter (panel is `character(0)` by default or the current code treats NULL as no-genes).

- [ ] **Step 3: Update `build_variant_table()` in `R/parse_smallvariants.R`**

Change the signature and the filter block. Find this section (around line 136–148):

```r
# Build the small-variant display table:
#   canonical rows from VEP text, per-caller VAFs joined on position.
build_variant_table = function(vep_data, caller_list, gene_panel) {
  if (is.null(vep_data) || nrow(vep_data) == 0) return(NULL)

  # Impact ranking for deduplication
  impact_rank = c(HIGH = 1L, MODERATE = 2L, LOW = 3L, MODIFIER = 4L)
  vep_data[, impact_rank := impact_rank[impact]]
  vep_data[is.na(impact_rank), impact_rank := 5L]

  # Filter to gene panel (by gene symbol or Ensembl ID fallback)
  if (!is.null(gene_panel) && length(gene_panel) > 0) {
    vep_data = vep_data[symbol %in% gene_panel | gene_id %in% gene_panel]
  }
  if (nrow(vep_data) == 0) return(data.table())
```

Replace with:

```r
# Build the small-variant display table:
#   canonical rows from VEP text, per-caller VAFs joined on position.
# gene_panel: character vector of HGNC symbols to keep, or NULL to return all variants.
build_variant_table = function(vep_data, caller_list, gene_panel = NULL) {
  if (is.null(vep_data) || nrow(vep_data) == 0) return(NULL)

  # Impact ranking for deduplication
  impact_rank = c(HIGH = 1L, MODERATE = 2L, LOW = 3L, MODIFIER = 4L)
  vep_data[, impact_rank := impact_rank[impact]]
  vep_data[is.na(impact_rank), impact_rank := 5L]

  # Filter to gene panel when provided; NULL returns all variants
  if (!is.null(gene_panel) && length(gene_panel) > 0) {
    vep_data = vep_data[symbol %in% gene_panel | gene_id %in% gene_panel]
  }
  if (nrow(vep_data) == 0) return(data.table())
```

- [ ] **Step 4: Run to confirm pass**

```bash
Rscript -e "testthat::test_dir('tests/testthat')"
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 6 ]`

- [ ] **Step 5: Commit**

```bash
git add R/parse_smallvariants.R tests/testthat/test-smallvariants.R
git commit -m "feat: make gene_panel optional in build_variant_table (NULL = all variants)"
```

---

## Task 4: Make `build_sv_table()` panel-optional

**Files:**
- Modify: `R/parse_severus.R`
- Create: `tests/testthat/test-severus.R`

When `gene_panel = NULL`, skip row explosion and return one row per SV with the raw `gene_hits` field intact.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-severus.R`:

```r
make_sv_tsv = function() {
  data.table(
    ID           = c("SV1", "SV2", "SV3"),
    SVTYPE       = c("DEL", "DUP", "INV"),
    DETAILED_TYPE = c("DEL", "DUP", "INV"),
    START_CHROM  = c("chr1", "chr2", "chr3"),
    START_POS    = c(1000L, 2000L, 3000L),
    END_CHROM    = c("chr1", "chr2", "chr3"),
    END_POS      = c(5000L, 6000L, 7000L),
    SV_LEN       = c(4000L, 4000L, 4000L),
    VAF          = c(0.4, 0.5, 0.3),
    NHL_GENE_HITS = c("MYC;BCL2", "TP53", NA_character_),
    gene_hits    = c("MYC;BCL2", "TP53", NA_character_)
  )
}

test_that("build_sv_table with NULL panel returns one row per SV", {
  sv = make_sv_tsv()
  result = build_sv_table(sv, gene_panel = NULL)
  expect_equal(nrow(result), 3L)
})

test_that("build_sv_table with panel filters and explodes rows", {
  sv = make_sv_tsv()
  result = build_sv_table(sv, gene_panel = c("MYC", "TP53"))
  # SV1 contributes 1 row (MYC matches; BCL2 does not), SV2 contributes 1 row (TP53)
  expect_equal(nrow(result), 2L)
  expect_true(all(result$gene %in% c("MYC", "TP53")))
})

test_that("build_sv_table with empty panel returns empty table", {
  sv = make_sv_tsv()
  result = build_sv_table(sv, gene_panel = character(0))
  expect_equal(nrow(result), 0L)
})
```

- [ ] **Step 2: Run to confirm failure**

```bash
Rscript -e "testthat::test_dir('tests/testthat')"
```

Expected: FAIL — NULL panel currently still runs the explosion + filter path.

- [ ] **Step 3: Update `build_sv_table()` in `R/parse_severus.R`**

Replace the function (lines 87–121) with:

```r
# Build the SV display table.
# gene_panel: character vector of HGNC symbols to keep, or NULL to return all SVs (one row each).
build_sv_table = function(sv_tsv, gene_panel = NULL) {
  if (is.null(sv_tsv) || nrow(sv_tsv) == 0) return(data.table())

  sv_tsv = copy(sv_tsv)

  # When no panel is supplied, return one row per SV without explosion
  if (is.null(gene_panel)) {
    display_cols = intersect(
      c("ID", "SVTYPE", "DETAILED_TYPE",
        "START_CHROM", "START_POS", "END_CHROM", "END_POS",
        "SV_LEN", "VAF", "NHL_GENE_HITS", "COSMIC_GENE_HITS",
        "NHL_NEAREST_GENE_HITS_1MBWINDOW"),
      names(sv_tsv)
    )
    return(sv_tsv[, ..display_cols])
  }

  # Panel-filtered path: explode multi-gene gene_hits, filter, return one row per gene×SV
  sv_tsv[, .ridx := .I]

  sv_long = sv_tsv[, {
    raw = as.character(gene_hits[1])
    genes = unique(trimws(unlist(strsplit(raw, "[;,]+"))))
    genes = genes[nchar(genes) > 0 & genes != "-" & toupper(genes) != "NA"]
    if (length(genes) == 0) genes = NA_character_
    list(gene = genes)
  }, by = .ridx]

  sv_long = merge(sv_long, sv_tsv, by = ".ridx")
  sv_long[, .ridx := NULL]
  sv_tsv[,  .ridx := NULL]

  if (length(gene_panel) > 0) {
    sv_long = sv_long[!is.na(gene) & gene %in% gene_panel]
  }
  if (nrow(sv_long) == 0) return(data.table())

  display_cols = intersect(
    c("gene", "ID", "SVTYPE", "DETAILED_TYPE",
      "START_CHROM", "START_POS", "END_CHROM", "END_POS",
      "SV_LEN", "VAF", "NHL_GENE_HITS", "COSMIC_GENE_HITS",
      "NHL_NEAREST_GENE_HITS_1MBWINDOW"),
    names(sv_long)
  )
  sv_long[, ..display_cols]
}
```

- [ ] **Step 4: Run to confirm pass**

```bash
Rscript -e "testthat::test_dir('tests/testthat')"
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 9 ]`

- [ ] **Step 5: Commit**

```bash
git add R/parse_severus.R tests/testthat/test-severus.R
git commit -m "feat: make gene_panel optional in build_sv_table (NULL = all SVs, no explosion)"
```

---

## Task 5: Update `bin/render_report.R` to load and pass all panels

**Files:**
- Modify: `bin/render_report.R`

- [ ] **Step 1: Load all panels and add them to the `quarto_render` call**

After the existing gene panel resolution block (around line 48–57), add:

```r
# Load all available built-in panels for the dropdown
all_panels = load_all_gene_panels(file.path(repo_dir, "assets"))
default_panel = if (file.exists(file.path(repo_dir, "assets", "gene_lists",
                                            paste0(gene_panel, ".tsv")))) {
  gene_panel
} else {
  if (length(all_panels) > 0) names(all_panels)[1] else "custom"
}
```

Then in the `quarto::quarto_render()` call, add to `execute_params`:

```r
execute_params = list(
  sample_id     = sample_id,
  sample_dir    = sample_dir,
  reference     = reference,
  sex           = sex,
  gene_panel    = gene_panel,
  default_panel = default_panel,
  all_panels    = all_panels,
  title         = title,
  repo_dir      = repo_dir,
  outputs       = outputs
)
```

- [ ] **Step 2: Verify the script still parses without error (dry run)**

```bash
Rscript bin/render_report.R --help
```

Expected: prints help text without error.

- [ ] **Step 3: Commit**

```bash
git add bin/render_report.R
git commit -m "feat: load all gene panels in render_report.R, pass to template"
```

---

## Task 6: Template data layer — embed all variants, update params

**Files:**
- Modify: `templates/per_sample.qmd`

This task changes the Quarto params and the R setup chunk. No JS yet.

- [ ] **Step 1: Add new params to the YAML front matter**

In the `params:` block, replace:

```yaml
params:
  sample_id:     "SAMPLE"
  sample_dir:    ""
  reference:     "t2t"
  sex:           "female"
  gene_panel:    "lymphoid"
  title:         "LRSomatic Report"
  repo_dir:      "."
  outputs:       NULL   # named list from locate_outputs()
```

With:

```yaml
params:
  sample_id:     "SAMPLE"
  sample_dir:    ""
  reference:     "t2t"
  sex:           "female"
  gene_panel:    "lymphoid"
  default_panel: "lymphoid"
  all_panels:    NULL   # named list from load_all_gene_panels()
  title:         "LRSomatic Report"
  repo_dir:      "."
  outputs:       NULL   # named list from locate_outputs()
```

- [ ] **Step 2: Update the setup chunk to embed all variants and compute default-panel counts**

In the `{r setup}` chunk, replace the three lines:

```r
variant_table = build_variant_table(vep_data, caller_list, panel_genes)
...
sv_table = build_sv_table(sev_tsv, panel_genes)
```

With:

```r
# Full unfiltered tables — panel filtering is done client-side via JS
variant_table = build_variant_table(vep_data, caller_list, gene_panel = NULL)
...
sv_table = build_sv_table(sev_tsv, gene_panel = NULL)
```

Also update the summary card counts so they still reflect the default panel:

```r
# Summary counts use the default panel for the header cards
n_snv        = if (!is.null(vep_data)) nrow(unique(vep_data[, .(chrom, pos, ref, alt)])) else NA_integer_
n_sv         = if (!is.null(sev_vcf)) nrow(sev_vcf$nontrans) + nrow(sev_vcf$translocations) else NA_integer_
n_panel_vars = if (!is.null(variant_table) && length(panel_genes) > 0)
                 sum(variant_table$symbol %in% panel_genes, na.rm = TRUE) else 0L
n_panel_svs  = if (!is.null(sv_table) && length(panel_genes) > 0 && "NHL_GENE_HITS" %in% names(sv_table))
                 sum(vapply(sv_table$NHL_GENE_HITS, function(h)
                   any(trimws(unlist(strsplit(as.character(h), "[;,]+"))) %in% panel_genes),
                   logical(1)), na.rm = TRUE) else 0L
```

- [ ] **Step 3: Verify template renders without error using a known sample (or at minimum that R sources load cleanly)**

```bash
Rscript -e "
  repo_dir <- '.'
  source('R/utils.R'); source('R/parse_smallvariants.R'); source('R/parse_severus.R')
  cat('Sources OK\n')
"
```

Expected: `Sources OK` with no errors.

- [ ] **Step 4: Commit**

```bash
git add templates/per_sample.qmd
git commit -m "feat: embed full unfiltered variant/SV tables in template; default-panel counts in header cards"
```

---

## Task 7: Template UI layer — dropdown + JS client-side filtering

**Files:**
- Modify: `templates/per_sample.qmd`

This task adds the dropdown control and the DataTables search extension.

- [ ] **Step 1: Add the JS panel serialisation block**

Insert a new R chunk immediately after the `{r setup}` chunk (before any rendered output). This chunk serialises panel gene lists to JavaScript:

````markdown
```{r panel-js-data, results='asis'}
# Serialise panel gene lists as a JS const for use by the dropdown filter
all_p = if (!is.null(params$all_panels) && length(params$all_panels) > 0)
          params$all_panels else list()
js_panels = paste0(
  "const GENE_PANELS = {",
  paste(vapply(names(all_p), function(nm) {
    genes_json = paste0('"', all_p[[nm]], '"', collapse = ", ")
    paste0('"', nm, '": new Set([', genes_json, '])')
  }, character(1)), collapse = ",\n"),
  "};\n",
  'const DEFAULT_PANEL = "', params$default_panel, '";\n'
)
cat("<script>\n", js_panels, "\n</script>\n", sep = "")
```
````

- [ ] **Step 2: Replace the `gene-filter-ui` chunk with the dropdown control**

Replace the entire existing `{r gene-filter-ui}` chunk with:

````markdown
```{r gene-filter-ui, results='asis'}
# Build dropdown options: all built-in panels + All genes + Custom
panel_opts = ""
all_p = if (!is.null(params$all_panels) && length(params$all_panels) > 0)
          params$all_panels else list()
for (nm in names(all_p)) {
  selected = if (nm == params$default_panel) " selected" else ""
  label = paste0(toupper(substr(nm, 1, 1)), substr(nm, 2, nchar(nm)))  # title-case
  panel_opts = paste0(panel_opts,
    '<option value="', nm, '"', selected, '>', label, '</option>\n')
}
panel_opts = paste0(panel_opts,
  '<option value="__all__">All genes</option>\n',
  '<option value="__custom__">Custom&hellip;</option>\n')

cat(paste0('
<div style="margin-bottom:12px;">
  <label for="panel-select" style="font-weight:600; margin-right:8px;">Gene panel:</label>
  <select id="panel-select" style="padding:4px 8px; border-radius:4px; border:1px solid #ced4da;">',
    panel_opts,
  '</select>
  <span id="panel-variant-count" style="margin-left:12px; font-size:0.9em; color:#6c757d;"></span>
</div>
<div id="custom-gene-panel" style="display:none; margin-bottom:12px;">
  <label style="font-weight:600; display:block; margin-bottom:4px;">
    Custom genes (one per line, or comma/space separated, case-insensitive):
  </label>
  <textarea id="custom-gene-input"
    rows="4" style="width:400px; font-family:monospace; padding:6px;
                    border:1px solid #ced4da; border-radius:4px;"></textarea>
</div>
'))
```
````

- [ ] **Step 3: Add `elementId` and `initComplete` to both DT calls**

In the `{r small-variant-table}` chunk, add `elementId = "snv-table"` to `DT::datatable()` and add `initComplete` to its options list:

```r
DT::datatable(
  dt_display,
  rownames    = FALSE,
  filter      = "top",
  elementId   = "snv-table",
  extensions  = c("Buttons", "Scroller"),
  options     = list(
    dom        = "Bfrtip",
    buttons    = c("copy", "csv", "excel"),
    scrollX    = TRUE,
    scrollY    = "400px",
    scroller   = TRUE,
    deferRender = TRUE,
    pageLength = 25,
    columnDefs = list(list(className = "dt-left", targets = "_all")),
    initComplete = JS("function() { window.snvTableElem = this.api().table().node(); }")
  )
)
```

In the `{r sv-table}` chunk, add `elementId = "sv-table"` and `initComplete`:

```r
DT::datatable(
  sv_table,
  rownames   = FALSE,
  filter     = "top",
  elementId  = "sv-table",
  extensions = c("Buttons", "Scroller"),
  options    = list(
    dom       = "Bfrtip",
    buttons   = c("copy", "csv", "excel"),
    scrollX   = TRUE,
    scrollY   = "350px",
    scroller  = TRUE,
    deferRender = TRUE,
    pageLength = 25,
    initComplete = JS("function() { window.svTableElem = this.api().table().node(); }")
  )
)
```

- [ ] **Step 4: Add the JS filtering block**

Insert a new raw HTML block at the very end of the `.qmd` (after the footer line), containing the panel filter logic:

````markdown
```{=html}
<script>
(function () {
  // Panel filter — registered once, applied on every DataTables draw
  // Only targets the SNV and SV tables (identified by stored DOM node references)
  let activePanel = DEFAULT_PANEL;
  let customGenes = new Set();

  $.fn.dataTable.ext.search.push(function (settings, data) {
    const node = settings.nTable;
    if (node !== window.snvTableElem && node !== window.svTableElem) return true;
    const gene = (data[0] || "").toUpperCase();
    if (activePanel === "__all__")    return true;
    if (activePanel === "__custom__") return customGenes.has(gene);
    const panel = GENE_PANELS[activePanel];
    return panel ? panel.has(gene) : true;
  });

  function parseCustomGenes(text) {
    return new Set(
      text.split(/[\s,\n]+/).map(g => g.trim().toUpperCase()).filter(g => g.length > 0)
    );
  }

  function updateCount() {
    const snvTbl = window.snvTableElem ? $(window.snvTableElem).DataTable() : null;
    if (!snvTbl) return;
    const n = snvTbl.rows({ search: "applied" }).count();
    const el = document.getElementById("panel-variant-count");
    if (el) el.textContent = n + " variant" + (n === 1 ? "" : "s") + " shown";
  }

  function redraw() {
    if (window.snvTableElem) $(window.snvTableElem).DataTable().draw();
    if (window.svTableElem)  $(window.svTableElem).DataTable().draw();
    updateCount();
  }

  $(document).ready(function () {
    $("#panel-select").on("change", function () {
      activePanel = this.value;
      const isCustom = activePanel === "__custom__";
      $("#custom-gene-panel").toggle(isCustom);
      redraw();
    });

    let debounceTimer;
    $("#custom-gene-input").on("input", function () {
      clearTimeout(debounceTimer);
      debounceTimer = setTimeout(function () {
        customGenes = parseCustomGenes(document.getElementById("custom-gene-input").value);
        redraw();
      }, 300);
    });

    // Trigger initial draw once both tables have initialised
    // (initComplete stores window.snvTableElem / window.svTableElem)
    const waitForTables = setInterval(function () {
      if (window.snvTableElem || window.svTableElem) {
        clearInterval(waitForTables);
        redraw();
      }
    }, 100);
  });
})();
</script>
```
````

- [ ] **Step 5: Update the SV info chunk to remove the old static panel label**

Replace the existing `{r sv-info}` chunk body with:

```r
if (is.null(sv_table) || nrow(sv_table) == 0) {
  tags$div(class = "alert alert-info",
    if (is.null(sev_tsv)) "Severus gene-annotated SV file not found."
    else "No structural variants detected.")
}
```

(The dropdown selector already shows the active panel; no need to repeat it.)

- [ ] **Step 6: Commit**

```bash
git add templates/per_sample.qmd
git commit -m "feat: interactive gene panel dropdown with JS client-side filtering"
```

---

## Task 8: Smoke test

This task has no automated test (it requires a browser). Run manually against a real sample.

- [ ] **Step 1: Render a report against the DLBCL3_pooled sample**

```bash
Rscript bin/render_report.R \
  --sample-dir /staging/leuven/stg_00096/home/averham/LR_SOMATIC_T2T/DLBCL3_pooled \
  --sample-id  DLBCL3_pooled \
  --sex        male \
  --reference  auto \
  --output     /tmp/DLBCL3_pooled_report.html
```

- [ ] **Step 2: Check HTML file size**

```bash
du -sh /tmp/DLBCL3_pooled_report.html
```

If the file exceeds 150 MB, the report may be slow to open. In that case, open a follow-up plan to cap "All genes" to the union of all built-in panels rather than all 135k VEP rows.

- [ ] **Step 3: Open in browser and verify**

Checklist:
- [ ] Page loads and circos plot, summary cards, and QC sections all render
- [ ] Default panel (lymphoid) is pre-selected and its variants are shown in the SNV table
- [ ] Switching to another built-in panel (if one exists) updates both tables
- [ ] Selecting "All genes" shows all variants in both tables
- [ ] Selecting "Custom…" shows the textarea; typing `MYC` filters to MYC rows only
- [ ] Custom input is case-insensitive (`myc` matches MYC)
- [ ] The mosdepth coverage table (a third DataTable on the page) is NOT affected by the dropdown
- [ ] The variant count badge next to the dropdown updates on every switch
- [ ] Download buttons (CSV / Excel) export only the currently visible rows

- [ ] **Step 4: Final commit (if any minor fixes were needed)**

```bash
git add -p
git commit -m "fix: smoke-test corrections to gene panel selector"
```
