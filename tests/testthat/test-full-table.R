# The full small-variant table embedded in the HTML (--embed-full-table): payload writer,
# size cap, controls markup and the display spec assets/js/full_table.js rebuilds rows from

full_fixture = function() {
  data.table(
    symbol = c("TP53", "KRAS", NA), chrom = "chr17", pos = c(7675088L, 25245350L, 100L),
    ref = c("C", "C", "AT"), alt = c("T", "A", "A"),
    consequence = c("missense_variant", "missense_variant,splice_region_variant", "intron_variant"),
    hgvsc = c("c.524G>A", "c.35G>T", NA), vaf = c(0.4567, 0.125, NA), dp = c(40L, NA, 12L),
    flag = c(TRUE, FALSE, NA)
  )
}

sidecar_of = function(vt) {
  f = tempfile(fileext = ".tsv.gz")
  write_full_variant_table(vt, f)
  f
}

test_that("the payload decodes to exactly the sidecar's bytes, and gunzips to the table", {
  vt = full_fixture()
  f = sidecar_of(vt)
  p = full_table_payload(f)
  expect_equal(p$status, "embedded")
  expect_equal(p$bytes, file.size(f))
  raw_in = readBin(f, "raw", file.size(f))
  expect_identical(base64enc::base64decode(p$b64), raw_in)

  back = tempfile(fileext = ".tsv.gz")
  writeBin(base64enc::base64decode(p$b64), back)
  expect_equal(readLines(gzfile(back)), readLines(gzfile(f)))
  expect_equal(fread(back, sep = "\t"), fread(f, sep = "\t"))
  expect_equal(fread(back, sep = "\t")$hgvsc[1:2], vt$hgvsc[1:2])
})

test_that("the committed fixture the node test reads round-trips too", {
  f = file.path(repo_root, "tests", "fixtures", "full_table", "fixture_small_variants.tsv.gz")
  p = full_table_payload(f)
  expect_equal(p$status, "embedded")
  expect_identical(base64enc::base64decode(p$b64), readBin(f, "raw", file.size(f)))
})

test_that("the size cap and the option switch embedding off", {
  f = sidecar_of(full_fixture())
  size = file.size(f)
  expect_equal(full_table_payload(f, max_bytes = size)$status, "embedded")
  big = full_table_payload(f, max_bytes = size - 1)
  expect_equal(big$status, "too_large")
  expect_null(big$b64)
  expect_equal(big$bytes, size)
  off = full_table_payload(f, enabled = FALSE)
  expect_equal(off$status, "off")
  expect_null(off$b64)
  expect_equal(full_table_payload(file.path(tempdir(), "nope.tsv.gz"))$status, "missing")
  expect_equal(full_table_payload(NULL)$status, "missing")
  expect_equal(FULL_TABLE_EMBED_MAX_BYTES, 5 * 1024^2)
})

test_that("controls carry the payload inertly, and only when it was embedded", {
  f = sidecar_of(full_fixture())
  p = full_table_payload(f)
  expect_null(full_table_controls(full_table_payload(f, enabled = FALSE), 3))
  expect_null(full_table_controls(full_table_payload(f, max_bytes = 1), 3))
  expect_null(full_table_controls(NULL, 3))

  html = as.character(full_table_controls(p, 12345))
  expect_match(html, '<script type="application/octet-stream" id="snv-full-table-data"', fixed = TRUE)
  expect_match(html, paste0('data-filename="', basename(f), '"'), fixed = TRUE)
  expect_match(html, paste0(">", p$b64, "</script>"), fixed = TRUE)
  expect_match(html, "Download full table (12,345 variants, TSV)", fixed = TRUE)
  expect_false(grepl("full-table__show", html, fixed = TRUE))

  html_all = as.character(full_table_controls(p, 12345, show_all = TRUE))
  expect_match(html_all, '<button type="button" class="full-table__btn full-table__show" hidden>',
               fixed = TRUE)
})

test_that("the notice names the sidecar in every case, and says why it was not embedded", {
  f = sidecar_of(full_fixture())
  txt = function(p) gsub("\\s+", " ", as.character(full_table_where(p, "S_small_variants.tsv.gz")))
  expect_match(txt(full_table_payload(f)), "embedded in this report", fixed = TRUE)
  expect_match(txt(full_table_payload(f, enabled = FALSE)),
               "The complete table is in <code>S_small_variants.tsv.gz</code> next to this report.",
               fixed = TRUE)
  big = txt(full_table_payload(f, max_bytes = 1))
  expect_match(big, "<code>S_small_variants.tsv.gz</code>", fixed = TRUE)
  expect_match(big, "above the", fixed = TRUE)
  over = list(status = "too_large", bytes = 6.3 * 1024^2, max_bytes = FULL_TABLE_EMBED_MAX_BYTES)
  expect_match(txt(over), "at 6.3 MB compressed it is above the 5 MB limit", fixed = TRUE)
  expect_match(txt(modifyList(over, list(max_bytes = 2.5 * 1024^2))), "above the 2.5 MB limit",
               fixed = TRUE)
})

test_that("the display table and its column spec agree with what the DT shows", {
  vt = full_fixture()
  d = snv_display_table(vt)
  expect_true("vaf" %in% names(vt))                       # not modified by reference
  expect_equal(names(d), sub("^vaf$", "VAF%", names(vt)))
  expect_equal(d[["VAF%"]], c(45.7, 12.5, NA))
  expect_null(snv_display_table(NULL))

  spec = snv_column_spec(d)
  expect_match(spec, '["VAF%","vaf","number",100,1]', fixed = TRUE)
  expect_match(spec, '["pos","pos","number",1,null]', fixed = TRUE)
  expect_match(spec, '["flag","flag","logical",1,null]', fixed = TRUE)
  expect_match(spec, '["hgvsc","hgvsc","string",1,null]', fixed = TRUE)
  # one entry per display column, in display order
  expect_equal(lengths(regmatches(spec, gregexpr("\\],\\[", spec))) + 1L, ncol(d))
  expect_true(startsWith(spec, '[["symbol"'))
})

test_that("the multi-caller tint covers the full table's combinations when Show all is offered", {
  full = data.table(callers = c("clairs", "clairs,deepsomatic", "clairs,deepsomatic", NA,
                                "deepsomatic,severus", ""), vaf = 0.2)
  slim = snv_display_table(full[1:3])
  expect_equal(snv_multi_caller_values(slim), "clairs,deepsomatic")
  expect_equal(snv_multi_caller_values(slim, full),
               c("clairs,deepsomatic", "deepsomatic,severus"))
  # same column under the display transform, so the values match the swapped-in cells
  expect_equal(snv_multi_caller_values(slim, snv_display_table(full)),
               snv_multi_caller_values(slim, full))
  expect_equal(snv_multi_caller_values(slim[, !"callers"], full), character(0))
  expect_equal(snv_multi_caller_values(NULL, full), character(0))
})

test_that("parse_positive_mb reads --embed-full-table-max-mb", {
  expect_equal(parse_positive_mb("5"), 5)
  expect_equal(parse_positive_mb(" 2.5 "), 2.5)
  expect_true(is.na(parse_positive_mb("0")))
  expect_true(is.na(parse_positive_mb("-1")))
  expect_true(is.na(parse_positive_mb("five")))
  expect_true(is.na(parse_positive_mb("Inf")))
  expect_true(is.na(parse_positive_mb(NULL)))
})

test_that("parse_yes_no reads the --embed-full-table values", {
  expect_true(parse_yes_no("yes"))
  expect_true(parse_yes_no(" TRUE "))
  expect_false(parse_yes_no("no"))
  expect_false(parse_yes_no("0"))
  expect_true(is.na(parse_yes_no("maybe")))
  expect_true(is.na(parse_yes_no(NULL)))
  expect_true(is.na(parse_yes_no(NA)))
})
