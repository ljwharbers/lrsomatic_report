# The CLI contract the pipeline depends on: `render_report.R --version` prints the bare
# version and nothing else, so Nextflow can use it directly in an `eval` version topic.
# It must answer without the R package stack, which is what makes it cheap per task.

render_report = file.path(repo_root, "bin", "render_report.R")

test_that("--version prints the bare VERSION and exits 0", {
  out = suppressWarnings(system2("Rscript", c(render_report, "--version"),
                                 stdout = TRUE, stderr = TRUE))
  expect_equal(attr(out, "status"), NULL)
  expect_equal(out, readLines(file.path(repo_root, "VERSION"), warn = FALSE)[1])
})

test_that("VERSION agrees with the conda recipe", {
  # The container workflow refuses to build when the git tag disagrees with VERSION;
  # this catches the other half, a recipe left behind by a release.
  version = readLines(file.path(repo_root, "VERSION"), warn = FALSE)[1]
  recipe  = readLines(file.path(repo_root, "recipe", "meta.yaml"), warn = FALSE)
  set_ver = grep('^\\{%\\s*set version', recipe, value = TRUE)
  expect_length(set_ver, 1)
  expect_equal(sub('.*"([^"]+)".*', "\\1", set_ver), version)
})

test_that("--help lists the flags the pipeline passes", {
  out = suppressWarnings(system2("Rscript", c(render_report, "--help"),
                                 stdout = TRUE, stderr = TRUE))
  help = paste(out, collapse = "\n")
  for (flag in c("--sample-dir", "--sample-id", "--sex", "--reference",
                 "--gene-panel", "--gene-lists-dir", "--output", "--embed", "--run-mode")) {
    expect_true(grepl(flag, help, fixed = TRUE), info = flag)
  }
  # The help names both accepted --embed values
  expect_true(grepl("slim | full", help, fixed = TRUE))
})

test_that("--embed rejects anything but slim or full, before any rendering", {
  d = tempfile("sample")
  dir.create(d)
  out = suppressWarnings(system2("Rscript",
    c(render_report, "--sample-dir", d, "--sex", "male", "--embed", "bogus"),
    stdout = TRUE, stderr = TRUE))
  expect_equal(attr(out, "status"), 1L)
  expect_true(any(grepl("--embed must be 'slim' or 'full'", out, fixed = TRUE)))
  expect_false(dir.exists(file.path(getwd(), "._render")))
})

test_that("--embed-full-table is documented and rejects anything but yes or no", {
  out = suppressWarnings(system2("Rscript", c(render_report, "--help"),
                                 stdout = TRUE, stderr = TRUE))
  help = paste(out, collapse = "\n")
  expect_true(grepl("--embed-full-table", help, fixed = TRUE))
  expect_true(grepl("yes | no", help, fixed = TRUE))

  d = tempfile("sample")
  dir.create(d)
  out = suppressWarnings(system2("Rscript",
    c(render_report, "--sample-dir", d, "--sex", "male", "--embed-full-table", "maybe"),
    stdout = TRUE, stderr = TRUE))
  expect_equal(attr(out, "status"), 1L)
  expect_true(any(grepl("--embed-full-table must be 'yes' or 'no'", out, fixed = TRUE)))
  expect_false(dir.exists(file.path(getwd(), "._render")))
})

test_that("--embed-full-table-max-mb is documented with its default and rejects a non-positive size", {
  out = suppressWarnings(system2("Rscript", c(render_report, "--help"),
                                 stdout = TRUE, stderr = TRUE))
  help = paste(out, collapse = " ")
  expect_true(grepl("--embed-full-table-max-mb", help, fixed = TRUE))
  expect_true(grepl("(default: 5)", gsub("\\s+", " ", help), fixed = TRUE))

  d = tempfile("sample")
  dir.create(d)
  for (bad in c("0", "lots")) {
    out = suppressWarnings(system2("Rscript",
      c(render_report, "--sample-dir", d, "--sex", "male", "--embed-full-table-max-mb", bad),
      stdout = TRUE, stderr = TRUE))
    expect_equal(attr(out, "status"), 1L)
    expect_true(any(grepl("--embed-full-table-max-mb must be a number above 0", out, fixed = TRUE)))
  }
  expect_false(dir.exists(file.path(getwd(), "._render")))
})

test_that("--run-mode rejects anything but paired or tumour_only, before any rendering", {
  d = tempfile("sample")
  dir.create(d)
  out = suppressWarnings(system2("Rscript",
    c(render_report, "--sample-dir", d, "--sex", "male", "--run-mode", "matched"),
    stdout = TRUE, stderr = TRUE))
  expect_equal(attr(out, "status"), 1L)
  expect_true(any(grepl("--run-mode must be 'paired' or 'tumour_only'", out, fixed = TRUE)))
  expect_false(dir.exists(file.path(getwd(), "._render")))
})

test_that("every R source file parses", {
  # bin/render_report.R is sourced by no unit test -- the helpers it calls are tested
  # directly -- so a syntax error in it survives a fully green suite and only surfaces at
  # render time. A v1.5.0 release candidate shipped exactly that: a top-level `if (...)`
  # with `else` on the following line, which R rejects.
  files = c(render_report,
            list.files(file.path(repo_root, "R"), pattern = "[.]R$",
                       recursive = TRUE, full.names = TRUE))
  expect_gt(length(files), 1)
  for (f in files) {
    expect_no_error(parse(f))
  }
})

test_that("--purity-spread-warn is documented with the report's default and validated", {
  out = suppressWarnings(system2("Rscript", c(render_report, "--help"), stdout = TRUE, stderr = TRUE))
  help = paste(out, collapse = " ")
  expect_true(grepl("--purity-spread-warn", help, fixed = TRUE))
  expect_true(grepl(paste0("(default: ", PURITY_SPREAD_WARN, ")"), help, fixed = TRUE))
  d = tempfile("sample")
  dir.create(d)
  for (bad in c("abc", "0", "-0.1")) {
    out = suppressWarnings(system2("Rscript",
      c(render_report, "--sample-dir", d, "--sex", "male", "--purity-spread-warn", shQuote(bad)),
      stdout = TRUE, stderr = TRUE))
    expect_equal(attr(out, "status"), 1L, info = bad)
    expect_true(any(grepl("--purity-spread-warn must be a positive number", out, fixed = TRUE)), info = bad)
  }
  expect_false(dir.exists(file.path(getwd(), "._render")))
})
