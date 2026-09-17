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
                 "--gene-panel", "--gene-lists-dir", "--output")) {
    expect_true(grepl(flag, help, fixed = TRUE), info = flag)
  }
})
