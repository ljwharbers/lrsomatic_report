# Fake wakhan/ trees in the two layouts the pipeline links in. Each solution's real directory
# holds its plots, and a rank symlink points at it, as Wakhan itself writes them.
#   0.4.x: <P>_<U>_<C>/<sample>_<P>_<U>_<C>_genome_copynumbers_breakpoints.html, solution_<n> ->
#   0.5.0: solution_<P>_<U>_<C>/integer_profile.html,                           solution_rank_<n> ->
wakhan_solutions_dt = function(repository_names) {
  data.table(solution_rank = seq_along(repository_names),
             ploidy = c(2.1, 3.4)[seq_along(repository_names)],
             cell_purity = c(0.82, 0.55)[seq_along(repository_names)],
             repository_name = repository_names)
}

write_wakhan_tree = function(layout, link = TRUE) {
  w = file.path(tempfile("wakhan"), "wakhan")
  dir.create(w, recursive = TRUE)
  puc = c("2.1_0.82_0.97", "3.4_0.55_0.91")
  for (i in seq_along(puc)) {
    if (layout == "0.4") {
      real = puc[i]
      plots = paste0("S1_", puc[i], c("_genome_copynumbers_breakpoints.html",
                                       "_genome_copynumbers_breakpoints_subclonal.html"))
      rank_link = paste0("solution_", i)
    } else {
      real = paste0("solution_", puc[i])
      plots = c("integer_profile.html", "subclonal_profile.html")
      rank_link = paste0("solution_rank_", i)
    }
    dir.create(file.path(w, real))
    for (p in plots) writeLines("<html></html>", file.path(w, real, p))
    if (link) file.symlink(real, file.path(w, rank_link))
  }
  repo = if (layout == "0.4") puc else paste0("solution_", puc)
  list(dir = w, solutions = wakhan_solutions_dt(repo))
}

test_that("locate_wakhan_cn_plots finds the 0.4.x plots through solution_<rank>", {
  t = write_wakhan_tree("0.4")
  got = locate_wakhan_cn_plots(t$dir, t$solutions)
  expect_length(got, 2)
  expect_equal(vapply(got, `[[`, numeric(1), "rank"), c(1, 2))
  expect_equal(basename(got[[1]]$plot), "S1_2.1_0.82_0.97_genome_copynumbers_breakpoints.html")
  expect_equal(basename(dirname(got[[1]]$plot)), "solution_1")
  expect_equal(basename(got[[2]]$plot), "S1_3.4_0.55_0.91_genome_copynumbers_breakpoints.html")
  expect_equal(got[[2]]$purity, 0.55)
  expect_equal(got[[2]]$ploidy, 3.4)
})

test_that("locate_wakhan_cn_plots finds the 0.5.0 integer_profile.html through solution_rank_<rank>", {
  # The 0.4.x lookup found nothing here and the CN tabs vanished without a notice
  t = write_wakhan_tree("0.5")
  got = locate_wakhan_cn_plots(t$dir, t$solutions)
  expect_length(got, 2)
  expect_equal(basename(got[[1]]$plot), "integer_profile.html")
  expect_equal(basename(dirname(got[[1]]$plot)), "solution_rank_1")
  expect_equal(basename(dirname(got[[2]]$plot)), "solution_rank_2")
  expect_equal(got[[1]]$purity, 0.82)
})

test_that("locate_wakhan_cn_plots falls back to repository_name when the rank links are missing", {
  for (layout in c("0.4", "0.5")) {
    t = write_wakhan_tree(layout, link = FALSE)
    got = locate_wakhan_cn_plots(t$dir, t$solutions)
    expect_length(got, 2)
    expect_equal(basename(dirname(got[[2]]$plot)), t$solutions$repository_name[2])
  }
})

test_that("locate_wakhan_cn_plots never picks a subclonal plot", {
  for (layout in c("0.4", "0.5")) {
    t = write_wakhan_tree(layout)
    plots = vapply(locate_wakhan_cn_plots(t$dir, t$solutions), `[[`, character(1), "plot")
    expect_false(any(grepl("subclonal", plots)))
  }
  # A directory holding only the subclonal plot yields no tab rather than the wrong plot
  t = write_wakhan_tree("0.5")
  unlink(file.path(t$dir, "solution_2.1_0.82_0.97", "integer_profile.html"))
  got = locate_wakhan_cn_plots(t$dir, t$solutions)
  expect_length(got, 1)
  expect_equal(got[[1]]$rank, 2)
})

test_that("locate_wakhan_cn_plots returns nothing for missing input", {
  expect_equal(locate_wakhan_cn_plots(NULL, wakhan_solutions_dt("x")), list())
  expect_equal(locate_wakhan_cn_plots(tempdir(), NULL), list())
  # Solutions listed but no directory for them
  w = tempfile("wakhan")
  dir.create(w)
  expect_equal(locate_wakhan_cn_plots(w, wakhan_solutions_dt("2.1_0.82_0.97")), list())
})
