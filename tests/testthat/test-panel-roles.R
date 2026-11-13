# The optional `role` column, the all-row coordinate vectors, and how both survive the
# execute_params YAML round trip (see load_gene_panel() in R/utils.R).

write_panel = function(lines) {
  f = tempfile(fileext = ".tsv")
  writeLines(lines, f)
  f
}

test_that("role values normalise case-insensitively to oncogene / tsg / both / blank", {
  f = write_panel(c(
    "# reference: hg38",
    "gene\tchrom\tstart\tend\trole",
    "MYC\tchr8\t127735434\t127742951\tOG",
    "BCL2\tchr18\t63123346\t63320128\toncogene",
    "TP53\tchr17\t7661779\t7687538\tTSG",
    "CDKN2A\tchr9\t21967752\t21995301\ttumour suppressor",
    "PTEN\tchr10\t87863113\t87971930\tTumor Suppressor",
    "NOTCH1\tchr9\t136494433\t136546048\tOG/TSG",
    "CREBBP\tchr16\t3725054\t3880120\ttsg/og",
    "EZH2\tchr7\t148807383\t148884321\tBoth",
    "IGH\tchr14\t105586437\t106879844\t"))
  p = load_gene_panel(f, "hg38")
  expect_true(p$has_roles)
  roles = setNames(p$role, p$role_gene)
  expect_equal(unname(roles[c("MYC", "BCL2")]), c("oncogene", "oncogene"))
  expect_equal(unname(roles[c("TP53", "CDKN2A", "PTEN")]), rep("tsg", 3))
  expect_equal(unname(roles[c("NOTCH1", "CREBBP", "EZH2")]), rep("both", 3))
  # A blank role is "no role", so the gene is simply not listed
  expect_false("IGH" %in% p$role_gene)
  expect_equal(length(p$role_gene), length(p$role))

  r = panel_roles(p)
  expect_equal(names(r), c("gene", "role"))
  expect_equal(nrow(r), 8)
})

test_that("an unrecognised role value is an error naming the gene", {
  f = write_panel(c("gene\trole", "MYC\tOG", "TP53\tdriver"))
  expect_error(load_gene_panel(f), "unrecognised role value.*driver.*TP53")
})

test_that("a gene listed twice keeps its first declared role, once", {
  f = write_panel(c("gene\trole", "MYC\t", "MYC\tOG", "MYC\tTSG"))
  p = load_gene_panel(f)
  expect_equal(p$role_gene, "MYC")
  expect_equal(p$role, "oncogene")
})

test_that("a panel without a role column loads exactly as before, with empty role fields", {
  f = write_panel(c("gene", "MYC", "TP53"))
  p = load_gene_panel(f)
  expect_false(p$has_roles)
  expect_length(p$role_gene, 0)
  expect_equal(nrow(panel_roles(p)), 0)
  expect_equal(p$genes, c("MYC", "TP53"))
  expect_false(p$has_coords)
  # Symbol-only: no coordinates on either path
  expect_equal(nrow(panel_gene_intervals(p, "sv")), 0)
  expect_equal(nrow(panel_gene_intervals(p, "all")), 0)
  expect_null(panel_intervals(p))
})

test_that("the all_* vectors carry every row with coordinates, whatever applies_to says", {
  f = write_panel(c(
    "# reference: hg38",
    "gene\tchrom\tstart\tend\tapplies_to",
    "MYC\tchr8\t127735434\t127742951\t",
    "IGH\t14\t105586437\t106879844\tsv",
    "TP53\tchr17\t7661779\t7687538\tsnv",
    "KMT2D\t\t\t\tsnv"))
  p = load_gene_panel(f, "hg38")
  # The SV-scoped interval vectors are unchanged
  expect_equal(p$interval_gene, c("MYC", "IGH"))
  sv = panel_gene_intervals(p, "sv")
  expect_equal(sv$gene, c("MYC", "IGH"))
  expect_equal(names(sv), c("chrom", "start", "end", "gene"))
  # "all" adds the snv-scoped row that carries coordinates; the blank one is skipped
  all = panel_gene_intervals(p, "all")
  expect_equal(all$gene, c("MYC", "IGH", "TP53"))
  expect_equal(all$chrom, c("chr8", "chr14", "chr17"))
  expect_type(all$start, "integer")
  expect_equal(all$end[3], 7687538L)
})

test_that("an snv-only coordinate panel has no SV intervals but still has all-row ones", {
  f = write_panel(c("gene\tchrom\tstart\tend\tapplies_to",
                    "TP53\tchr17\t7661779\t7687538\tsnv"))
  p = load_gene_panel(f)
  expect_false(p$has_coords)
  expect_equal(nrow(panel_gene_intervals(p, "sv")), 0)
  expect_equal(panel_gene_intervals(p, "all")$gene, "TP53")
})

test_that("roles and all-row intervals survive the execute_params YAML round trip", {
  f = write_panel(c(
    "# reference: hg38",
    "gene\tchrom\tstart\tend\trole",
    "MYC\tchr8\t127735434\t127742951\tOG",
    "TP53\tchr17\t7661779\t7687538\tTSG"))
  p  = load_gene_panel(f, "hg38")
  rt = yaml::yaml.load(yaml::as.yaml(list(all_panels = list(test = p))))$all_panels$test
  expect_equal(panel_roles(rt), panel_roles(p))
  expect_equal(panel_gene_intervals(rt, "all"), panel_gene_intervals(p, "all"))
  expect_equal(panel_gene_intervals(rt, "sv"), panel_gene_intervals(p, "sv"))

  # Empty vectors come back as list(), which the accessors read as "none"
  q  = load_gene_panel(write_panel(c("gene", "MYC")))
  rq = yaml::yaml.load(yaml::as.yaml(list(q = q)))$q
  expect_equal(nrow(panel_roles(rq)), 0)
  expect_equal(nrow(panel_gene_intervals(rq, "all")), 0)
})

test_that("the bundled panels still load for both references", {
  for (ref in c("hg38", "t2t")) {
    ps = load_all_gene_panels(gene_lists_root, ref)
    expect_true(all(c("lymphoid", "sarcoma") %in% names(ps)))
    for (p in ps) {
      all = panel_gene_intervals(p, "all")
      expect_gte(nrow(all), nrow(panel_gene_intervals(p, "sv")))
    }
  }
})
