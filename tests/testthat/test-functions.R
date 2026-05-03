library(cellchatr)

# ── rust_subset_genes ────────────────────────────────────────

test_that("rust_subset_genes: normal case returns intersection", {
  expr <- c("TGFB1", "VEGFA", "CD44", "GAPDH")
  db   <- c("TGFB1", "VEGFA", "MYC")
  result <- rust_subset_genes(expr, db)
  expect_equal(sort(result), c("TGFB1", "VEGFA"))
})

test_that("rust_subset_genes: empty expr_genes returns empty", {
  result <- rust_subset_genes(character(0), c("TGFB1", "VEGFA"))
  expect_equal(result, character(0))
})

test_that("rust_subset_genes: empty db_genes returns empty", {
  result <- rust_subset_genes(c("TGFB1", "VEGFA"), character(0))
  expect_equal(result, character(0))
})

test_that("rust_subset_genes: no overlap returns empty", {
  result <- rust_subset_genes(c("GENE1", "GENE2"), c("GENE3", "GENE4"))
  expect_equal(result, character(0))
})

test_that("rust_subset_genes: complete overlap returns all", {
  genes  <- c("TGFB1", "VEGFA", "CD44")
  result <- rust_subset_genes(genes, genes)
  expect_equal(sort(result), sort(genes))
})

test_that("rust_subset_genes: duplicates in expr handled", {
  result <- rust_subset_genes(c("TGFB1", "TGFB1", "VEGFA"), c("TGFB1"))
  expect_true(all(result == "TGFB1"))
})

# ── rust_wilcoxon_filter ─────────────────────────────────────

test_that("rust_wilcoxon_filter: clearly overexpressed genes detected", {
  set.seed(42)
  counts <- matrix(c(
    rep(10, 50), rep(0, 50),  # gene 1 — high in A
    rep(1,  50), rep(1, 50)   # gene 2 — flat
  ), nrow = 2, byrow = TRUE)
  labels     <- c(rep("A", 50), rep("B", 50))
  gene_names <- c("GENE1", "GENE2")
  result <- rust_wilcoxon_filter(counts, labels, gene_names, 0.05)
  expect_true("GENE1" %in% result)
  expect_false("GENE2" %in% result)
})

test_that("rust_wilcoxon_filter: all zeros returns empty", {
  counts     <- matrix(0, nrow = 10, ncol = 50)
  labels     <- c(rep("A", 25), rep("B", 25))
  gene_names <- paste0("GENE", 1:10)
  result <- rust_wilcoxon_filter(counts, labels, gene_names, 0.05)
  expect_equal(result, character(0))
})

test_that("rust_wilcoxon_filter: pval_threshold 0 returns nothing", {
  set.seed(42)
  counts     <- matrix(rnorm(200), nrow = 10)
  labels     <- c(rep("A", 10), rep("B", 10))
  gene_names <- paste0("GENE", 1:10)
  result <- rust_wilcoxon_filter(counts, labels, gene_names, 0.0)
  expect_equal(result, character(0))
})

test_that("rust_wilcoxon_filter: pval_threshold 1 returns all genes", {
  set.seed(42)
  counts     <- matrix(rnorm(1000), nrow = 50)
  labels     <- c(rep("A", 10), rep("B", 10))
  gene_names <- paste0("GENE", 1:50)
  result <- rust_wilcoxon_filter(counts, labels, gene_names, 1.0)
  expect_equal(length(result), 50)
})

# ── rust_match_lr_pairs ──────────────────────────────────────

test_that("rust_match_lr_pairs: correct interactions returned", {
  oe      <- c("TGFB1", "VEGFA")
  ligs    <- c("TGFB1", "MYC",  "VEGFA")
  recs    <- c("TGFBR1","MYCBP","FLT1")
  names   <- c("TGFB1_TGFBR1", "MYC_MYCBP", "VEGFA_FLT1")
  result  <- rust_match_lr_pairs(oe, ligs, recs, names)
  expect_equal(sort(result), c("TGFB1_TGFBR1", "VEGFA_FLT1"))
})

test_that("rust_match_lr_pairs: empty overexpressed returns empty", {
  result <- rust_match_lr_pairs(
    character(0),
    c("TGFB1"), c("TGFBR1"), c("TGFB1_TGFBR1")
  )
  expect_equal(result, character(0))
})

test_that("rust_match_lr_pairs: no matches returns empty", {
  result <- rust_match_lr_pairs(
    c("GENE1"),
    c("TGFB1"), c("TGFBR1"), c("TGFB1_TGFBR1")
  )
  expect_equal(result, character(0))
})

test_that("rust_match_lr_pairs: receptor match also works", {
  result <- rust_match_lr_pairs(
    c("TGFBR1"),          # receptor, not ligand
    c("TGFB1"), c("TGFBR1"), c("TGFB1_TGFBR1")
  )
  expect_equal(result, "TGFB1_TGFBR1")
})

test_that("rust_match_lr_pairs: all interactions match", {
  oe    <- c("TGFB1", "VEGFA")
  ligs  <- c("TGFB1", "VEGFA")
  recs  <- c("TGFBR1", "FLT1")
  names <- c("INT1", "INT2")
  result <- rust_match_lr_pairs(oe, ligs, recs, names)
  expect_equal(length(result), 2)
})
