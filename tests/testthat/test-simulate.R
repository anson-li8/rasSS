test_that("sim_genotypes returns integer genotypic dosages in valid clamped range", {
  set.seed(1)
  X <- sim_genotypes(20, diag(3))
  expect_equal(dim(X), c(20, 3))
  expect_true(all(X %in% 0:2))
})

test_that("prune_ld eliminates only high-LD redundant SNPs", {
  R <- matrix(c(1, 0.9, 0.1, 0.9, 1, 0.1, 0.1, 0.1, 1), 3, 3)
  expect_equal(prune_ld(R, 0.2), c(TRUE, FALSE, TRUE))
})

test_that("make_pheno is reproducible and accurate type-wise", {
  X <- matrix(rnorm(60), 20, 3)
  pc <- make_pheno(X, rep(0, 3), "continuous", 42)
  pb <- make_pheno(X, rep(0, 3), "binary", 42)
  expect_identical(pc, make_pheno(X, rep(0, 3), "continuous", 42))
  expect_length(pc, 20)
  expect_true(all(pb %in% 0:1))
})