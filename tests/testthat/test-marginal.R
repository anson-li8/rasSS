test_that("get_marginal_stats matches lm() (linear regression) on a single SNP", {
  X <- matrix(c(0, 1, 2, 1), 4, 1)
  y <- 2 * X[, 1] + c(0.1, -0.1, 0.1, -0.1)
  out <- get_marginal_stats(X, y)
  fit <- summary(lm(y ~ X))
  expect_equal(out$beta[1], fit$coefficients[2, 1], tolerance = 1e-8)
  expect_equal(out$z[1], fit$coefficients[2, 3], tolerance = 1e-6)
})

test_that("get_marginal_stats_bin is equal to manually-computed Rao score", {
  X <- matrix(c(0, 1, 0, 1), 4, 1)
  y <- c(1, 0, 1, 1)
  out <- get_marginal_stats_bin(X, y)
  # p0=0.75, U=-0.5, V=0.1875 (manually calculated)
  expect_equal(out$beta[1], -0.5 / 0.1875, tolerance = 1e-8)
  expect_equal(out$z[1], -0.5 / sqrt(0.1875), tolerance = 1e-8)
})