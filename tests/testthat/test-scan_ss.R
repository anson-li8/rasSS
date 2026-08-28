test_that("scan_ss returns the exact manually-calculated profile", {
  n <- 5
  b <- c(1, 0, 0, 0, 0)
  z <- c(3, 0, 0, 0, 0)
  out <- scan_ss(b, z, diag(n), rep(TRUE, n),
                 skip1 = 1, skip2 = 1, min_window_size = 1, max_window_size = 1)
  expect_equal(out$x, 1:5)
  # t = 3 in both windows -> p = 2*pnorm(-3)
  expect_equal(out$y[1], -log10(2 * pnorm(-3)), tolerance = 1e-12)
  # all-zero window -> t = 0 -> p = 1 -> y = 0
  expect_equal(out$y[3], 0, tolerance = 1e-12)
})

test_that("scan_ss follows the LD-pruning mask", {
  n <- 3
  b <- rep(1, n); z <- rep(3, n)
  full <- scan_ss(b, z, diag(n), rep(TRUE, n),
                  skip1 = 1, skip2 = 1, min_window_size = 1, max_window_size = 1)
  masked <- scan_ss(b, z, diag(n), c(TRUE, FALSE, TRUE),
                    skip1 = 1, skip2 = 1, min_window_size = 1, max_window_size = 1)
  # expect a shrunken window statistic after pruning
  expect_gt(full$y[2], masked$y[2])
})