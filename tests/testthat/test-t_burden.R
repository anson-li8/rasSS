test_that("t_burden matches manual calculation", {
  w <- c(1, 2)
  Z <- c(0.5, -1)
  R <- matrix(c(1, 0.2, 0.2, 1), 2, 2)
  # w'Rw = 5.8, w'Z = -1.5 (manually calculated)
  expect_equal(t_burden(w, Z, R), -1.5 / sqrt(5.8), tolerance = 1e-12)
})

test_that("t_burden reduces to Z when w and R are trivial", {
  expect_equal(t_burden(1, 2, matrix(1)), 2)
})

test_that("t_burden is safe against invalid denominators", {
  expect_equal(t_burden(c(0, 0), c(1, 1), diag(2)), 0)
})