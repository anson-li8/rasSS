test_that("simple test to ensure quiet returns appropriate value and silences output", {
  expect_equal(quiet(2 + 2), 4)
  expect_silent(quiet({ cat("x"); message("y") }))
})