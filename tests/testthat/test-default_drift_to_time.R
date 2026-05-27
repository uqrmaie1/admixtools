test_that("default_drift_to_time is identity on a vector", {
  x <- c(0.05, 0.10, 0.20)
  expect_equal(default_drift_to_time(x), x)
})

test_that("default_drift_to_time ignores the second argument", {
  x <- c(1, 2, 3)
  expect_equal(
    default_drift_to_time(x, segment_vec = c("A", "B", "C")),
    x
  )
})
