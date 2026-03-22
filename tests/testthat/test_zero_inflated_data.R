test_that("zero inflated data structure", {
  dat <- simulate_zero_inflated_data(seed = 1)

  expect_s3_class(dat, "data.frame")
  expect_true(all(c("Y", "Env", "Treatment", "Group", "ZeroInflated") %in% names(dat)))
})

test_that("returns components", {
  sim <- simulate_zero_inflated_data(return_components = TRUE, seed = 1)

  expect_true(is.list(sim))
  expect_true(all(c("data", "mu", "zi_prob") %in% names(sim)))
})

test_that("reproducibility", {
  d1 <- simulate_zero_inflated_data(seed = 42)
  d2 <- simulate_zero_inflated_data(seed = 42)

  expect_equal(d1, d2)
})

test_that("zero inflation increases zeros", {
  d1 <- simulate_zero_inflated_data(zi_intercept = -5, seed = 1)
  d2 <- simulate_zero_inflated_data(zi_intercept = 2, seed = 1)

  expect_true(sum(d2$Y == 0) > sum(d1$Y == 0))
})
