test_that("simulate_longitudinal_data structure", {
  dat <- simulate_longitudinal_data(seed = 1)

  expect_s3_class(dat, "data.frame")
  expect_named(dat, c("ID", "Group", "Time", "Treatment", "Response"))
})

test_that("reproducibility", {
  d1 <- simulate_longitudinal_data(seed = 99)
  d2 <- simulate_longitudinal_data(seed = 99)

  expect_equal(d1, d2)
})

test_that("correct number of observations", {
  dat <- simulate_longitudinal_data(n_individuals = 10, n_timepoints = 5, seed = 1)

  expect_equal(nrow(dat), 10 * 5)
})

test_that("time effect influences response", {
  dat <- simulate_longitudinal_data(time_slope = 2, seed = 1)

  expect_true(cor(dat$Time, dat$Response) > 0)
})
