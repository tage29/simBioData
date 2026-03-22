test_that("simulate_binomial_data individual output", {
  dat <- simulate_binomial_data(seed = 1)

  expect_s3_class(dat, "data.frame")
  expect_true("Outcome" %in% names(dat))
})

test_that("simulate_binomial_data aggregated output", {
  dat <- simulate_binomial_data(aggregated = TRUE, seed = 1)

  expect_true(all(c("Successes", "Failures", "Trials") %in% names(dat)))
})

test_that("reproducibility works", {
  d1 <- simulate_binomial_data(seed = 42)
  d2 <- simulate_binomial_data(seed = 42)

  expect_equal(d1, d2)
})

test_that("treatment effect changes probability", {
  d1 <- simulate_binomial_data(treatment_effect = 0, seed = 1)
  d2 <- simulate_binomial_data(treatment_effect = 2, seed = 1)

  p1 <- mean(d1$Outcome)
  p2 <- mean(d2$Outcome)

  expect_true(p2 > p1)
})
