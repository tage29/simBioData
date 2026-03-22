test_that("multivariate traits structure", {
  dat <- simulate_multivariate_traits(seed = 1)

  expect_s3_class(dat, "data.frame")
  expect_true(ncol(dat) >= 4) # ID, Group, Env, Treatment + traits
})

test_that("returns covariance when requested", {
  sim <- simulate_multivariate_traits(return_cov = TRUE, seed = 1)

  expect_true(is.list(sim))
  expect_true(all(c("data", "covariance") %in% names(sim)))
})

test_that("reproducibility", {
  d1 <- simulate_multivariate_traits(seed = 123)
  d2 <- simulate_multivariate_traits(seed = 123)

  expect_equal(d1, d2)
})

test_that("treatment effect influences traits", {
  d1 <- simulate_multivariate_traits(treatment_effects = c(0, 0, 0), seed = 1)
  d2 <- simulate_multivariate_traits(treatment_effects = c(2, 2, 2), seed = 1)

  expect_true(mean(d2[, 5:ncol(d2)]) > mean(d1[, 5:ncol(d1)]))
})
