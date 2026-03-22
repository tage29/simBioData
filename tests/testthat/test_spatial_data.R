test_that("spatial data structure", {
  dat <- simulate_spatial_data(seed = 1)

  expect_s3_class(dat, "data.frame")
  expect_named(dat, c("X", "Y", "Env", "Treatment", "Group", "Response"))
})

test_that("returns components when requested", {
  sim <- simulate_spatial_data(return_components = TRUE, seed = 1)

  expect_true(is.list(sim))
  expect_true(all(c("data", "spatial_effect", "covariance_matrix") %in% names(sim)))
})

test_that("reproducibility", {
  d1 <- simulate_spatial_data(seed = 123)
  d2 <- simulate_spatial_data(seed = 123)

  expect_equal(d1, d2)
})

test_that("spatial effect induces variation", {
  dat <- simulate_spatial_data(spatial_sd = 5, seed = 1)

  expect_true(sd(dat$Response) > 1)
})
