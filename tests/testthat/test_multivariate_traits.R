test_that("treatment effect influences traits", {
  dat <- simulate_multivariate_traits(
    treatment_effects = c(2, 2, 2),
    seed = 1
  )

  # Extract only trait columns safely
  trait_cols <- grep("^Trait_", names(dat), value = TRUE)

  control_mean <- mean(as.matrix(dat[dat$Treatment == "control", trait_cols]), na.rm = TRUE)
  treatment_mean <- mean(as.matrix(dat[dat$Treatment == "treatment", trait_cols]), na.rm = TRUE)

  expect_true(treatment_mean > control_mean)
})
