test_that("simulate_seed_data_nb returns correct structure", {
  dat <- simulate_seed_data_nb(seed = 1)

  expect_s3_class(dat, "data.frame")
  expect_true(nrow(dat) > 0)
  expect_named(dat, c("Seed_count", "Elevation", "Pollination", "Species"))
})

test_that("simulate_seed_data_nb reproducibility", {
  d1 <- simulate_seed_data_nb(seed = 29)
  d2 <- simulate_seed_data_nb(seed = 29)

  expect_equal(d1, d2)
})

test_that("treatment effect influences counts", {
  d_low <- simulate_seed_data_nb(treatment_effect = log(0.5), seed = 1)
  d_high <- simulate_seed_data_nb(treatment_effect = log(2), seed = 1)

  m_low <- mean(d_low$Seed_count[d_low$Pollination == "hand"])
  m_high <- mean(d_high$Seed_count[d_high$Pollination == "hand"])

  expect_true(m_high > m_low)
})

test_that("single species works", {
  dat <- simulate_seed_data_nb(species = "sp1", seed = 1)
  expect_true(all(dat$Species == "sp1"))
})
