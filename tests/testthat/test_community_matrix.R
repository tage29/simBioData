test_that("community matrix default output", {
  mat <- simulate_community_matrix(seed = 1)

  expect_true(is.matrix(mat))
  expect_true(nrow(mat) > 0)
})

test_that("long format output", {
  df <- simulate_community_matrix(return_long = TRUE, seed = 1)

  expect_s3_class(df, "data.frame")
  expect_true(all(c("Site", "Species", "Abundance") %in% names(df)))
})

test_that("reproducibility", {
  m1 <- simulate_community_matrix(seed = 42)
  m2 <- simulate_community_matrix(seed = 42)

  expect_equal(m1, m2)
})

test_that("zero inflation increases zeros", {
  m1 <- simulate_community_matrix(zero_inflation = 0, seed = 1)
  m2 <- simulate_community_matrix(zero_inflation = 0.8, seed = 1)

  expect_true(sum(m2 == 0) > sum(m1 == 0))
})
