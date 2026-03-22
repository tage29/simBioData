test_that("simulate_phylo_traits returns list with tree", {
  sim <- simulate_phylo_traits(seed = 1)

  expect_true(is.list(sim))
  expect_true(all(c("data", "tree") %in% names(sim)))
})

test_that("data frame structure", {
  sim <- simulate_phylo_traits(return_tree = FALSE, seed = 1)

  expect_s3_class(sim, "data.frame")
  expect_true("Species" %in% names(sim))
})

test_that("multiple traits work", {
  sim <- simulate_phylo_traits(n_traits = 3, return_tree = FALSE, seed = 1)

  expect_equal(ncol(sim), 1 + 3)
})

test_that("reproducibility", {
  s1 <- simulate_phylo_traits(seed = 123)
  s2 <- simulate_phylo_traits(seed = 123)

  expect_equal(s1$data, s2$data)
})
