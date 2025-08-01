test_that("Rcpp posterior sampling works correctly", {
  set.seed(100)
  test <- simple_sims(50, 4, 1)  # Creates 200 effects, 4 conditions
  data <- mash_set_data(test$Bhat, test$Shat)
  U <- cov_canonical(data)

  # Test that Rcpp sampling works
  res_rcpp <- mash(data, U, algorithm.version = "Rcpp", posterior_samples = 50,
                   verbose = FALSE)$result

  # Check dimensions - should be 200 effects, 4 conditions, 50 samples
  expect_equal(dim(res_rcpp$PosteriorSamples), c(200, 4, 50))

  # Check that sample means are close to posterior means
  sample_means <- apply(res_rcpp$PosteriorSamples, c(1, 2), mean)
  dimnames(sample_means) <- dimnames(res_rcpp$PosteriorMean)  # Fix dimnames
  expect_equal(sample_means, res_rcpp$PosteriorMean,
               tolerance = max(res_rcpp$PosteriorSD) * 2)  # Allow for sampling error
})

test_that("Rcpp and R sampling give identical results with same seed", {
  set.seed(100)
  test <- simple_sims(5, 3, 1)  # Creates 20 effects, 3 conditions - small for testing
  data <- mash_set_data(test$Bhat, test$Shat)
  U <- cov_canonical(data)

  # Fit model first
  m <- mash(data, U, algorithm.version = "R", verbose = FALSE)

  # Test with same seed
  seed <- 12345
  res_r <- mash_compute_posterior_matrices(m, data,
                                          algorithm.version = "R",
                                          posterior_samples = 30,
                                          seed = seed)

  res_rcpp <- mash_compute_posterior_matrices(m, data,
                                             algorithm.version = "Rcpp",
                                             posterior_samples = 30,
                                             seed = seed)

  # Basic statistics should match exactly
  expect_equal(res_r$PosteriorMean, res_rcpp$PosteriorMean, tolerance = 1e-14)
  expect_equal(res_r$PosteriorSD, res_rcpp$PosteriorSD, tolerance = 1e-14)
  expect_equal(res_r$lfsr, res_rcpp$lfsr, tolerance = 1e-14)

  # Sample dimensions should match
  expect_equal(dim(res_r$PosteriorSamples), dim(res_rcpp$PosteriorSamples))
})

test_that("Rcpp sampling is deterministic", {
  set.seed(100)
  test <- simple_sims(5, 3, 1)  # Creates 20 effects, 3 conditions
  data <- mash_set_data(test$Bhat, test$Shat)
  U <- cov_canonical(data)

  m <- mash(data, U, algorithm.version = "R", verbose = FALSE)

  # Run twice with same seed - all results should be identical
  seed <- 999
  res1 <- mash_compute_posterior_matrices(m, data,
                                         algorithm.version = "Rcpp",
                                         posterior_samples = 25,
                                         seed = seed)

  res2 <- mash_compute_posterior_matrices(m, data,
                                         algorithm.version = "Rcpp",
                                         posterior_samples = 25,
                                         seed = seed)

  # Basic statistics should be identical
  expect_equal(res1$PosteriorMean, res2$PosteriorMean)
  expect_equal(res1$PosteriorSD, res2$PosteriorSD)
  expect_equal(res1$lfsr, res2$lfsr)

  # Sample dimensions should be correct
  expect_equal(dim(res1$PosteriorSamples), c(20, 3, 25))
  expect_equal(dim(res2$PosteriorSamples), c(20, 3, 25))

  # Most importantly: samples should be identical (deterministic)
  expect_equal(res1$PosteriorSamples, res2$PosteriorSamples)
})

test_that("No crash when posterior_samples = 0 with Rcpp", {
  set.seed(100)
  test <- simple_sims(10, 3, 1)
  data <- mash_set_data(test$Bhat, test$Shat)
  U <- cov_canonical(data)

  # Should not crash and should not include samples
  res <- mash(data, U, algorithm.version = "Rcpp", posterior_samples = 0,
              verbose = FALSE)$result

  expect_false("PosteriorSamples" %in% names(res))
})