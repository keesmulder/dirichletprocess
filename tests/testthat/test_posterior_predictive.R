context("Posterior predictive")

test_that("Sample from posterior predictive distribution", {

  dp <- Fit(DirichletProcessGaussian(rnorm(100)), 10)

  expect_is(DrawFromDataDistribution(dp, c(1, 10)), "numeric")
  expect_is(PosteriorPredictiveDraw(dp), "numeric")
  expect_is(PosteriorPredictiveDraw(dp, 3), "numeric")


})

test_that("Posterior predictive interval", {

  dp <- Fit(DirichletProcessGaussian(rnorm(100)), 10)

  expect_is(PosteriorPredictiveInterval(dp, 0, 1, 4), "numeric")
  expect_is(PosteriorPredictiveInterval(dp, 3, 0, 4), "numeric")
  expect_length(PosteriorPredictiveInterval(dp, 3, 4), 3)


})
