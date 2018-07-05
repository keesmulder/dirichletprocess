context("Von Mises Mixing Distribution Functions")

data_test <-  rvmc(10, 3, 2)
priorParameters_test = matrix(c(1, 0, 2), ncol = 3)

vonmises_object_test = MixingDistribution("vonmises", priorParameters_test, "conjugate")

test_that("Von Mises Class", {

  expect_is(vonmises_object_test, c("list", "vonmises", "conjugate"))

})

test_that("Von Mises Likelihood", {

  cluster_params_test = list(array(c(1,1), dim = c(1,1,2)),
                             array(c(1,2), dim = c(1,1,2)))

  expect_equal(Likelihood(vonmises_object_test, 0, list(0,1)), dvm(0, 0, 1))
})

test_that("Von Mises Prior Draw", {

  expect_length((PriorDraw(vonmises_object_test, 1)), 2)
  expect_equal(dim(PriorDraw(vonmises_object_test, 10)[[1]]), c(1,1,10))
  expect_equal(dim(PriorDraw(vonmises_object_test, 10)[[2]]), c(1,1,10))
})

test_that("Von Mises Posterior Draw", {

  expect_length((PosteriorDraw(vonmises_object_test, data_test, 1)), 2)
  expect_equal(dim(PosteriorDraw(vonmises_object_test, data_test, 10)[[1]]), c(1,1,10))
  expect_equal(dim(PosteriorDraw(vonmises_object_test, data_test, 10)[[2]]), c(1,1,10))

})

test_that("Von Mises Predictive", {


  pred_array <- Predictive(vonmises_object_test, matrix(data_test, ncol = 1))



  expect_equal(length(Predictive(vonmises_object_test, data_test)), length(data_test))
})

test_that("Development tests", {
  skip("Skip dev tests")

  # PRIOR DRAW

  # General test
  PriorDraw.vonmises(vonmises_object_test, 10)


  # Draw picture of marginal prior for comparison.
  prior_vm <- function(mu, kp, mu_0, R_0, n_0) {
    logprob <- kp * R_0 * cos(mu - mu_0) - n_0 * (log(2 * pi) + logBesselI(kp, 0))
    exp(logprob)
  }
  margprior_kp <- Vectorize(function(kp, mu_0, R_0, n_0) {
    integrate(function(mu) prior_vm(mu, kp, mu_0 = mu_0, R_0 = R_0, n_0 = n_0), 0, 2*pi)$value
  }, "kp")
  margprior_mu <- Vectorize(function(mu, mu_0, R_0, n_0) {
    integrate(function(kp) prior_vm(mu, kp, mu_0 = mu_0, R_0 = R_0, n_0 = n_0), 0, Inf)$value
  }, "mu")

  # Make easy numerically marginalized functions.
  mu_0 <- 1; R_0 <- 1.8; n_0 <- 2
  mykpmarg_unnormalized <- Vectorize(function(kp) margprior_kp(kp, mu_0, R_0, n_0))
  mymumarg_unnormalized <- Vectorize(function(mu) margprior_mu(mu, mu_0, R_0, n_0))

  nc_kp <- integrate(mykpmarg_unnormalized, 0, Inf)$value
  nc_mu <- integrate(mykpmarg_unnormalized, 0, 2*pi)$value

  mykpmarg <- function(kp) mykpmarg_unnormalized(kp) / nc_kp
  mymumarg <- function(mu) mymumarg_unnormalized(mu) / nc_mu

  # Test curves
  curve(mykpmarg, 0, 10)
  curve(mymumarg, 0, 2*pi)

  # Draw for comparison
  n_draws <- 100000
  vonmises_object_test$priorParameters <- matrix(c(mu_0, R_0, n_0), ncol = 3)
  prior_draws <- PriorDraw.vonmises(vonmises_object_test, n_draws, nsamp = 15)

  pd_df <- as.data.frame(cbind(mu = prior_draws$mu[1, 1, ], kp = prior_draws$kp[1, 1, ]))

  ggplot(pd_df, aes(mu)) +
    geom_histogram(aes(y = ..density..), binwidth = .02) +
    stat_function(fun = mymumarg)

  ggplot(pd_df, aes(kp)) +
    geom_histogram(aes(y = ..density..), binwidth = .05) +
    stat_function(fun = mykpmarg) + xlim(0, 5)



  # POSTERIOR DRAW
  PosteriorDraw.vonmises(vonmises_object_test, data_test, n = 10)


  # POSTERIOR PARAMETERS
  x <- 1:3
  PosteriorParameters.vonmises(vonmises_object_test, x)
  PosteriorParameters.vonmises(vonmises_object_test, 1)

  # Works without new data as well.
  expect_equal(vonmises_object_test$priorParameters,
               PosteriorParameters.vonmises(vonmises_object_test, numeric(0)))


  # PREDICTIVE
  Predictive(vonmises_object_test, data_test)


  # REAL DATA
  x <- c(rvmc(30, 0, 2), rvmc(50, 2, 15), rvmc(80, 3.6, 20)) %% (2*pi)
  hist(x, breaks = 100)

  priPar     <- matrix(c(0, 0, 1), ncol = 3)
  vonmisesMD <-  MixingDistribution("vonmises", priPar, "conjugate")

  # curve(dvm(1, priPar[1], priPar[2] * x), 0, 100)

  dpvm     <- DirichletProcessCreate(x, vonmisesMD)
  str(dpvm)

  dpvmintd <- Initialise(dpvm)
  str(dpvmintd)

  dpfit <- Fit(dpvmintd, 10)

  graph <- plot(dpfit, likelihood = TRUE, single = TRUE); graph
})
