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




test_that("Constructor function", {

  n <- 100
  y <- rnorm(n) %% (2*pi)

  dp <- DirichletProcessVonMises(y, c(0, 0, 1), c(2, 3))
  expect_s3_class(dp, "dirichletprocess")

  dp <- Fit(dp, 4)

  expect_s3_class(dp, "dirichletprocess")
})





test_that("Marginalized prior sampling/computation", {

  expect_null(vonMisesMixtureCreate(c(0,0,1))$muMargSample)
  expect_true(!is.null(vonMisesMixtureCreate(c(NA,0,1))$muMargSample))

  expect_is(PriorDraw.vonmises(vonMisesMixtureCreate(c(0,0,1))), "list")

  expect_is(PriorDraw.vonmises(vonMisesMixtureCreate(c(NA,0,1))), "list")

  n <- 100
  y <- rnorm(n) %% (2*pi)

  dp <- DirichletProcessVonMises(y, c(NA, 0, 1), c(2, 3))
  expect_s3_class(dp, "dirichletprocess")

  dp <- Fit(dp, 4)

  expect_s3_class(dp, "dirichletprocess")
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

  dirichletprocess:::Predictive.vonmises(vonmisesMD, x)
  # curve(dvm(1, priPar[1], priPar[2] * x), 0, 100)

  dpvm     <- DirichletProcessCreate(x, vonmisesMD)

  dpvmintd <- Initialise(dpvm)

  dpfit <- Fit(dpvmintd, 30)

  graph <- plot(dpfit, likelihood = TRUE, single = TRUE); graph
  # graph
  # plot(density(x))




  # Elliptic integral tests.

  rn <- function(x, C = 4, S = 3, r0 = 3.5) sqrt((C + r0 * cos(x))^2 + (S + r0 * sin(x))^2)

  C = 4; S = 3; r0 = 3.5

  myrn <- function(x) rn(x, C, S, r0)
  curve(myrn, 0, 2*pi)
  integ_approx <- integrate(myrn, 0, 2*pi)$value / (2*pi)

  R <- sqrt(C^2 + S^2)

  ellip_approx <- (2/pi) * (r0 + R) * pracma::ellipke(4 * (r0 * R) / (r0 + R)^2)$e

  # They are equal.
  round(ellip_approx, 6) == round(integ_approx, 6)


  # Which is faster?
  ellip_fun <- function(tol = .00001) (2/pi) * (r0 + R) * pracma::ellipke(4 * (r0 * R) / (r0 + R)^2, tol = tol)$e
  integ_fun <- function() integrate(myrn, 0, 2*pi)$value / (2*pi)
  ellip_fun()
  integ_fun()

  microbenchmark(
    ellip_fun(.00001),
    ellip_fun(.01),
    ellip_fun(.1),
    integ_fun(),
    times = 10000
  )
  # only a little bit faster.


  # Tolerance doesnt matter???
  system.time(replicate(50000, ellip_fun(.00001)))
  system.time(replicate(50000, ellip_fun(.001)))
  system.time(replicate(50000, ellip_fun(.1)))



  # Mu_0 Marginalization
  mu_sam <- as.numeric(PriorDraw.vonmises(vonMisesMixtureCreate(c(NA,.8,1)), 1000)[[1]])
  hist(mu_sam)


  th <- c(0, pi/2)
  th_bar <- pi/4


  # R_n = sqrt(2) but, m = 3. Also mu_n = th_bar
  PosteriorParameters(vonMisesMixtureCreate(c(0, 0, 1)), th)

  # mu_0 doesnt matter if R_0 = 0.  mu_n = th_bar.
  PosteriorParameters(vonMisesMixtureCreate(c(NA,0,1)), th)

  # If R_0 > 0, mu_n != th_bar.
  PosteriorParameters(vonMisesMixtureCreate(c(0,0.5,1)), th)

  # Setting mu_0 == NA thus helps if we have R_0 > 0 to get more concentrated
  # priors with still R_n = th_bar.
  PosteriorParameters(vonMisesMixtureCreate(c(NA, 1,1)), th)


  # This is even more pronounced if a cluster only gets assigned n = 1, which
  # happens often in the DPM model.
  th <- pi/4
  PosteriorParameters(vonMisesMixtureCreate(c(0, 0, 1)), th)

  PosteriorParameters(vonMisesMixtureCreate(c(NA, 1, 1)), th)

  # Can we set R_0 > n_0???
  PosteriorParameters(vonMisesMixtureCreate(c(NA, sqrt(2), 1)), th)
  PosteriorParameters(vonMisesMixtureCreate(c(NA, 2000, 1)), th)
  # This seems dangerous.


  # Mu_0 Marginalization
  hist(PosteriorDraw.vonmises(vonMisesMixtureCreate(c(NA,.8,1)), th, 1000)[[1]])
  hist(PosteriorDraw.vonmises(vonMisesMixtureCreate(c(NA,.8,1), muMargMethod = "sample"), th, 1000)[[1]])

})


