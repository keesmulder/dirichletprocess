context("Von Mises Interval Censored Non Conjugate Functions")





test_that("Interval Censored works", {

  skip("Dev tests")

  vmicMd <- MixingDistribution(distribution = "vonmises_ic",
                             priorParameters = c(.9, 1),
                             conjugate = "nonconjugate",
                             mhStepSize = c(.5, .5))


  y <- c(rnorm(40, -2, .2), rnorm(40, 0, .5), rnorm(40, 1.2, .1)) %% (2*pi)

  y_cens <- cbind(y, y + c(0,1)) %% (2 * pi)

  hist(y, breaks = 60)

  dp <- DirichletProcessCreate(y_cens, vmicMd)
  dp <- Initialise(dp)
  res <- Fit(dp, 1000)



  plot(res)

  library(ggplot2)


  # Marginalized kappa prior based on the conjugate where mu_0 is integrated out
  unnorm_marg_kp_prior <- Vectorize(function(kp, R0 = sqrt(1/2), n0 = 1) {
    exp(logBesselI(R0 * kp, 0) - n0 * logBesselI(kp, 0))
  })

  marg_kp_prior <- function(kp, R0 = sqrt(1/2), n0 = 1) {
    nc <- integrate(unnorm_marg_kp_prior, 0, Inf, R0 = R0, n0 = n0)$value
    unnorm_marg_kp_prior(kp, R0, n0) / nc
  }


  # Approximation through gamma function
  marg_kp_approx <- function(kp, R0 = sqrt(1/2), n0 = 1) {
    dgamma(kp, shape = (n0 + 1) / 2, rate = n0 - R0)
  }

  comparison_plot <- function(R0 = sqrt(1/2), n0 = 1, xmax = 15) {
    this_args <- list(R0 = R0, n0 = n0)
    ggplot(data.frame(x = c(0, xmax)), aes(x)) +
      stat_function(fun = marg_kp_prior, col = "darkolivegreen", args = this_args) +
      stat_function(fun = marg_kp_approx, col = "tomato", args = this_args)
  }

  # These all do well
  comparison_plot()
  comparison_plot(0, 1)
  comparison_plot(.95, 1, 100)
  comparison_plot(8, 10)

  # Unrealistic priors for DP models do a bit worse, which is not a problem
  comparison_plot(.08, 1)
  comparison_plot(80, 100)

})
