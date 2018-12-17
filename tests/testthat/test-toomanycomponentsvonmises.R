context("Too many components in von Mises models")

test_that("Predictive", {

  skip("Dev tests")

  mdobj <- MixingDistribution("normal", c(0, 1, 1, 1), "conjugate")
  dirichletprocess:::Predictive.normal(mdobj, 1:10)


  dpexp <- DirichletProcessExponential(1:10)
  dirichletprocess:::Predictive.exponential(dpexp$mixingDistribution, 1:10)


  # Is the predictive always the same?
  dp <- DirichletProcessVonMises(rnorm(100))
  dirichletprocess:::Predictive.vonmises(dp$mixingDistribution, 1:10)

  csq <- seq(0, 2*pi, length.out = 10)
  plot(csq, dirichletprocess:::Predictive.vonmises(DirichletProcessVonMises(
    rnorm(100), g0Priors = c(3, 5, 10))$mixingDistribution, csq))
  plot(csq, dirichletprocess:::Predictive.vonmises(DirichletProcessVonMises(
    rnorm(100), g0Priors = c(1, 0, 10))$mixingDistribution, csq))
  plot(csq, dirichletprocess:::Predictive.vonmises(DirichletProcessVonMises(
    rnorm(100), g0Priors = c(1, 5, 10))$mixingDistribution, csq))
  plot(csq, dirichletprocess:::Predictive.vonmises(DirichletProcessVonMises(
    rnorm(100), g0Priors = c(1, 0, 1))$mixingDistribution, csq))
  plot(csq, dirichletprocess:::Predictive.vonmises(DirichletProcessVonMises(
    rnorm(100), g0Priors = c(NA, 0, 1), priorMeanMethod = "datadependent")$mixingDistribution, csq))
  plot(csq, dirichletprocess:::Predictive.vonmises(DirichletProcessVonMises(
    rnorm(100), g0Priors = c(NA, 0, 1), priorMeanMethod = "integrate")$mixingDistribution, csq))
  plot(csq, dirichletprocess:::Predictive.vonmises(DirichletProcessVonMises(
    rnorm(100), g0Priors = c(NA, 0, 1), priorMeanMethod = "sample")$mixingDistribution, csq))
  plot(csq, dirichletprocess:::Predictive.vonmises(DirichletProcessVonMises(
    rnorm(100), g0Priors = c(NA, 0.8, 1), priorMeanMethod = "datadependent")$mixingDistribution, csq))
  plot(csq, dirichletprocess:::Predictive.vonmises(DirichletProcessVonMises(
    rnorm(100), g0Priors = c(NA, 0.8, 1), priorMeanMethod = "integrate")$mixingDistribution, csq))
  plot(csq, dirichletprocess:::Predictive.vonmises(DirichletProcessVonMises(
    rnorm(100), g0Priors = c(NA, 0.8, 1), priorMeanMethod = "sample")$mixingDistribution, csq))




  # In what way does the sampler get too many componets with simple data?
  th <- c(circglmbayes::rvmc(200, 3, 5), circglmbayes::rvmc(200, 5, 5))

  dp <- DirichletProcessVonMises(th)

  dp_fit <- Fit(dp, 1000)

  DiagnosticPlots(dp_fit)
  plot(dp, data_method = "hist")

})
