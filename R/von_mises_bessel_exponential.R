#'Create a Dirichlet Mixture of Von Mises distributions
#'
#'This is the constructor function to produce a \code{dirichletprocess} object
#'with a Von Mises mixture kernel with unknown mean and unknown concentration
#'kappa. The base measure is conjugate to the posterior distribution.
#'
#'
#'The base measure is \eqn{G_0 (\mu, \kappa \mid \gamma) = I_0(R \kappa)^{- n_0}
#'\exp(R_0 \kappa \cos(\mu - \mu_0))}.
#'
#'
#'@param y Data
#'@param g0Priors Base Distribution Priors \eqn{\gamma = (\mu_0, R_0, n_0)}
#'@param alphaPriors Alpha prior parameters. See \code{\link{UpdateAlpha}}.
#'@return Dirichlet process object
#'@export
DirichletProcessVonMises <- function(y, g0Priors = c(0, 1, 1),
                                     alphaPriors = c(2, 4)) {

  mdobj <- vonMisesMixtureCreate(g0Priors)
  dpobj <- DirichletProcessCreate(y, mdobj, alphaPriors)
  dpobj <- Initialise(dpobj)
  return(dpobj)
}


#'Create a von Mises mixing distribution
#'
#'@param priorParameters Prior parameters for the base measure which are, in
#'  order, (mu_0, R_0, c).
#'@param muMargMethod Character; Method to deal with with marginalizing out the
#'  prior mean mu_0. If \code{"marginal"}, the posterior parameters are computed
#'  averaged over a uniform distribution on the prior mean. If \code{"sample"},
#'  a prior mean is sampled uniformly at each iteration.
#'@param n_samp Number of samples to generate from the Gibbs sampler before
#'  accepting the sample as an i.i.d. sample from the posterior. Higher values
#'  mean we can be more certain that the sampler works correctly, at the cost of
#'  computational time. The sampler mixes very fast, but different values can be
#'  attempted by means of sensitivity.
#'
#'@return Mixing distribution object
#'@export
vonMisesMixtureCreate <- function(priorParameters,
                                  muMargMethod = "marginal",
                                  n_samp = 3) {
  mdobj <- MixingDistribution("vonmises", priorParameters, "conjugate")

  mdobj$n_samp <- n_samp

  # If the prior mean direction mu_0 should be treated as unknown, add the
  # method to deal with this to the mixing distribution object.
  if (is.na(priorParameters[1])) {
    if (muMargMethod == "sample") {
      mdobj$muMargSample <- TRUE
    } else if (muMargMethod == "marginal") {
      mdobj$muMargSample <- FALSE
    } else {
      stop(paste("Unknown method for marginalizing out the prior",
                 "mean direction. Select 'sample' or 'marginal'."))
    }
  }

  return(mdobj)
}

# Should solve this at some point, I can't very well depend on all my wonky
# packages.
library(flexcircmix)
library(circglmbayes)

# Helper function, log of a bessel function.
logBesselI <- function(x, nu) log(besselI(x, nu, expon.scaled = TRUE)) + x

# von Mises
dvm <- Vectorize(function(x, mu, kp) {
  logpdf <- kp * cos(x - mu) - log(2*pi) - logBesselI(kp, 0)
  return(exp(logpdf))
})



dbesselexp2 <- function(kp, mu, mu_n, R_n, n_n) {
  eta <- n_n
  g <- -R_n * cos(mu - mu_n) / n_n
  dbesselexp(kp, eta, g)
}
rbesselexp2 <- function(n, mu, mu_n, R_n, n_n) {
  eta <- n_n
  g <- -R_n * cos(mu - mu_n) / n_n
  rbesselexp(n, eta, g)
}



PosteriorParameters.vonmises <- function(mdobj, x) {
  priorParameters <- mdobj$priorParameters

  # If the mean direction is given, compute the posterior parameters as usual.
  # Else, if mu_0 is NA, we want to marginalize over it.
  if (!is.na(priorParameters[1])) {
    mu_0 <- priorParameters[1]

  } else if (mdobj$muMargSample) {

    # If we wish to sample it, sample a uniform prior mean and continue with the
    # usual computation.
    mu_0 <- runif(1, 0, 2*pi)

  } else {
      # When we don't sample the mean direction, we wish to marginalize by
      # taking the expectation over the uniform distribution on mu_0 for R_n.
      R_0  <- priorParameters[2]
      n_0  <- priorParameters[3]

      # Approximate the correct R_n with the complete elliptic integral of the
      # second kind.
      C <- sum(cos(x))
      S <- sum(sin(x))
      R <- sqrt(C^2 + S^2)
      R_n <- (2/pi) * (R_0 + R) * pracma::ellipke(4 * R_0 * R / (R_0 + R)^2)$e

      n_n <- n_0 + length(x)

      # In this case, mu_n is the mean direction of the data.
      mu_n <- atan2(S, C)
      PosteriorParameters <- matrix(c(mu_n, R_n, n_n), ncol = 3)
      return(PosteriorParameters)

  }

  # If we don't marginalize or have sampled mu_0, compute the posterior
  # parameters in the usual way for conjugate von Mises models.
  R_0  <- priorParameters[2]
  n_0  <- priorParameters[3]

  n_n <- n_0 + length(x)
  C_n <- sum(cos(x)) + R_0 * cos(mu_0)
  S_n <- sum(sin(x)) + R_0 * sin(mu_0)

  R_n  <- sqrt(C_n^2 + S_n^2)
  mu_n <- atan2(S_n, C_n)

  PosteriorParameters <- matrix(c(mu_n, R_n, n_n), ncol = 3)

  return(PosteriorParameters)
}



Likelihood.vonmises <- function(mdobj, x, theta) dvm(x, theta[[1]], theta[[2]])



PriorDraw.vonmises <- function(mdobj, n = 1) {

  priorParameters <- mdobj$priorParameters
  R_0  <- priorParameters[2]
  n_0  <- priorParameters[3]

  # Random starting value.
  mu <- runif(n, -pi, pi)

  if (n_0 < 0 | R_0 < -1) stop("Prior parameters out of bounds.")

    # If mu_0 is NA, sample uniform prior mean mu_0.
  if (is.na(priorParameters[1])) {

    mu_0 <- runif(n, 0, 2*pi)

    for (i in 1:mdobj$n_samp) {

      kp <- vapply(1:n, function(i) {
        rbesselexp2(1, mu[i], mu_0[i], R_0, n_0)
      }, FUN.VALUE = 0)
      mu <- vapply(1:n, function(i) {
        rvmc(1, mu_0[i], R_0 * kp[i])
      }, FUN.VALUE = 0)
    }
  } else {

    # Standard algorithm with given mu_0.
    mu_0 <- priorParameters[1]

    for (i in 1:mdobj$n_samp) {

      kp <- vapply(mu, function(mui) {
        rbesselexp2(1, mui, mu_0, R_0, n_0)
      }, FUN.VALUE = 0)
      mu <- vapply(kp, function(kpi) {
        rvmc(1, mu_0, R_0 * kpi)
      }, FUN.VALUE = 0)
    }
  }


  theta <- list(mu = array(mu %% (2*pi), dim = c(1, 1, n)),
                kp = array(kp, dim = c(1, 1, n)))
  return(theta)
}


PosteriorDraw.vonmises <- function(mdobj, x, n = 1) {

  PosteriorParameters_calc <- PosteriorParameters(mdobj, x)

  mu_n <- PosteriorParameters_calc[1]
  R_n  <- PosteriorParameters_calc[2]
  n_n  <- PosteriorParameters_calc[3]

  # if (n_n < 0 | R_n < -1) {
  #   stop("Posterior parameters out of bounds.")
  # }
  if (!isTRUE(n_n >  0)) n_n <- 0
  if (!isTRUE(R_n > -1)) R_n <- -1

  # Random starting value.
  mu <- runif(n, -pi, pi)
  for (i in 1:mdobj$n_samp) {

    kp <- vapply(mu, function(mui) {
      rbesselexp2(1, mui, mu_n, R_n, n_n)
    }, FUN.VALUE = 0)
    mu <- vapply(kp, function(kpi) {
      rvmc(1, mu_n, R_n * kpi)
    }, FUN.VALUE = 0)
  }

  theta <- list(mu = array(mu, dim = c(1, 1, n)),
                kp = array(kp, dim = c(1, 1, n)))
  return(theta)
}



Predictive.vonmises <- function(mdobj, x) {

  # If uninformative prior, only do this once because the predictive will be the
  # same for all data points.
  # if (mdobj$priorParameters[2] == 0) x <- x[1]

  predictiveArray <- numeric(length(x))
  for (i in seq_along(x)) {

    PosteriorParameters_calc <- PosteriorParameters(mdobj, x[i])
    R_n  <- PosteriorParameters_calc[2]
    n_n  <- PosteriorParameters_calc[3]

    # The posterior predictive density with the mean direction integrated out.
    marglikfun <- function(kp) exp(logBesselI(R_n * kp, 0) -
                                     n_n * logBesselI(kp, 0))

    predictiveArray[i] <- integrate(marglikfun, 0, Inf)$value / (2*pi)
  }

  return(predictiveArray)
}
