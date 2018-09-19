#'Create a von Mises mixing distribution
#'
#'
#'@param priorParameters Prior parameters for the base measure which are, in
#'  order, (mu_0, R_0, c).
#'
#'@return Mixing distribution object
#'@export
vonMisesMixtureCreate <- function(priorParameters){
  mdobj <- MixingDistribution("vonmises", priorParameters, "conjugate")
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

  mu_0 <- priorParameters[1]
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



PriorDraw.vonmises <- function(mdobj, n = 1, nsamp = 3) {

  priorParameters <- mdobj$priorParameters
  mu_0 <- priorParameters[1]
  R_0  <- priorParameters[2]
  n_0  <- priorParameters[3]

  if (n_0 < 0 | R_0 < -1) stop("Prior parameters out of bounds.")

  # Random starting value.
  mu <- runif(n, -pi, pi)
  for (i in 1:nsamp) {

    kp <- vapply(mu, function(mui) {
      rbesselexp2(1, mui, mu_0, R_0, n_0)
    }, FUN.VALUE = 0)
    mu <- vapply(kp, function(kpi) {
      rvmc(1, mu_0, R_0 * kpi)
    }, FUN.VALUE = 0)
  }

  theta <- list(mu = array(mu, dim = c(1, 1, n)),
                kp = array(kp, dim = c(1, 1, n)))
  return(theta)
}


PosteriorDraw.vonmises <- function(mdobj, x, n = 1, nsamp = 3) {

  PosteriorParameters_calc <- PosteriorParameters(mdobj, x)

  mu_n <- PosteriorParameters_calc[1]
  R_n  <- PosteriorParameters_calc[2]
  n_n  <- PosteriorParameters_calc[3]

  if (n_n < 0 | R_n < -1) stop("Posterior parameters out of bounds.")

  # Random starting value.
  mu <- runif(n, -pi, pi)
  for (i in 1:nsamp) {

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
  if (mdobj$priorParameters[2] == 0) x <- x[1]

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
