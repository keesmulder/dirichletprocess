#'Create a von Mises mixing distribution
#'
#'
#'@param priorParameters Prior parameters for the base measure which are, in
#'  order, (mu_0, R_0, c).
#'
#'@return Mixing distribution object
#'@export
vonMisesMixtureCreate <- function(priorParameters){
  mdobj <- MixingDistribution("normal", priorParameters, "conjugate")
  return(mdobj)
}

# Should solve this at some point, I can't very well depend on all my wonky
# packages.
library(flexcircmix)
library(circglmbayes)

# Helper function, log of a bessel function.
logBesselI <- function(x, nu) log(besselI(x, nu, expon.scaled = TRUE)) + x

# von Mises
dvm <- function(x, mu, kp) {
  logpdf <- kp * cos(x - mu) - log(2*pi) - logBesselI(kp, 0)
  return(exp(logpdf))
}



dbesselexp2 <- function(kp, mu, mu_n, R_n, n_n) {
  eta <- n_n
  g <- -R_n * cos(mu - mu_n) / n_n
  dbesselexp(kp, eta, g)
}
rbesselexp2 <- function(kp, mu, mu_n, R_n, n_n) {
  eta <- n_n
  g <- -R_n * cos(mu - mu_n) / n_n
  rbesselexp(kp, eta, g)
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
  mu <- runif(1, -pi, pi)
  for (i in 1:nsamp) {

    kp <- rbesselexp2(n, mu, mu_0, R_0, n_0)
    mu <- vapply(kp, function(kpi) {
      rvmc(1, mu_0, R_0 * kpi)
    }, FUN.VALUE = 0)
  }

  theta <- list(array(mu, dim = c(1, 1, n)),
                array(kp, dim = c(1, 1, n)))
  return(theta)
}

