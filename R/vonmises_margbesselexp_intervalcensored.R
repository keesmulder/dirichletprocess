#'Create a von Mises mixing distribution
#'
#'
#'@param priorParameters Prior parameters for the base measure which are, in
#'  order, (R_0, n_0).
#'
#'@return Mixing distribution object
#'@export
vonMisesICMixtureCreate <- function(priorParameters){
  mdobj <- MixingDistribution("vonmises_ic", priorParameters, "nonconjugate")
  return(mdobj)
}


# Helper function, log of a bessel function.
logBesselI <- function(x, nu) log(besselI(x, nu, expon.scaled = TRUE)) + x

# von Mises
dvm <- Vectorize(function(x, mu, kp) {
  logpdf <- kp * cos(x - mu) - log(2*pi) - logBesselI(kp, 0)
  return(exp(logpdf))
})

# Vectorized cdf
pvm <- Vectorize(circular:::pvonmises)


Likelihood.vonmises_ic <- function(mdobj, x, theta) {
  # If the data is not interval-censored, just use the regular density.
  if (identical(x[1], x[2])) {
    prob <- dvm(x[1], theta[[1]], theta[[2]])

    # For interval-censored data, the likelihood is given roughly by the formula
    # F(end | params) - F(start | params) / (end - start).
  } else {
    prob <- suppressWarnings(
      (pvm(x[2], theta[[1]], theta[[2]]) - pvm(x[1], theta[[1]], theta[[2]])) /
        (x[2] - x[1])
      )
  }
  return(prob)
}


# The prior is dgamma(x, shape = (n0 + 1)/2, rate = n0 - R0)
# if (n_0 < 0 | R_0 < -1) stop("Prior parameters out of bounds.")
PriorDraw.vonmises_ic <- function(mdobj, n = 1) {

  priorParameters <- mdobj$priorParameters
  R_0  <- priorParameters[1]
  n_0  <- priorParameters[2]

  mu <- runif(1, 0, 2*pi)
  kp <- rgamma(1, shape = (n_0 + 1)/2, rate = n_0 - R_0)

  theta <- list(mu = array(mu, dim = c(1, 1, n)),
                kp = array(kp, dim = c(1, 1, n)))
  return(theta)
}

PriorDensity.vonmises_ic <- function(mdobj, theta) {

  priorParameters <- mdobj$priorParameters
  R_0  <- priorParameters[1]
  n_0  <- priorParameters[2]
  as.numeric(dgamma(theta[[2]], shape = (n_0 + 1)/2, rate = n_0 - R_0))
}


MhParameterProposal.vonmises_ic <- function(mdObj, old_params) {

  mhStepSize <- mdObj$mhStepSize
  new_params <- old_params
  new_params[[1]] <- abs(old_params[[1]] + mhStepSize[1] * rnorm(1))
  new_params[[2]] <- abs(old_params[[2]] + mhStepSize[2] * rnorm(1))

  return(new_params)
}


