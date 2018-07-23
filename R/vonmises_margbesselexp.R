#'Create a von Mises mixing distribution
#'
#'
#'@param priorParameters Prior parameters for the base measure which are, in
#'  order, (R_0, c).
#'
#'@return Mixing distribution object
#'@export
vonMisesMixtureCreate <- function(priorParameters){
  mdobj <- MixingDistribution("normal", priorParameters, "nonconjugate")
  return(mdobj)
}


# Helper function, log of a bessel function.
logBesselI <- function(x, nu) log(besselI(x, nu, expon.scaled = TRUE)) + x

# von Mises
dvm <- Vectorize(function(x, mu, kp) {
  logpdf <- kp * cos(x - mu) - log(2*pi) - logBesselI(kp, 0)
  return(exp(logpdf))
})

Likelihood.vonmises <- function(mdobj, x, theta) dvm(x, theta[[1]], theta[[2]])


# The prior is dgamma(x, shape = (n0 + 1)/2, rate = n0 - R0)
# if (n_0 < 0 | R_0 < -1) stop("Prior parameters out of bounds.")
PriorDraw.vonmises <- function(mdobj, n = 1) {

  priorParameters <- mdobj$priorParameters
  R_0  <- priorParameters[1]
  n_0  <- priorParameters[2]

  mu <- runif(1, 0, 2*pi)
  kp <- rgamma(1, shape = (n_0 + 1)/2, rate = n_0 - R_0)

  theta <- list(mu = array(mu, dim = c(1, 1, n)),
                kp = array(kp, dim = c(1, 1, n)))
  return(theta)
}

PriorDensity.vonmises <- function(mdobj, theta) {

  priorParameters <- mdobj$priorParameters
  R_0  <- priorParameters[1]
  n_0  <- priorParameters[2]
  as.numeric(dgamma(theta[[2]], shape = (n_0 + 1)/2, rate = n_0 - R_0))
}



MhParameterProposal.vonmises <- function(mdObj, old_params) {

  mhStepSize <- mdObj$mhStepSize
  new_params <- old_params
  new_params[[1]] <- abs(old_params[[1]] + mhStepSize[1] * rnorm(1))
  new_params[[2]] <- abs(old_params[[2]] + mhStepSize[2] * rnorm(1))

  return(new_params)
}


