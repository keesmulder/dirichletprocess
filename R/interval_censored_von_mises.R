#' Check if a value is in a circular interval.
#'
#' @param x  Value to check from 0 to 2pi.
#' @param lb Lower bound from 0 to 2pi.
#' @param ub Upper bound from 0 to 2pi.
#'
#' @return Logical; whether x lies between lb and ub on the cirle.
#'
#' @examples
#' circ_in_interval(3, 2, 5)
#' circ_in_interval(3, 6, 4)
#' circ_in_interval(3, 6, 2)
#'
circ_in_interval <- function(x, lb, ub) {
  (lb < ub && (x > lb & x < ub)) || (lb > ub && (x > lb || x < ub))
}


#' Create a von Mises interval censored mixing distribution
#'
#'@param priorParameters Prior parameters for the base measure which are, in
#'  order, (mu_0, R_0, c).
#'
#'@return Mixing distribution object
#'@export
vonMisesICMixtureCreate <- function(priorParameters){
  mdobj <- MixingDistribution("vonmises", priorParameters, "ic_conjugate")
  return(mdobj)
}

#' Sample from the von Mises distribution with interval constraints.
#'
#' @param n Number of required samples.
#' @param mu Mean direction of the von Mises distribution.
#' @param kp Concentration parameter of the von Mises distribution.
#' @param lb Lower bound.
#' @param ub Upper bound.
#'
#' @return A numeric vector of sampled values.
#' @export
#'
#' @examples
#' x <- rvm_ic_reject(10000, 0, 4, 5, 1)
#' hist(x, breaks = 100)
rvm_ic_reject <- function(n, mu = 0, kp = 1, lb = 0, ub = 2*pi, max_attempts = 1000) {

  lb <- lb %% (2*pi)
  ub <- ub %% (2*pi)


  replicate(n, {
    attempts <- 0

    while (TRUE) {

      # Sample proposal.
      x_prop <- rvmc(1, mu, kp) %% (2*pi)

      # Check if it satisfies interval constraint.
      if (circ_in_interval(x_prop, lb, ub)) break

      if (attempts >= max_attempts) {
        x_prop <- NA
        break
      }
    }
    x_prop
  })
}




rvm_ic_envelope <- function(n, mu = 0, kp = 1, lb = 0, ub = 2*pi, max_pdf, max_attempts = 1000) {


  lb_pdf <- dvm(lb, mu, kp)
  ub_pdf <- dvm(ub, mu, kp)
  mu_pdf <- dvm(0, 0, kp)


  mu_in_interval <- circ_in_interval(mu, lb, ub)

  # If mu is in the interval, it is the highest point. Otherwise, either the
  # lower bound or the upper bound is the maximum.
  if (mu_in_interval) {
    max_pdf <- mu_pdf
  } else {
    max_pdf <- max(lb_pdf, ub_pdf)
  }

  lb <- lb %% (2*pi)
  ub <- ub %% (2*pi)


  replicate(n, {
    attempts <- 0

    while (TRUE) {
      attempts <- attempts + 1

      x_prop <- runif(1, lb, ifelse(lb > ub, ub + 2*pi, ub))

      prop_prob <- dvm(x_prop, mu, kp)

      if (max_pdf * runif(1) < prop_prob) {
        break
      }
      if (attempts >= max_attempts) {
        x_prop <- NA
        break
      }
    }
    x_prop
  }) %% (2*pi)
}

# Von Mises cdf approximation by normal.
pvm_normal_approx <- function(mu, kp, lb, ub) {

  # Ensure the right format for the following calculations.
  rot_lb <- (lb - mu + pi) %% (2*pi)
  rot_ub <- (ub - mu + pi) %% (2*pi)

  # Estimate of the acceptance rate for the rejection algorithm.
  if (rot_lb < rot_ub) {
    est_acc <- pnorm(rot_ub, mean = pi, sd = 1/sqrt(kp)) - pnorm(rot_lb, mean = pi, sd = 1/sqrt(kp))
  } else {
    est_acc <- 1 - pnorm(rot_ub,  mean = pi, sd = 1/sqrt(kp)) - pnorm(rot_lb,  mean = pi, sd = 1/sqrt(kp))
  }
  est_acc
}


#' Sample from the von Mises distribution with interval constraints.
#'
#' @param n Number of required samples.
#' @param mu Mean direction of the von Mises distribution.
#' @param kp Concentration parameter of the von Mises distribution.
#' @param lb Lower bound.
#' @param ub Upper bound.
#'
#' @return A numeric vector of sampled values.
#' @export
#'
#' @examples
#' x <- rvm_ic(10000, 0, 4, 5, 1)
#' hist(x, breaks = 100)
rvm_ic <- function(n, mu = 0, kp = 1, lb = 0, ub = 2*pi,
                   method = "adaptive", adaptive_cutoff = .1,
                   max_attempts = 1000) {

  if (method == "reject") {
    rvm_ic_reject(n, mu, kp, lb, ub, max_attempts)
  } else if (method == "envelope") {
    rvm_ic_envelope(n, mu, kp, lb, ub, max_attempts)
  } else if (method == "adaptive") {

    # Estimated acceptance rate of the rejection algorithm.
    estimated_acceptance <- pvm_normal_approx(mu, kp, lb, ub)

    # If the acceptance rate is too low,
    if (estimated_acceptance < adaptive_cutoff) {
      rvm_ic_envelope(n, mu, kp, lb, ub, max_attempts)
    } else {
      rvm_ic_reject(n, mu, kp, lb, ub, max_attempts)
    }
  }
}




IntervalCensoredDraw <- function(dpobj) UseMethod("IntervalCensoredDraw", dpobj)

IntervalCensoredDraw.vonmises <- function(dpobj,
                                          method = "adaptive",
                                          adaptive_cutoff = .1,
                                          max_attempts = 10000) {

  y <- dpobj$data
  y_imp <- y[, 1]

  clusterLabels <- dpobj$clusterLabels
  clusterParams <- dpobj$clusterParameters

  cens_idx   <- which(y[, 1] != y[,2])
  cens_cl    <- clusterLabels[cens_idx]

  for (cens_yi in cens_idx) {

    clust_i <- clusterLabels[cens_yi]
    y_imp[cens_yi] <- rvm_ic(1,
                             mu = clusterParams[[1]][, , clust_i, drop = FALSE],
                             kp = clusterParams[[2]][, , clust_i, drop = FALSE],
                             lb = y[cens_yi, 1], ub = y[cens_yi, 2],
                             method = method, adaptive_cutoff = adaptive_cutoff,
                             max_attempts = max_attempts)
  }

  return(matrix(y_imp))
}



