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


dvm_approx <- function(x, mu= 0, kp = 1) sqrt(kp) * exp(kp * (cos(x - mu) - 1)) / sqrt(2*pi)
dvm_taylor <- function(x, mu= 0, kp = 1) {
  (exp(kp) -
     exp(kp) * kp * (x - mu)^2 / 2 +
     exp(kp) * kp * (3*kp + 1) * (x - mu)^4 / 24 -
     exp(kp) * kp * (15*kp^2 + 15*kp + 1) * (x - mu)^6 / 720

  ) * exp(-kp) * sqrt(kp) / sqrt(2*pi)
}

dvm_taylor2 <- function(x, mu= 0, kp = 1) {
                                         * exp(-kp) * sqrt(kp) / sqrt(2*pi)
}

curve(dvm(x, 3, 1), 0, 2*pi, ylim = c(0, 1))
curve(dvm_approx(x, 3, 1), 0, 2*pi, add = TRUE, col = "blue")
curve(dvm_taylor(x, 3, 1), 0, 2*pi, add = TRUE, col = "red")

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

  attempts <- 0

  replicate(n, {

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
rvm_ic <- function(n, mu = 0, kp = 1, lb = 0, ub = 2*pi, method = "reject", max_attempts = 1000) {

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

  if (method == "reject") {

  } else if (method == "envelope") {

  } else if (method == "adaptive") {
    env_sample <- rvm_ic_envelope(n, mu, kp, lb, ub, max_pdf)
  }

}




IntervalCensoredDraw <- function(dpobj) UseMethod("IntervalCensoredDraw", dpobj)

IntervalCensoredDraw.vonMises <- function(dpobj) {

  y <- dpobj$data
  y_imp <- y[, 1]

  clusterLabels <- dpobj$clusterLabels
  clusterParams <- dpobj$clusterParameters

  cens_idx   <- which(y[, 1] != y[,2])
  cens_cl    <- clusterLabels[cens_idx]
  cens_clust <- unique(cens_cl)


  for (clust_i in cens_clust) {
    y_this_clust <- which(clusterLabels == clust_i)
    y_imp[y_this_clust] <- rvmc(length(y_this_clust),
                                clusterParams[[1]][, , clust_i, drop = FALSE],
                                clusterParams[[2]][, , clust_i, drop = FALSE]
    )
  }

  return(matrix(y_imp))
}



