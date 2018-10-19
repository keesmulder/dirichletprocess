PosteriorPredictiveDraw <- function(x, ...) {
  UseMethod("PosteriorPredictiveDraw", x)
}


#' Draw from the posterior predictive distribution of a dirichlet process model.
#'
#' By drawing from the posterior predictive distribution, we generate possible
#' future observations, taking into account the uncertainty in the parameters.
#' Thus, this function is useful for prediction as well as decision making.
#'
#' @param dpobj A dirichletprocess object.
#'
#' @return A numeric vector of samples from the posterior predictive
#'   distribution.
#' @export
#'
#' @examples
#' dp <- Fit(DirichletProcessGaussian(c(rnorm(100), rnorm(100, 3))), 10)
#' PosteriorPredictiveDraw(dp, 10)
#'
PosteriorPredictiveDraw.dirichletprocess <- function(dpobj, n = 1) {

  chainlen <- length(dpobj$likelihoodChain)

  # Each time, sample an iteration of the posterior sample, sample a cluster
  # according to the stick-breaking weights, then sample a value from that
  # cluster.
  replicate(n, {

    sampled_chain   <- sample(1:chainlen, size = 1)
    cluster_weights <- dpobj$weightsChain[[sampled_chain]]
    number_clusters <- length(cluster_weights)
    sampled_cluster <- sample(1:number_clusters, size = 1, prob = cluster_weights)
    sampled_cluster

    chosen_params <- vapply(dpobj$clusterParametersChain[[sampled_chain]],
                            function(parray) {
                              parray[1, 1, sampled_cluster]
                            }, FUN.VALUE = 0)

    DrawFromDataDistribution(dpobj, chosen_params, n = 1)
  })
}

#' Draw data according to the distribution of a dirichletprocess
#'
#' This function is simply an S3 wrapper around the \code{r____()} function of
#' each mixing distribution. This construction means that
#' \code{DrawFromDataDistribution.____} can be easily added for each mixing
#' distribution, so that we can use \code{PosteriorPredictiveDraw}, for example.
#'
#' @param dpobj A dirichletprocess object.
#' @param params A parameter vector.
#' @param n Number of required samples
#' @param ...
#'
#' @return n sampled values.
#'
#' @examples
#' dp <- Fit(DirichletProcessGaussian(c(rnorm(100), rnorm(100, 3))), 10)
#' hist(DrawFromDataDistribution(dp, c(1, 5), 1000))
#'
DrawFromDataDistribution <- function(dpobj, ...) {
  UseMethod("DrawFromDataDistribution", dpobj$mixingDistribution)
}

DrawFromDataDistribution.normal <- function(dpobj, params, n = 1, ...) {
  rnorm(n, params[1], params[2])
}



#' Obtain an estimate of the proportion of future data in a certain interval
#'
#' Using the posterior predictive distribution, we can obtain an estimate of the
#' proportion of data in a certain interval, as well as the uncertainty around
#' this estimate across the posterior.
#'
#' @param dpobj A dirichletprocess object.
#' @param from Start of the interval.
#' @param to End of the interval.
#' @param n_reps Nnumber of posterior samples to use to estimate the interval.
#' @param ci_size Width of the credible interval around the probability
#'   estimate.
#' @param ... further arguments.
#'
#' @return The median probabiliity along with the credible interval.
#' @export
#'
#' @examples
#' dp <- Fit(DirichletProcessGaussian(c(rnorm(100), rnorm(100, 3))), 10)
#' PosteriorPredictiveInterval(dp, 0, 3)
#'
PosteriorPredictiveInterval <- function(dpobj, ...) {
  UseMethod("PosteriorPredictiveInterval", dpobj)
}


PosteriorPredictiveInterval.dirichletprocess <- function(dpobj,
                                                         from, to,
                                                         n_reps = 100,
                                                         ci_size = .05,
                                                         ...) {
  # Chain size.
  nits <- length(dpobj$likelihoodChain)

  prob_sample <- replicate(n_reps, {
    integrate(PosteriorFunction(dpobj, sample(1:nits, 1)), from, to)$value
  })

  # Return the quantile
  quantile(prob_sample, c(ci_size / 2, .5, 1 - ci_size / 2))
}
