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

