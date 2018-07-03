context("Von Mises Mixing Distribution Functions")

data_test = rvmc(10, 3, 2)
priorParameters_test = matrix(c(1,1,2), ncol = 3)

vonmises_object_test = MixingDistribution("vonmises", priorParameters_test, "conjugate")

test_that("Von Mises Class", {

  expect_is(vonmises_object_test, c("list", "vonmises", "conjugate"))

})

test_that("Von Mises Likelihood", {

  cluster_params_test = list(array(c(1,1), dim = c(1,1,2)),
                             array(c(1,2), dim = c(1,1,2)))

  expect_equal(Likelihood(vonmises_object_test, 0, list(0,1)), dvm(0, 0, 1))
})

test_that("Von Mises Prior Draw", {

  expect_length((PriorDraw(vonmises_object_test, 1)), 2)
  expect_equal(dim(PriorDraw(vonmises_object_test, 10)[[1]]), c(1,1,10))
  expect_equal(dim(PriorDraw(vonmises_object_test, 10)[[2]]), c(1,1,10))
})

test_that("Von Mises Posterior Draw", {

  expect_length((PosteriorDraw(vonmises_object_test, data_test, 1)), 2)
  expect_equal(dim(PosteriorDraw(vonmises_object_test, data_test, 10)[[1]]), c(1,1,10))
  expect_equal(dim(PosteriorDraw(vonmises_object_test, data_test, 10)[[2]]), c(1,1,10))

})

test_that("Von Mises Predictive", {


  pred_array <- Predictive(vonmises_object_test, matrix(data_test, ncol = 1))



  expect_equal(length(Predictive(vonmises_object_test, data_test)), length(data_test))
})

test_that("Development tests", {
  skip("Skip dev testsS")

  # PRIOR DRAW
  PriorDraw.vonmises(vonmises_object_test, 10)

  # POSTERIOR DRAW


  # POSTERIOR PARAMETERS
  x <- 1:3
  PosteriorParameters.vonmises(vonmises_object_test, x)
  PosteriorParameters.vonmises(vonmises_object_test, 1)

  # Works without new data as well.
  expect_equal(vonmises_object_test$priorParameters,
               PosteriorParameters.vonmises(vonmises_object_test, numeric(0)))


})
