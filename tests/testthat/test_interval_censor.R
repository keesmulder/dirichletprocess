context("Interval censored models work")


test_that("pvm approximation is accurate", {

  # test pvm approximation
  expect_true(abs(pvm_normal_approx(pi, 2, 4, 5) -
                integrate(dvm, 4, 5, mu = pi, kp = 2)$value) < .1)

  expect_true(abs(pvm_normal_approx(pi, 4, 5, 1) -
                (1 - integrate(dvm, 1, 5, mu = pi, kp = 4)$value)) < .1)

  expect_true(abs(pvm_normal_approx(-pi, 40, 5, 6) -
                (integrate(dvm, 5, 6, mu = -pi, kp = 40)$value)) < .1)

})



test_that("Interval censoring works", {
  x <- rvm_ic(1000, 2, 4, 1, 1.5)

  expect_false(any(x < 1))
  expect_false(any(x > 1.5))

})

test_that("Interval censoring in dp object", {

  n <- 100
  y_ast <- rnorm(n) %% (2*pi)
  #generate sample data
  y <- cbind(y_ast - sample(0:1, n, replace = TRUE) * abs(rnorm(n)),
             y_ast + sample(0:1, n, replace = TRUE) * abs(rnorm(n))) %% (2*pi)

#
#   regvm <- DirichletProcessCreate(y_ast, vonMisesMixtureCreate(c(0, 0, 1)))
#   regfit <- Fit(Initialise(regvm), 100)
#   plot(regfit)

  md <- vonMisesICMixtureCreate(c(0, 0, 1))

  dpvm     <- DirichletProcessCreate(y, md)

  # We can intialize the ic data.
  DrawInitialIntervalCensored(dpvm)


  dpvmintd <- Initialise(dpvm)
  dpfit <- Fit(dpvmintd, 100)


  plot(dpfit)

})
