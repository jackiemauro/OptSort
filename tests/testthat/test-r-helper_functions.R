context("test-r-helper_functions")
require(testthat)

test_that("expit no greater than 1", {
  expect_true(expit(100)<=1)
})

test_that('no a throws error', {
  x = data.frame(X1 = rbeta(1000,1,2), X2 = rbeta(1000,1,2))
  as = sample(c(1:4), 1000, replace = T)
  y = rbinom(1000,1,as/5)
  df = data.frame(y = y, A = as.factor(as), x)
  expect_error(optSort(df, '~/temporaryoutputs/', optimizers = c("unconstrained","constrained"), sections = list(c(1,2),c(3,4))))
})

test_that('no y throws error', {
  x = data.frame(X1 = rbeta(1000,1,2), X2 = rbeta(1000,1,2))
  as = sample(c(1:4), 1000, replace = T)
  y = rbinom(1000,1,as/5)
  df = data.frame(Y = y, a = as.factor(as), x)
  expect_error(optSort(df, '~/temporaryoutputs/', optimizers = c("unconstrained","constrained"), sections = list(c(1,2),c(3,4))))
})

test_that('constrained opt works in sections', {
  require(SuperLearner)
  x = data.frame(X1 = rbeta(1000,1,2), X2 = rbeta(1000,1,2))
  as = sample(c(1:4), 1000, replace = T)
  y = rbinom(1000,1,as/5)
  df = data.frame(y = y, a = as.factor(as), x)
  out = optSort(df, '~/temporaryoutputs/', optimizers = c("constrained"), sections = list(c(1,2),c(3,4)))
})

test_that('unconstrained opt works in sections', {
  require(SuperLearner)
  x = data.frame(X1 = rbeta(1000,1,2), X2 = rbeta(1000,1,2))
  as = sample(c(1:4), 1000, replace = T)
  y = rbinom(1000,1,as/5)
  df = data.frame(y = y, a = as.factor(as), x)
  out = optSort(df, '~/temporaryoutputs/', optimizers = c("unconstrained"), sections = list(c(1,2),c(3,4)))
})

test_that('approximate constrained opt works in sections', {
  require(SuperLearner)
  x = data.frame(X1 = rbeta(1000,1,2), X2 = rbeta(1000,1,2))
  as = sample(c(1:4), 1000, replace = T)
  y = rbinom(1000,1,as/5)
  df = data.frame(y = y, a = as.factor(as), x)
  out = optSort(df, '~/temporaryoutputs/', optimizers = c("approximate"), sections = list(c(1,2),c(3,4)))
})

test_that('we get about the result we expect', {
  # mean is .3 if at optimal a, .5 otherwise
  # under constraint, optimal is ~.34
  require(SuperLearner)
  n = 5000
  x = matrix(rnorm(2*n), ncol = 2)
  p1 = 0.5*(x[,1] > 1) + 0.2; p2 = 0.5*(x[,1] < -1) + 0.2; p3 = 0.6*( (x[,1]> -1)*(x[,1] < 1)) + .1
  a = apply(cbind(p1,p2,p3), 1, function(k) sample(c(1:3), 1, replace = F, prob = k))

  p1s = 0.5*(x[,1] > 0) + 0.2; p2s = 0.5*(x[,1] < -1/2) + 0.2; p3s = 0.6*( (x[,1]> -1/2)*(x[,1] < 0)) + .1
  astar = apply(cbind(p1s,p2s,p3s), 1, function(k) sample(c(1:3), 1, replace = F, prob = k))

  mu = runif(n,.25,.75) - 0.2*(a==astar)
  y = rbinom(n,1,mu)
  df = data.frame(y = y, a = as.factor(a), x)
  out = optSort(df, '~/temporaryoutputs/', optimizers = c("unconstrained", "constrained", "approximate"))
  expect_less_than(out$unconstrained$results$Estimate[2],.35)
  expect_less_than(out$constrained$results$Estimate[2],.4)
  expect_less_than(out$approximate$results$Estimate[2],.4)
})
