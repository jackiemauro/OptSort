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

test_that('opt works in sections', {
  require(SuperLearner)
  N = 500
  x = data.frame(X1 = rbeta(N,1,2), X2 = rbeta(N,1,2))
  as = sample(c(1:4), N, replace = T)
  y = rbinom(N,1,as/5)
  df = data.frame(y = y, a = as.factor(as), x)
  out = optSort(df, '~/temporaryoutputs/', optimizers = c("approximate","unconstrained","constrained"), sections = list(c(1,2),c(3,4)))
})

test_that('we get about the result we expect on an easier simulation', {
  #a=1 is optimal for everyone
  #if a=1, mu is .35. otherwise it's .65
  #prob getting a=1 27% so constrained min is 0.65*(1-.27) + 0.35*0.27 = 0.56--essentially you can do no better

  require(SuperLearner)
  n = 5000
  x = matrix(rnorm(2*n), ncol = 2)
  p1 = 0.5*(x[,1] > 1) + 0.2; p2 = 0.5*(x[,1] < -1) + 0.2; p3 = 0.6*( (x[,1]> -1)*(x[,1] < 1)) + .1
  a = apply(cbind(p1,p2,p3), 1, function(k) sample(c(1:3), 1, replace = F, prob = k))

  mu = runif(n,.5,.8) - 0.3*(a==1)
  y = rbinom(n,1,mu)
  df = data.frame(y = y, a = as.factor(a), x)
  rm("a","x","y")

  out = optSort(df, '~/temporaryoutputs/', optimizers = c("unconstrained", "constrained","approximate"))
  #out = optSort(df, '~/temporaryoutputs/', optimizers = c("unconstrained", "constrained", "approximate"))
  expect_less_than(out$unconstrained$results$Estimate[2],.4)
  expect_less_than(out$constrained$results$Estimate[2],.6)
  expect_less_than(out$approximate$results$Estimate[2],.5)
})

test_that('we get about the result we expect on an slightly tougher simulation', {
  #a=1 is optimal for everyone, but best if you have high values of x2
  #if a=1 and x2>0, mu = .25
  #if a=1 and x2<0, mu = .5
  #if a!=1, mu = .75
  #unconstrained: .5*.25 + .5*.5 = .375
  #constrained: top 27% w x2>0 goes to a, gets .25; everyone goes a!=1 gets .75 = .61

  require(SuperLearner)
  n = 5000
  x = matrix(rnorm(2*n), ncol = 2)
  p1 = 0.5*(x[,1] > 1) + 0.2; p2 = 0.5*(x[,1] < -1) + 0.2; p3 = 0.6*( (x[,1]> -1)*(x[,1] < 1)) + .1
  a = apply(cbind(p1,p2,p3), 1, function(k) sample(c(1:3), 1, replace = F, prob = k))

  mu = (runif(n,.3,.7) - 0.25*(x[,2]>0))*(a==1) + runif(n,.55,.95)*(a!=1)
  y = rbinom(n,1,mu)
  df = data.frame(y = y, a = as.factor(a), x)
  rm("a","x","y")

  out = optSort(df, '~/temporaryoutputs/', optimizers = c("unconstrained", "constrained", "approximate"))
  expect_less_than(out$unconstrained$results$Estimate[2],.4)
  expect_less_than(out$constrained$results$Estimate[2],.7)
  expect_less_than(out$approximate$results$Estimate[2],.7)
})

test_that('we get about the result we expect', {
  # making astar less noisy
  # mean is .3 if at optimal a, .5 otherwise
  # under constraint, optimal is ~.34
  require(SuperLearner)
  n = 10000
  x = matrix(rnorm(2*n), ncol = 2)
  p1 = 0.5*(x[,1] > 1) + 0.2; p2 = 0.5*(x[,1] < -1) + 0.2; p3 = 0.6*( (x[,1]> -1)*(x[,1] < 1)) + .1
  a = apply(cbind(p1,p2,p3), 1, function(k) sample(c(1:3), 1, replace = F, prob = k))

  p1s = 0.8*(x[,1] > 0) + 0.1; p2s = 0.8*(x[,1] < -1/2) + 0.1; p3s = 0.8*( (x[,1]> -1/2)*(x[,1] < 0))
  astar = apply(cbind(p1s,p2s,p3s), 1, function(k) sample(c(1:3), 1, replace = F, prob = k))

  mu = runif(n,.25,.75) - 0.2*(a==astar)
  y = rbinom(n,1,mu)
  df = data.frame(y = y, a = as.factor(a), x)
  truth = .5 - .2*(a==astar)
  rm("a","x","y")

  out = optSort(df, '~/temporaryoutputs/', optimizers = c("unconstrained", "constrained", "approximate"))
  expect_less_than(out$unconstrained$results$Estimate[2],.35)
  expect_less_than(out$constrained$results$Estimate[2],.4)
  expect_less_than(out$approximate$results$Estimate[2],.4)
})

test_that('we get about the result we expect', {
  # mean is .3 if at optimal a, .5 otherwise
  # under constraint, optimal is ~.34

  # this can't estimate low enough muhats because astar is too noisy?
  require(SuperLearner)
  n = 10000
  x = matrix(rnorm(2*n), ncol = 2)
  p1 = 0.5*(x[,1] > 1) + 0.2; p2 = 0.5*(x[,1] < -1) + 0.2; p3 = 0.6*( (x[,1]> -1)*(x[,1] < 1)) + .1
  a = apply(cbind(p1,p2,p3), 1, function(k) sample(c(1:3), 1, replace = F, prob = k))

  p1s = 0.5*(x[,1] > 0) + 0.2; p2s = 0.5*(x[,1] < -1/2) + 0.2; p3s = 0.6*( (x[,1]> -1/2)*(x[,1] < 0)) + .1
  astar = apply(cbind(p1s,p2s,p3s), 1, function(k) sample(c(1:3), 1, replace = F, prob = k))

  mu = runif(n,.25,.75) - 0.2*(a==astar)
  y = rbinom(n,1,mu)
  df = data.frame(y = y, a = as.factor(a), x)
  truth = .5 - .2*(a==astar)
  rm("a","x","y")

  out = optSort(df, '~/temporaryoutputs/', optimizers = c("unconstrained", "constrained", "approximate"))

  # get truth at estimated fhat (psi_tilde)
  true_mumat <- sapply(sort(unique(df$a)), function(k) .5-.2*(astar==k))
  true_pimat <- cbind(p1,p2,p3)
  fhat.mat <- sapply(sort(unique(df$a)), function(a) as.numeric(out$unconstrained$assig.vec == a))
  muhat = diag(true_mumat %*% t(fhat.mat))
  pihat = diag(true_pimat %*% t(fhat.mat))
  if_tilde = (as.numeric(df$a == out$unconstrained$assig.vec)/pihat)*(df$y - muhat) + muhat
  psi_tilde = mean(if_tilde)

  expect_less_than((out$unconstrained$results$Estimate[2] - psi_tilde)^2,1)
  expect_less_than(out$constrained$results$Estimate[2],mean(df$y))
  expect_less_than(out$approximate$results$Estimate[2],mean(df$y))
})
