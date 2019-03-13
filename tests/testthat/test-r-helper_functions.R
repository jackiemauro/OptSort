context("test-r-helper_functions")
require(testthat)

test_that("expit no greater than 1", {
  expect_true(expit(100)<=1)
})


test_that('constrained opt works', {
  require(SuperLearner)
  x = data.frame(X1 = rbeta(1000,1,2), X2 = rbeta(1000,1,2))
  as = sample(c(1:4), 1000, replace = T)
  y = rbinom(1000,1,as/5)
  df = data.frame(y = y, a = as.factor(as), x)
  sl.lib = c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.lda","SL.xgboost","SL.polymars")
  out = optSort(df, '~/temporaryoutputs/', optimizers = c("unconstrained","constrained"), sections = list(c(1,2),c(3,4)))
})
