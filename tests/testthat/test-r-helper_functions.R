context("test-r-helper_functions")
require(testthat)

test_that("expit no greater than 1", {
  expect_true(expit(100)<=1)
})


# test_that('pi mat works', {
#   require(SuperLearner)
#   x = data.frame(X1 = rbeta(1000,1,2), X2 = rbeta(1000,1,2))
#   avals = c(1:4)
#   as = as.factor(sample(avals, 1000, replace = T))
#   xtest = data.frame(X1 = rbeta(1000,1,2), X2 = rbeta(1000,1,2))
#   sl.lib.pi = c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.lda","SL.xgboost","SL.polymars")
#   mat = make_pi_matrix(avals = c(1:4), train_a = as, train_x = x, test_x = xtest, epsilon = 1e-3, sl.lib.pi=sl.lib.pi)
#   expect_true(all(mat<1) & all(mat>0))
# })

# test_that('mu mat works', {
#   require(SuperLearner)
#   x = data.frame(X1 = rbeta(1000,1,2), X2 = rbeta(1000,1,2))
#   avals = c(1:3)
#   as = sample(avals, 1000, replace = T)
#   a = as.factor(as)
#   y = rbinom(1000,1,as/4)
#   xtest = data.frame(X1 = rbeta(1000,1,2), X2 = rbeta(1000,1,2))
#   sl.lib = c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.lda","SL.xgboost","SL.polymars")
#   make_mu_matrix(train_y=y, train_a = a, train_x=x, test_x=xtest, avals=avals, sl.lib=sl.lib)
# })

# test_that('nuis est works', {
#   require(SuperLearner)
#   x = data.frame(X1 = rbeta(1000,1,2), X2 = rbeta(1000,1,2))
#   as = sample(c(1:3), 1000, replace = T)
#   y = rbinom(1000,1,as/4)
#   df = data.frame(y = y, a = as.factor(as), x)
#   sl.lib = c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.lda","SL.xgboost","SL.polymars")
#   nuisance_est(df, output_root = "~/temporaryoutputs/", fudge = .05, nsplits = 2,epsilon = 1e-3)
# })

# test_that('nuis est works w subsets', {
#   require(SuperLearner)
#   x = data.frame(X1 = rbeta(1000,1,2), X2 = rbeta(1000,1,2))
#   as = sample(c(1:4), 1000, replace = T)
#   y = rbinom(1000,1,as/5)
#   s = sample(rep(1:2,ceiling(1000/2))[1:1000])
#   df = data.frame(y = y, a = as.factor(as), x)
#   sl.lib = c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.lda","SL.xgboost","SL.polymars")
#   out = nuisance_est(df, s,output_root = "~/temporaryoutputs/", fudge = .05, nsplits = 2,epsilon = 1e-3, sections = list(c(1,2),c(3,4)),
#                      sl.lib,sl.lib)
# })

test_that('constrained opt works', {
  require(SuperLearner)
  x = data.frame(X1 = rbeta(1000,1,2), X2 = rbeta(1000,1,2))
  as = sample(c(1:4), 1000, replace = T)
  y = rbinom(1000,1,as/5)
  s = sample(rep(1:2,ceiling(1000/2))[1:1000])
  df = data.frame(y = y, a = as.factor(as), x)
  sl.lib = c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.lda","SL.xgboost","SL.polymars")
  out = optSort(df, '~/temporaryoutputs/', optimizers = "constrained")
})
