#' Single Shift Estimator Across Delta Range
#'
#' @description Estimates the effects of shifting up
#' across a range of values. Includes the multiplier bootstrap
#' for calculations of uniform confidence bands
#'
#' @param df full dataframe with the outcome labeled "y" and the treatment labeled "a"
#' @param output_root file path to put temporary output files
#' @param fudge degree to increase constrained over observed levels; default is 5%
#' @param nsplits number of data splits; default is 2
#' @param epsilon cutoff for extreme pi values; default is 1e-3
#' @param sections list of subsets of a values; defaults to empty. If it is not empty,
#' units with observed a values in a subset may only be sorted into the same subset.
#' @param opt_out_name name of csv file to be sent to matlab for optimization
#' @param optimizers which optimizers will be run; default is c("unconstrained", "constrained", "approximate")
#' @param sl.lib default: c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.earth","SL.xgboost","SL.glm")
#' @param sl.lib.pi default: c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.lda","SL.xgboost","SL.polymars")
#'
#' @return a list including an estimate of the effect for each delta value and its standard deviation,


optSort <- function(df, output_root = "~", fudge = .05, nsplits = 2,epsilon = 1e-3, sections = list(),
                    optimizers = c("unconstrained", "constrained", "approximate"), script_name = "prison_assignment_sl_nm.m",
                    sl.lib = c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.earth","SL.xgboost","SL.glm"),
                    sl.lib.pi = c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.lda","SL.xgboost","SL.polymars")){
  n = dim(df)[1]; p = length(unique(df$a))
  s = sample(rep(1:nsplits,ceiling(n/nsplits))[1:n])

  etas = nuisance_est(df, s, output_root, fudge, nsplits ,epsilon, sections, sl.lib, sl.lib.pi)

  muhat.mat = etas$muhat.mat
  pihat.mat = etas$pihat.mat
  avals = unique(df$a)

  if("unconstrained" %in% tolower(optimizers)){ unconstrained.output = unconstrained_opt(df,muhat.mat,avals) }
  #if("constrained" %in% tolower(optimizers)){ constrained.output = constrained_opt(df, muhat.mat, pihat.mat, s, output_root, nsplits, fudge, avals)}
  if("constrained" %in% tolower(optimizers)){ constrained.output = constrained_opt(df, muhat.mat, pihat.mat, output_root, fudge, avals)}
  if("approximate" %in% tolower(optimizers)){ approx.output = approximate_opt()}
}
