#' Causal estimat of optimal treatment across a constrained or unconstrained objective
#'
#' @description Estimates the effect of an optimal treatment
#'
#' @param df full dataframe with the outcome labeled "y" and the treatment (factor) labeled "a"
#' @param output_root file path to put temporary output files
#' @param fudge degree to increase constrained over observed levels; default is 5%
#' @param nsplits number of data splits; default is 3 (one to train assignment vector model), must be >=3
#' @param epsilon cutoff for extreme pi values; default is 1e-3
#' @param sections list of subsets of a values; defaults to empty. If it is not empty,
#' units with observed a values in a subset may only be sorted into the same subset.
#' @param opt_out_name name of csv file to be sent to matlab for optimization
#' @param optimizers which optimizers will be run; default is c("unconstrained", "constrained", "approximate")
#' @param sl.lib default: c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.earth","SL.xgboost","SL.glm")
#' @param sl.lib.pi default: c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.lda","SL.xgboost","SL.polymars")
#'
#' @return a list including an estimate of the effect for each delta value and its standard deviation,


optSort <- function(df, output_root = "~", fudge = .05, nsplits = 3,epsilon = 1e-3, sections = list(),
                    optimizers = c("unconstrained", "constrained", "approximate"), script_name = "prison_assignment_sl_nm.m",
                    sl.lib = c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.lda","SL.earth","SL.xgboost","SL.polymars","SL.glm"),
                    sl.lib.pi = c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.lda","SL.earth","SL.xgboost","SL.polymars","SL.glm")){
  n = dim(df)[1]; p = length(unique(df$a))
  s = sample(rep(1:nsplits,ceiling(n/nsplits))[1:n])
  avals = sort(unique(df$a))

  if(!("a"%in% names(df)) | !("y"%in% names(df)) ){stop("Treatment needs to be named 'a'; outcome needs to be named 'y'")}
  if(!(is.factor(df$a))){df$a = as.factor(df$a); warning("Converting a to factor")}
  if(fudge>.5){warning("Fudge factor over 0.5, likely too high")}
  if(nsplits<3){stop("must have at least 3 folds to train assignment vector")}

  etas = nuisance_est(df, s, output_root, fudge, nsplits ,epsilon, sections, sl.lib, sl.lib.pi,avals)

  write.csv(etas$muhat.mat, paste(output_root,"saved_muhat_mat.csv", sep = ""))
  write.csv(etas$pihat.mat, paste(output_root,"saved_pihat_mat.csv", sep = ""))

  unconstrained.output <- constrained.output <- approx.output <- list()
  if("unconstrained" %in% tolower(optimizers)){ unconstrained.output = unconstrained_opt(df,avals,output_root) }
  if("constrained" %in% tolower(optimizers)){ constrained.output = constrained_opt(df,output_root, fudge, avals)}
  if("approximate" %in% tolower(optimizers)){ approx.output = approximate_opt(df,output_root,fudge,sections,sl.lib,sl.lib.pi,s,nsplits)}

  return(list(unconstrained = unconstrained.output, constrained = constrained.output, approximate = approx.output, etas = etas))
}
