#' Helper Functions
#'
#' @description Functions to estimate muhats and pihats
#'
#' @return matrices of predicted pi's or mu's

globalVariables(package = "SuperLearner")
globalVariables('root')

expit <- function(x) {exp(x)/(1 + exp(x))}
output_constraints <- function(df, output_root, fudge = .05, suffix = ""){
  write.csv(round(table(df$a) * (1+fudge)), paste(output_root, 'contraint',suffix, '.csv', sep = ""))
  }

make_mu_matrix<-function(train_y, train_a, train_x, test_x, avals, sl.lib){
  requireNamespace("SuperLearner")
  train_df = data.frame(A = as.numeric(train_a), train_x)
  mu.model = SuperLearner(Y = train_y, X = train_df, family = binomial(), SL.library = sl.lib)
  muhat = sapply(avals, function(k) predict_mu(a=k, avals = avals, mu.model = mu.model, test_x = test_x)  )
  return(muhat)
}

predict_mu <- function(a, avals, mu.model, test_x){
  test_x$A = as.numeric(a)
  return(c(predict.SuperLearner(object = mu.model, newdata = test_x, onlySL = T)$pred))
}

make_pi_matrix <- function(avals, train_a, train_x, test_x, epsilon = 1e-3, sl.lib.pi){
  pihat = sapply(avals, function(k) predict_pi(a = k, train_a, train_x, test_x, sl.lib.pi = sl.lib.pi))

  # truncate extreme pi values
  print(paste(sum(pihat < epsilon | pihat > 1-epsilon), 'extreme pi values'))
  pihat[pihat < epsilon] = epsilon
  pihat[pihat > 1-epsilon] = 1-epsilon
  return(pihat)
}

predict_pi <- function(a, train_a, train_x, test_x, sl.lib.pi){
  #for a single a value, create a model of prob(A = a) and predict on new data#
  requireNamespace("SuperLearner")
  pi_model = SuperLearner(Y = as.numeric(train_a == a), X = train_x, family = binomial(), SL.library = sl.lib.pi)
  return(c(predict.SuperLearner(object = pi_model, newdata = test_x, onlySL = T)$pred))
}

nuisance_est <- function(df, s, output_root, fudge, nsplits ,epsilon, sections,sl.lib,sl.lib.pi){
  n = dim(df)[1]; p = length(unique(df$a))

  muhat.mat <- pihat.mat <- matrix(rep(NA, n*p), ncol = p)

  for(vfold in 1:nsplits){
    train = (s==vfold); test = (s!=vfold)
    output_constraints(df[test,], output_root, fudge, suffix = vfold)

    if(length(sections)==0){
      train_y = df$y[train]; train_a = df$a[train];
      train_x = subset(df, select = -c(a, y))[train,]; test_x = subset(df, select = -c(a, y))[test,]
      avals = sort(unique(train_a))
      muhat.mat[test,] = make_mu_matrix(train_y, train_a, train_x, test_x, avals, sl.lib)
      pihat.mat[test,] = make_pi_matrix(avals, train_a, train_x, test_x, epsilon, sl.lib.pi)
    }

    else{
      for(sec in sections){
        avals = sort(sec)
        train_cond = (train & (df$a %in% sec))
        test_cond = (test & (df$a %in% sec))

        train_y = df$y[train_cond]; train_a = df$a[train_cond];
        train_x = subset(df, select = -c(a, y))[train_cond,]; test_x = subset(df, select = -c(a, y))[test_cond,]
        muhat.mat[test_cond, levels(df$a) %in% sec] = make_mu_matrix(train_y, train_a, train_x, test_x, avals, sl.lib)
        pihat.mat[test_cond, levels(df$a) %in% sec] = make_pi_matrix(avals, train_a, train_x, test_x, epsilon, sl.lib.pi)

        # set to inf or epsilon if wrong sec level
        pihat.mat[test_cond, !(levels(df$a) %in% sec)] = epsilon
        muhat.mat[test_cond, !(levels(df$a) %in% sec)] = 1e10
      }
    }
  }
  return(list(muhat.mat = muhat.mat, pihat.mat = pihat.mat))
}

unconstrained_opt <- function(df,muhat.mat,avals){
  # takes in muhat and pihat
  # finds unconstrained optimal assignment based on muhat
  # outputs plug in and if-based estimates

  fU = apply(muhat.mat,1,which.min)
  pluginU = mean(apply(muhat.mat,1,min))
  fU.mat = sapply(c(1:length(avals)), function(a) as.numeric(fU == a))
  pihatU = diag(pihat.mat %*% t(fU.mat))
  muhatU = diag(muhat.mat %*% t(fU.mat))
  ifU = (as.numeric(df$a == avals[fU])/pihatU)*(df$y - muhatU) + muhatU
  psiU = mean(ifU)
  sdU = sd(ifU)/sqrt(length(muhatU))
  plugin.sdU = sd(muhatU)/sqrt(length(muhatU))

  return(list(assig.vec = fU, plugin = pluginU, plugin.sd = plugin.sdU, infl.func = ifU, psi = psiU, sd = sdU))
}

constrained_opt_split <- function(df, muhat.mat, pihat.mat, s, output_root, nsplits, fudge, avals){
  # not sure if this should split across samples...
  fC <- muhatC <- ifC <- rep(NA, dim(df)[1])
  for(vfold in 1:nsplits){
    train = (s==vfold); test = (s!=vfold)
    write.csv(round(c(table(df[test,]$a) * (1+fudge))), paste(output_root, 'constraint.csv', sep = ""))
    write.csv(muhat.mat[test,], paste(output_root,"muhat_mat.csv", sep = ""))
    run_matlab_script(paste(output_root,'constrained_optimizer_code.m',sep= ""))

    fC.mat = read.csv(paste(output_root,"constrained_optimized_assignment.csv", sep = ""), header = F)

    fC[test] = apply(fC.mat,1,which.max)
    pihatC = diag(pihat.mat[test,] %*% t(fC.mat))
    muhatC[test] = diag(muhat.mat[test,] %*% t(fC.mat))
  }
  ifC = (as.numeric(df$a == avals[fC])/pihatC)*(df$y - muhatC) + muhatC
  pluginC = mean(muhatC)
  plugin.sdC = sd(muhatC)/sqrt(length(muhatC))
  psiC = mean(ifC)
  sdC = sd(ifC)/sqrt(length(muhatC))

  return(list(assig.vec = fC, plugin = pluginC, plugin.sd = plugin.sdC, infl.func = ifC, psi = psiC, sd = sdC))
}

constrained_opt <- function(df, muhat.mat, pihat.mat, output_root, fudge, avals){
  # not doing sample splitting
  # takes in muhat and pihat estimates
  # runs matlab code to find assignment that minimizes total muhat
  # calculates plugin and if-based estimates given optimal constrained assignment

  write.csv(round(c(table(df$a) * (1+fudge))), paste(output_root, 'constraint.csv', sep = ""))
  write.csv(muhat.mat, paste(output_root,"muhat_mat.csv", sep = ""))
  run_matlab_script(paste(output_root,'constrained_optimizer_code.m',sep= ""))

  fC.mat = read.csv(paste(output_root,"constrained_optimized_assignment.csv", sep = ""), header = F)

  fC = apply(fC.mat,1,which.max)
  pihatC = diag(pihat.mat %*% t(fC.mat))
  muhatC = diag(muhat.mat %*% t(fC.mat))

  ifC = (as.numeric(df$a == avals[fC])/pihatC)*(df$y - muhatC) + muhatC
  pluginC = mean(muhatC)
  plugin.sdC = sd(muhatC)/sqrt(length(muhatC))
  psiC = mean(ifC)
  sdC = sd(ifC)/sqrt(length(muhatC))

  return(list(assig.vec = fC, plugin = pluginC, plugin.sd = plugin.sdC, infl.func = ifC, psi = psiC, sd = sdC))
}







