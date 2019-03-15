#' Helper Functions
#'
#' @description Functions to estimate muhats and pihats
#'
#' @return matrices of predicted pi's or mu's

globalVariables(package = "SuperLearner")
globalVariables('root')

expit <- function(x) {exp(x)/(1 + exp(x))}

make_mu_matrix<-function(train_y, train_a, train_x, test_x, avals, sl.lib){
  requireNamespace("SuperLearner")
  muhat = sapply(avals, function(k) predict_mu(train_y, train_a, k, train_x, sl.lib, test_x))
  return(muhat)
}

predict_mu <- function(train_y, train_a, a, train_x, sl.lib, test_x){
  mu.model = SuperLearner(Y = train_y[train_a==a], X = train_x[train_a==a,], family = binomial(), SL.library = sl.lib, verbose = FALSE)
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
  pi_model = invisible(SuperLearner(Y = as.numeric(train_a == a), X = train_x, family = binomial(), SL.library = sl.lib.pi))
  return(c(predict.SuperLearner(object = pi_model, newdata = test_x, onlySL = T)$pred))
}

nuisance_est <- function(df, s, output_root, fudge, nsplits, epsilon, sections, sl.lib, sl.lib.pi){

  # splits into 3 samples (fixed):
  # sample 1: train mu and pi
  # sample 2: train mu2 for f
  # sample 3: predict muhat and pihat from sample1; predict f = argmin(muhat2) from sample 2

  n = dim(df)[1]; p = length(unique(df$a))
  muhat.mat <- pihat.mat <- muhat2.mat <- matrix(rep(NA, n*p), ncol = p); assig.vec <- rep(NA,n)

  for(vfold in 1:3){

    train = (s==vfold)
    train2 = sapply(s, function(k) ifelse(vfold<nsplits, k==vfold+1, k==1) )
    test = !(train | train2)

    if(length(sections)==0){
      train_y = df$y[train]; train_a = df$a[train];
      train_x = subset(df, select = -c(a, y))[train,]; test_x = subset(df, select = -c(a, y))[test,]
      avals = sort(unique(train_a))
      muhat.mat[test,] = make_mu_matrix(train_y, train_a, train_x, test_x, avals, sl.lib)
      pihat.mat[test,] = make_pi_matrix(avals, train_a, train_x, test_x, epsilon, sl.lib.pi)

      train_2y = df$y[train2]; train_2a = df$a[train2];train_2x = subset(df, select = -c(a, y))[train2,]
      muhat2.mat[test,] = make_mu_matrix(train_2y, train_2a, train_2x, test_x, avals, sl.lib)
      assig.vec[test] = apply(muhat2.mat[test,],1,which.min)
    }

    else{
      for(sec in sections){
        avals = sort(sec)
        train_cond = (train & (df$a %in% sec))
        train2_cond = (train2 & (df$a %in% sec))
        test_cond = (test & (df$a %in% sec))

        # estimate nuisance parameters within sections
        train_y = df$y[train_cond]; train_a = df$a[train_cond];
        train_x = subset(df, select = -c(a, y))[train_cond,]; test_x = subset(df, select = -c(a, y))[test_cond,]
        muhat.mat[test_cond, levels(df$a) %in% sec] = make_mu_matrix(train_y, train_a, train_x, test_x, avals, sl.lib)
        pihat.mat[test_cond, levels(df$a) %in% sec] = make_pi_matrix(avals, train_a, train_x, test_x, epsilon, sl.lib.pi)

        # replace NA's (wrong section) with large/small values
        muhat.mat[test_cond, !(levels(df$a) %in% sec)] = 1e10
        pihat.mat[test_cond, !(levels(df$a) %in% sec)] = 1e-3

        # estimate assignment vector
        train_2y = df$y[train2_cond]; train_2a = df$a[train2_cond];train_2x = subset(df, select = -c(a, y))[train2_cond,]
        muhat2.mat = make_mu_matrix(train_2y, train_2a, train_2x, test_x, avals, sl.lib)
        assig.vec[test_cond] = avals[apply(muhat2.mat,1,which.min)]
      }
    }
  }
  write.csv(muhat.mat, paste(output_root,"muhat_mat.csv", sep = ""))
  write.csv(muhat2.mat, paste(output_root,"muhat2_mat.csv", sep = ""))
  write.csv(pihat.mat, paste(output_root,"pihat_mat.csv", sep = ""))
  write.csv(assig.vec, paste(output_root,"assig_vec.csv", sep = ""))
  return(list(muhat.mat = muhat.mat, pihat.mat = pihat.mat, assig.vec = assig.vec, muhat2.mat = muhat2.mat))
}

unconstrained_opt <- function(df,avals,output_root){
  print("Estimating Unconstrained Optimal")
  # takes in muhat and pihat
  # finds unconstrained optimal assignment based on muhat
  # outputs plug in and if-based estimates

  muhat.mat = as.matrix(read.csv(paste(output_root,"muhat_mat.csv", sep = "")))[,-1]
  pihat.mat = as.matrix(read.csv(paste(output_root,"pihat_mat.csv", sep = "")))[,-1]
  fU = read.csv(paste(output_root,"assig_vec.csv", sep = ""))[,-1]

  pluginU = mean(apply(muhat.mat,1,min))
  fU.mat = sapply(c(1:length(avals)), function(a) as.numeric(fU == a))
  pihatU = diag(pihat.mat %*% t(fU.mat))
  muhatU = diag(muhat.mat %*% t(fU.mat))
  ifU = (as.numeric(df$a == avals[fU])/pihatU)*(df$y - muhatU) + muhatU
  psiU = mean(ifU)
  sdU = sd(ifU)/sqrt(length(muhatU))
  plugin.sdU = sd(muhatU)/sqrt(length(muhatU))

  res = data.frame(Estimate = c(pluginU, psiU), SD = c(plugin.sdU, sdU),
                   lower = c(pluginU - 1.96*plugin.sdU, psiU - 1.96*sdU),
                   upper = c(pluginU + 1.96*plugin.sdU, psiU + 1.96*sdU))
  rownames(res) = c("Plug in", "IF-based")
  print(res)

  return(list(results = res, assig.vec = fU, infl.func = ifU))
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

constrained_opt <- function(df, output_root, fudge, avals){
  print("Estimating Constrained Optimal")
  # not doing sample splitting
  # takes in muhat and pihat estimates
  # runs matlab code to find assignment that minimizes total muhat
  # calculates plugin and if-based estimates given optimal constrained assignment

  muhat.mat = as.matrix(read.csv(paste(output_root,"muhat_mat.csv", sep = "")))[,-1]
  pihat.mat = as.matrix(read.csv(paste(output_root,"pihat_mat.csv", sep = "")))[,-1]
  write.csv(round(c(table(df$a) * (1+fudge))), paste(output_root, 'constraint.csv', sep = ""))
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

  res = data.frame(Estimate = c(pluginC, psiC), SD = c(plugin.sdC, sdC),
                   lower = c(pluginC - 1.96*plugin.sdC, psiC - 1.96*sdC),
                   upper = c(pluginC + 1.96*plugin.sdC, psiC + 1.96*sdC))
  rownames(res) = c("Plug in", "IF-based")
  print(res)

  return(list(results = res, assig.vec = fC, infl.func = ifC))
}

approximate_opt <- function(df, output_root, fudge, sections,sl.lib = sl.lib, sl.lib.pi = sl.lib.pi){
  print("Estimating Approximate Constrained Optimal")
  # get appromximate constrained regression-based estimates
  n = dim(df)[1]; p = length(unique(df$a))
  fhat <- rep(NA,n)
  s = sample(rep(1:2,ceiling(n/2))[1:n])

  muhat.mat <- pihat.mat <- matrix(rep(NA, n*p), ncol = p)
  for(vfold in 1:2){
    # step 1: train and predict mu model on training data for each a
    train = (s==vfold); test = (s!=vfold)

    if(length(sections)==0){
      train_y = df$y[train]; train_a = df$a[train]; avals = sort(unique(train_a))
      train_x = subset(df, select = -c(a, y))[train,]; test_x = subset(df, select = -c(a, y))[test,]
      mu.model = SuperLearner(Y = train_y, X = data.frame(a = as.numeric(train_a), train_x), family = binomial(), SL.library = sl.lib)
      muhat.mat.train = sapply(avals, function(k) predict_mu(a=k, avals = avals, mu.model = mu.model, test_x = train_x)  )
      muhat.mat[test,] = sapply(avals, function(k) predict_mu(a=k, avals = avals, mu.model = mu.model, test_x = test_x)  )
      pihat.mat[test,] = make_pi_matrix(avals, train_a, train_x, test_x, epsilon=1e-3, sl.lib.pi)
    }
    else{
      muhat.mat.train <- matrix(rep(NA,sum(train)*p), ncol = p)
      for(sec in sections){
        avals = sort(sec); train_cond = (train & (df$a %in% sec));  test_cond = (test & (df$a %in% sec))

        # estimate nuisance parameters within sections
        train_y = df$y[train_cond]; train_a = df$a[train_cond]
        train_x = subset(df, select = -c(a, y))[train_cond,]
        test_x = subset(df, select = -c(a, y))[test_cond,]
        train_df = data.frame(a = as.numeric(train_a), train_x)
        mu.model = SuperLearner(Y = train_y, X = train_df, family = binomial(), SL.library = sl.lib)
        muhat.mat.train[df[train,]$a %in% sec,levels(df$a) %in% sec] = sapply(avals, function(k) predict_mu(a=k, avals = avals, mu.model = mu.model, test_x = train_x)  )
        muhat.mat[test_cond,levels(df$a) %in% sec] = sapply(avals, function(k) predict_mu(a=k, avals = avals, mu.model = mu.model, test_x = test_x)  )
        pihat.mat[test_cond,] = make_pi_matrix(avals, train_a, train_x, test_x, epsilon=1e-3, sl.lib.pi)

        # replace NA's (wrong section) with large values so optimization doesn't assign them
        muhat.mat.train[df[train,]$a %in% sec, !(levels(df$a) %in% sec)] = 1e10
        muhat.mat[test_cond, !(levels(df$a) %in% sec)] = 1e10
      }
    }
    # step 2: get constrained assignment based on step muhats
    write.csv(round(c(table(df[train,]$a) * (1+fudge))), paste(output_root, 'constraint.csv', sep = ""))
    write.csv(muhat.mat.train, paste(output_root,"muhat_mat.csv", sep = ""))
    run_matlab_script(paste(output_root,'constrained_optimizer_code.m',sep= ""))
    fC.mat = read.csv(paste(output_root,"constrained_optimized_assignment.csv", sep = ""), header = F)
    fC = apply(fC.mat,1,which.max)

    # step 3: train the model E(f|muhat) on the training data
    class.df = data.frame(a = as.factor(fC), muhat.mat.train)
    f.model = ranger::ranger(a~., data = class.df, write.forest = TRUE)

    # step 4: fhat on testing data
    fhat[test] = predict(f.model, data.frame(muhat.mat[test,]), type='response')$pre
  }
  res = get_res(df, fhat, pihat.mat, muhat.mat)

  return(list(results = res$res, assig.vec = fhat, infl.func = res$infl.func))
}

get_res <- function(df, fhat, pihat.mat, muhat.mat){
  avals = sort(unique(df$a))
  fhat.mat <- sapply(avals, function(a) as.numeric(fhat == a))
  pihat = diag(pihat.mat %*% t(fhat.mat))
  muhat = diag(muhat.mat %*% t(fhat.mat))
  plugin = mean(muhat)
  plugin.sd = sd(muhat)/sqrt(length(muhat))
  ifA = (as.numeric(df$a == avals[fhat])/pihat)*(df$y - muhat) + muhat
  psi = mean(ifA)
  sd = sd(ifA)/sqrt(length(muhat))

  res = data.frame(Estimate = c(plugin, psi), SD = c(plugin.sd, sd),
                   lower = c(plugin - 1.96*plugin.sd, psi - 1.96*sd),
                   upper = c(plugin + 1.96*plugin.sd, psi + 1.96*sd))
  rownames(res) = c("Plug in", "IF-based")
  print(res)

  return(list(res = res, infl.func = ifA))
}



