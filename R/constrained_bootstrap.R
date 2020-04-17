# playing around with bootstrap
rm(list = ls())
library(ranger)

# required functions
simFunc <- function(n = 1000){
  x = matrix(rnorm(2*n), ncol = 2)
  p1 = 0.5*(x[,1] > 1) + 0.2; p2 = 0.5*(x[,1] < -1) + 0.2; p3 = 0.6*( (x[,1]> -1)*(x[,1] < 1)) + .1
  a = apply(cbind(p1,p2,p3), 1, function(k) sample(c(1:3), 1, replace = F, prob = k))

  p1s = 0.8*(x[,1] > 0) + 0.1; p2s = 0.8*(x[,1] < -1/2) + 0.1; p3s = 0.8*( (x[,1]> -1/2)*(x[,1] < 0))
  astar = apply(cbind(p1s,p2s,p3s), 1, function(k) sample(c(1:3), 1, replace = F, prob = k))

  mu = runif(n,.25,.75) - 0.2*(a==astar)
  y = rbinom(n,1,mu)
  df = data.frame(y = y, a = as.factor(a), x)
  list(df = df, astar = astar)
}
constraint_f<- function(i,j){(abs(dat[i,'X1']) <= abs(dat[j,'X1'])*2) & (abs(dat[i,'X1']) >= abs(dat[j,'X1'])/2)}
train_rf <- function(dat){return(ranger::ranger(y ~ ., data = dat))}
make_test_dat <- function(drows){
  start <- end <- drows
  end[1,'a'] <- start[2,'a']; end[2,'a'] <- start[1,'a']
  rbind(start,end)
}

# make data, constraint rule and data copy to start
simOut <- simFunc()
dat <- simOut$df


constrained_bootstrap <- function(rounds, dat){
  dat_copy <- dat

  # run over rounds
  switch_decision = rep(NA, rounds)
  pb <- txtProgressBar(min = 0, max = rounds, style = 3)
  for(round in 1:rounds){
    setTxtProgressBar(pb, round)

    # pick rows
    row_num <- sample(1:dim(dat)[1],1)
    constraint_set <- sapply(c(1:dim(dat)[1]), function(k) constraint_f(row_num,k))
    row_num2 <- sample(c(1:dim(dat)[1])[constraint_set],1)

    # make sure they have different initial treatments
    counter = 0
    while((dat_copy[row_num, 'a'] == dat_copy[row_num2, 'a']) & (counter < 5) ){
      counter = counter + 1
      row_num <- sample(1:dim(dat)[1],1)
      constraint_set <- sapply(c(1:dim(dat)[1]), function(k) constraint_f(row_num,k))
      row_num2 <- sample(c(1:dim(dat)[1])[constraint_set],1)
    }

    # train model on original data (excluding rows 1 and 2)
    train_dat <- dat[-c(row_num, row_num2),]
    model <- train_rf(train_dat)
    predictions <- predict(model, make_test_dat(dat_copy[c(row_num,row_num2),]))$pred
    switch_decision[round] <- sum(predictions[1:2]) > sum(predictions[3:4])
    dat_copy[row_num, 'a'] <- ifelse(switch_decision[round],  dat_copy[row_num2, 'a'], dat_copy[row_num, 'a'])
    dat_copy[row_num2, 'a'] <- ifelse(switch_decision[round],  dat_copy[row_num, 'a'], dat_copy[row_num2, 'a'])
  }
  return(list(final_data = dat_copy, switch_history = switch_decision))
}

temp = constrained_bootstrap(1000,dat)

table(dat$a, temp$final_data$a)
table(dat$a, simOut$astar)
table(temp$final_data$a, simOut$astar)

ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}
plot(ma(temp$switch_history, n = 100), type = 'l')

