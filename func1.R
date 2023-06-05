#####simulation setting####################

########################################################
####################setting 1###########################
########################################################
###生成真实数据
# install.packages("bindata")
library(bindata)
library(survRM2)
library(MASS)
library(randomForest)
library(survival)
library(rpart)
library(e1071)
library(neuralnet)
library(WeightIt)
library(nnet)
library(Rdonlp2)


###########################################################
####################平均生存时间###########################
###########################################################
true_race <- function(
  n.sample = n.sample,
  ua=ua,
  ub=ub,
  time = time,
  seed = seed){
  set.seed(seed)
  
  n.sample <- n.sample 

  mean <- c(0, 0)
  sigma1 <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)
  X1.5 <- mvrnorm(n.sample, mean, sigma1)
  X1 <- X1.5[,1]
  X5 <- X1.5[,2]

  mean <- c(0, 0)
  sigma2 <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)
  X2.6 <- mvrnorm(n.sample, mean, sigma2)
  X2 <- X2.6[,1]
  X6 <- X2.6[,2]

  mean <- c(0, 0)
  sigma3 <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)
  X3.8 <- mvrnorm(n.sample, mean, sigma3)
  X3 <- X3.8[,1]
  X8 <- X3.8[,2]

  mean <- c(0, 0)
  sigma4 <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)
  X4.9 <- mvrnorm(n.sample, mean, sigma4)
  X4 <- X4.9[,1]
  X9 <- X4.9[,2]
  
  
  X7 <- rnorm(n.sample,mean = 0,sd = 1)
  X10 <- rnorm(n.sample,mean = 0,sd = 1)
  
  
  a0 <- log(1.25)
  a1 <- log(1.25)
  a2 <- log(1.5)
  a3 <- log(1.75)
  a4 <- log(2)
  L1 <- -a0 + a1*X1 + a2*X2 + a3*X3 + a4*X4 + a1*X8 + a2*X9 + a3*X10
  L0 <-  a1*X1 + a2*X2 + a3*X3 + a4*X4 + a1*X8 + a2*X9 + a3*X10
  
  u <- runif(n.sample, ua,ub)
  
  Y1 <- (-log(u)/(0.0001*exp(L1)))^(1/3)  
  Y0 <- (-log(u)/(0.0001*exp(L0)))^(1/3)
  
  Y1 <- ifelse( Y1 >= time,time, Y1) 
  Y0 <- ifelse( Y0 >= time,time, Y0) 
  
  A <- mean(Y1) - mean(Y0)
  return(A)
}


#####################################################################
####################IPW rasce########################################
#####################################################################
###single method#####
sample_race <- function(
  n.sample = n.sample,
  ua=ua,ub=ub,
  model_type = model_type,
  time = time,  
  seed = seed,
  censor = censor){
  set.seed(seed)
  
  n.sample <- n.sample 
  
  
  mean <- c(0, 0)
  sigma1 <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)
  X1.5 <- mvrnorm(n.sample, mean, sigma1)
  X1 <- X1.5[,1]
  X5 <- X1.5[,2]
  
  mean <- c(0, 0)
  sigma2 <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)
  X2.6 <- mvrnorm(n.sample, mean, sigma2)
  X2 <- X2.6[,1]
  X6 <- X2.6[,2]
  

  mean <- c(0, 0)
  sigma3 <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)
  X3.8 <- mvrnorm(n.sample, mean, sigma3)
  X3 <- X3.8[,1]
  X8 <- X3.8[,2]
  
 
  mean <- c(0, 0)
  sigma4 <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)
  X4.9 <- mvrnorm(n.sample, mean, sigma4)
  X4 <- X4.9[,1]
  X9 <- X4.9[,2]
  
 
  X7 <- rnorm(n.sample,mean = 0,sd = 1)
  X10 <- rnorm(n.sample,mean = 0,sd = 1)
  
  b1 <- 0.8
  b2 <- -0.25
  b3 <- 0.6
  b4 <- -0.4
  b5 <- -0.8
  b6 <- -0.5
  b7 <- 0.7
  
  #true ps
  X.sum <- (X1 + X2 + X3 + X4 + X5 + X6 + X7)/7 
  Z <- rbinom( n.sample, 1, exp( X.sum )/( 1 + exp(X.sum)));
  ####simulate surivival data##########################################################
  u <- runif(n.sample,ua, ub)
  L <- rep(NA, n.sample)
  
  a0 <- log(1.25)
  a1 <- log(1.25)
  a2 <- log(1.5)
  a3 <- log(1.75)
  a4 <- log(2)
  
  L[Z==1] <- -a0 + a1*X1[Z==1] + a2*X2[Z==1] + a3*X3[Z==1] + a4*X4[Z==1] + a1*X8[Z==1] + a2*X9[Z==1] + a3*X10[Z==1]
  L[Z==0] <- a1*X1[Z==0] + a2*X2[Z==0] + a3*X3[Z==0] + a4*X4[Z==0] + a1*X8[Z==0] + a2*X9[Z==0] + a3*X10[Z==0]
  
  Y <- rep(NA, n.sample);
  Y[Z==1] <- ((-log(u[Z==1]))/(0.0001*exp(L[Z==1])))^(1/3) ; 
  Y[Z==0] <- ((-log(u[Z==0]))/(0.0001*exp(L[Z==0])))^(1/3) ;
  
  
  censor.time <- runif(n.sample, 0, censor);
  delta <- as.numeric( Y <= censor.time ); 
  Y[delta == 0] <- censor.time[delta == 0]
  
  data <- as.data.frame(cbind(X1, X2, X3, X4, X5,X6, 
                              X7,X8,X9,X10, Z, Y, delta))
  

  X4PS <- c(paste0("X",1:10))
  form <- as.formula( paste( "as.factor(Z) ~ ", paste( X4PS, collapse = " + " ), sep="" ) )

  out.ps <- ps.model(data = data, form = form, model_type = model_type)
  ps <- out.ps$ps.hat
  data$weight <- data$Z/ps + (1-data$Z)/(1-ps) 
  

  A <- akm_rmst(time = data$Y, 
                status = data$delta, 
                group = as.factor(data$Z), 
                weight=data$weight, 
                tau= time) 
  return(A)
}


###############################################################################
#####model averaging race######################################################
###############################################################################
avg_race <-function(
  n.sample = n.sample,
  ua=ua,ub=ub,
  time = time,
  seed = seed,
  censor = censor,
  method_smd = method_smd,
  model_type1 = model_type1,
  model_type2 = model_type2){
  set.seed(seed)
  
  n.sample <- n.sample 
  

  mean <- c(0, 0)
  sigma1 <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)
  X1.5 <- mvrnorm(n.sample, mean, sigma1)
  X1 <- X1.5[,1]
  X5 <- X1.5[,2]
  
  mean <- c(0, 0)
  sigma2 <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)
  X2.6 <- mvrnorm(n.sample, mean, sigma2)
  X2 <- X2.6[,1]
  X6 <- X2.6[,2]
  

  mean <- c(0, 0)
  sigma3 <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)
  X3.8 <- mvrnorm(n.sample, mean, sigma3)
  X3 <- X3.8[,1]
  X8 <- X3.8[,2]
  
  mean <- c(0, 0)
  sigma4 <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)
  X4.9 <- mvrnorm(n.sample, mean, sigma4)
  X4 <- X4.9[,1]
  X9 <- X4.9[,2]
  

  X7 <- rnorm(n.sample,mean = 0,sd = 1)
  X10 <- rnorm(n.sample,mean = 0,sd = 1)
  

  b1 <- 0.8
  b2 <- -0.25
  b3 <- 0.6
  b4 <- -0.4
  b5 <- -0.8
  b6 <- -0.5
  b7 <- 0.7
  

  X.sum <- (X1 + X2 + X3 + X4 + X5 + X6 + X7)/7 
  Z <- rbinom( n.sample, 1, exp( X.sum )/( 1 + exp(X.sum))); 
  
  u <- runif(n.sample,ua, ub)
  L <- rep(NA, n.sample)
  
  a0 <- log(1.25)
  a1 <- log(1.25)
  a2 <- log(1.5)
  a3 <- log(1.75)
  a4 <- log(2)
  
  L[Z==1] <- -a0 + a1*X1[Z==1] + a2*X2[Z==1] + a3*X3[Z==1] + a4*X4[Z==1] + a1*X8[Z==1] + a2*X9[Z==1] + a3*X10[Z==1]
  L[Z==0] <- a1*X1[Z==0] + a2*X2[Z==0] + a3*X3[Z==0] + a4*X4[Z==0] + a1*X8[Z==0] + a2*X9[Z==0] + a3*X10[Z==0]
  
  Y <- rep(NA, n.sample);
  Y[Z==1] <- ((-log(u[Z==1]))/(0.0001*exp(L[Z==1])))^(1/3) ; 
  Y[Z==0] <- ((-log(u[Z==0]))/(0.0001*exp(L[Z==0])))^(1/3) ;
  

  censor.time <- runif(n.sample, 0, censor);
  delta <- as.numeric( Y <= censor.time ); 
  Y[delta == 0] <- censor.time[delta == 0]
  
  data <- as.data.frame(cbind(X1, X2, X3, X4, X5,X6, 
                              X7,X8,X9,X10, Z, Y, delta))
  
 
  X4PS <- c(paste0("X",1:10))
  form <- as.formula( paste( "as.factor(Z) ~ ", paste( X4PS, collapse = " + " ), sep="" ) )
  
  out.ps1 <- ps.model(data = data, form = form,model_type = model_type1)
  ps1 <- out.ps1$ps.hat
  
  out.ps2 <- ps.model(data = data, form = form,model_type = model_type2)
  ps2 <- out.ps2$ps.hat
  
  
  save(ps1, file = "ps.1.Rdata")
  save(ps2, file = "ps.2.Rdata")
  save(data, file ="data.Rdata")
  if(method_smd == "mean_smd"){
    p = c(0,0)
    par.l = c(0,0); par.u = c(1,1); 
    A <- matrix(c(1,1),1, byrow = T)
    lin.l = 1;lin.u = 1
    
    ret = donlp2(p,f1,par.u = par.u, par.l = par.l,
                 A,lin.l = lin.l,lin.u =lin.u)
    
    lamda <- ret$par
    print(lamda)
    ps <- lamda[1] * ps1 + lamda[2] * ps2 
    data$weight <- data$Z/ps + (1-data$Z)/(1-ps)

    A <- akm_rmst(time = data$Y, 
                  status = data$delta, 
                  group = as.factor(data$Z), 
                  weight=data$weight, 
                  tau= time) #tau
    return(A)
    
  }else if(method_smd == "median_smd"){
    p = c(0,0)
    par.l = c(0,0); par.u = c(1,1); 
    A <- matrix(c(1,1),1, byrow = T)
    lin.l = 1;lin.u = 1
    
    ret = donlp2(p,f2,par.u = par.u, par.l = par.l,
                 A,lin.l = lin.l,lin.u =lin.u)
    
    lamda <- ret$par
    print(lamda)
    ps <- lamda[1] * ps1 + lamda[2] * ps2 
    
    data$weight <- data$Z/ps + (1-data$Z)/(1-ps) 
    
    A <- akm_rmst(time = data$Y, 
                  status = data$delta, 
                  group = as.factor(data$Z), 
                  weight=data$weight, 
                  tau= time) #tau
    return(A)
    
  }else if(method_smd == "max_smd") {
    p = c(0,0)
    par.l = c(0,0); par.u = c(1,1); 
    A <- matrix(c(1,1),1, byrow = T)
    lin.l = 1;lin.u = 1
    
    ret = donlp2(p,f3,par.u = par.u, par.l = par.l,
                 A,lin.l = lin.l,lin.u =lin.u)
    
    lamda <- ret$par
    print(lamda)
    ps <- lamda[1] * ps1 + lamda[2] * ps2 
    data$weight <- data$Z/ps + (1-data$Z)/(1-ps) 
    
    A <- akm_rmst(time = data$Y, 
                  status = data$delta, 
                  group = as.factor(data$Z), 
                  weight=data$weight, 
                  tau= time) #tau
    return(A)
  }else{
    p = c(0,0)
    par.l = c(0,0); par.u = c(1,1); 
    A <- matrix(c(1,1),1, byrow = T)
    lin.l = 1;lin.u = 1
    
    ret = donlp2(p,f4,par.u = par.u, par.l = par.l,
                 A,lin.l = lin.l,lin.u =lin.u)
    
    lamda <- ret$par
    print(lamda)
    ps <- lamda[1] * ps1 + lamda[2] * ps2 
    data$weight <- data$Z/ps + (1-data$Z)/(1-ps) 

    A <- akm_rmst(time = data$Y, 
                  status = data$delta, 
                  group = as.factor(data$Z), 
                  weight=data$weight, 
                  tau= time) #tau
    return(A)
  }
}


###################################
#########model averaging 7 models##
###################################
avg_race_m7 <- function(
  n.sample = n.sample,
  ua=ua,ub=ub,
  time = time,
  seed = seed,
  censor = censor,
  method_smd = method_smd,
  model_type1 = model_type1,
  model_type2 = model_type2,
  model_type3 = model_type3,
  model_type4 = model_type4,
  model_type5 = model_type5,
  model_type6 = model_type6,
  model_type7 = model_type7){
  set.seed(seed)
  
  n.sample <- n.sample 
  
  mean <- c(0, 0)
  sigma1 <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)
  X1.5 <- mvrnorm(n.sample, mean, sigma1)
  X1 <- X1.5[,1]
  X5 <- X1.5[,2]
  
  mean <- c(0, 0)
  sigma2 <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)
  X2.6 <- mvrnorm(n.sample, mean, sigma2)
  X2 <- X2.6[,1]
  X6 <- X2.6[,2]

  mean <- c(0, 0)
  sigma3 <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)
  X3.8 <- mvrnorm(n.sample, mean, sigma3)
  X3 <- X3.8[,1]
  X8 <- X3.8[,2]
  
  mean <- c(0, 0)
  sigma4 <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)
  X4.9 <- mvrnorm(n.sample, mean, sigma4)
  X4 <- X4.9[,1]
  X9 <- X4.9[,2]
  
  X7 <- rnorm(n.sample,mean = 0,sd = 1)
  X10 <- rnorm(n.sample,mean = 0,sd = 1)
  
  b1 <- 0.8
  b2 <- -0.25
  b3 <- 0.6
  b4 <- -0.4
  b5 <- -0.8
  b6 <- -0.5
  b7 <- 0.7
  
  X.sum <- (X1 + X2 + X3 + X4 + X5 + X6 + X7)/7  
  Z <- rbinom( n.sample, 1, exp( X.sum )/( 1 + exp(X.sum))); 
  
  u <- runif(n.sample,ua, ub)
  L <- rep(NA, n.sample)
  
  a0 <- log(1.25)
  a1 <- log(1.25)
  a2 <- log(1.5)
  a3 <- log(1.75)
  a4 <- log(2)
  
  L[Z==1] <- -a0 + a1*X1[Z==1] + a2*X2[Z==1] + a3*X3[Z==1] + a4*X4[Z==1] + a1*X8[Z==1] + a2*X9[Z==1] + a3*X10[Z==1]
  L[Z==0] <- a1*X1[Z==0] + a2*X2[Z==0] + a3*X3[Z==0] + a4*X4[Z==0] + a1*X8[Z==0] + a2*X9[Z==0] + a3*X10[Z==0]
  
  Y <- rep(NA, n.sample);
  Y[Z==1] <- ((-log(u[Z==1]))/(0.0001*exp(L[Z==1])))^(1/3) ; 
  Y[Z==0] <- ((-log(u[Z==0]))/(0.0001*exp(L[Z==0])))^(1/3) ;
  
  censor.time <- runif(n.sample, 0, censor);
  delta <- as.numeric( Y <= censor.time ); 
  Y[delta == 0] <- censor.time[delta == 0]
  
  data <- as.data.frame(cbind(X1, X2, X3, X4, X5,X6, 
                              X7,X8,X9,X10, Z, Y, delta))

  X4PS <- c(paste0("X",1:10))
  form <- as.formula( paste( "as.factor(Z) ~ ", paste( X4PS, collapse = " + " ), sep="" ) )
  
  out.ps1 <- ps.model(data = data, form = form,model_type = model_type1)
  ps1 <- out.ps1$ps.hat
  
  out.ps2 <- ps.model(data = data, form = form,model_type = model_type2)
  ps2 <- out.ps2$ps.hat
  
  out.ps3 <- ps.model(data = data, form = form,model_type = model_type3)
  ps3 <- out.ps3$ps.hat
  
  out.ps4 <- ps.model(data = data, form = form,model_type = model_type4)
  ps4 <- out.ps4$ps.hat
  
  out.ps5 <- ps.model(data = data, form = form,model_type = model_type5)
  ps5 <- out.ps5$ps.hat
  
  out.ps6 <- ps.model(data = data, form = form,model_type = model_type6)
  ps6 <- out.ps6$ps.hat
  
  out.ps7 <- ps.model(data = data, form = form,model_type = model_type7)
  ps7 <- out.ps7$ps.hat
  
  save(ps1, file = "ps.1.Rdata")
  save(ps2, file = "ps.2.Rdata")
  save(ps3, file = "ps.3.Rdata")
  save(ps4, file = "ps.4.Rdata")
  save(ps5, file = "ps.5.Rdata")
  save(ps6, file = "ps.6.Rdata")
  save(ps7, file = "ps.7.Rdata")
  
  save(data, file ="data.Rdata")
  if(method_smd == "mean_smd"){
    p = c(rep(0,7))
    par.l = c(rep(0,7)); par.u = c(rep(1,7)); 
    A <- matrix(c(rep(1,7)),1, byrow = T)
    lin.l = 1;lin.u = 1
    
    ret = donlp2(p,s1,par.u = par.u, par.l = par.l,
                 A,lin.l = lin.l,lin.u =lin.u)
    
    lamda <- ret$par
    
      ps <- lamda[1] * ps1 + lamda[2]* ps2 + lamda[3] * ps3 +
        lamda[4] * ps4 + lamda[5]* ps5 + lamda[6] * ps6+ lamda[7] * ps7
      data$weight <- data$Z/ps + (1-data$Z)/(1-ps) #产生治疗组和对照组的权重
      
      ######计算线性平均生存时间#####
      A <- akm_rmst(time = data$Y, 
                    status = data$delta, 
                    group = as.factor(data$Z), 
                    weight=data$weight, 
                    tau= time) #tau
    
  }else if(method_smd == "median_smd"){
    p = c(rep(0,7))
    par.l = c(rep(0,7)); par.u = c(rep(1,7)); 
    A <- matrix(c(rep(1,7)),1, byrow = T)
    lin.l = 1;lin.u = 1
    
    ret = donlp2(p,s2,par.u = par.u, par.l = par.l,
                 A,lin.l = lin.l,lin.u =lin.u)
    
    lamda <- ret$par
      ps <- lamda[1] * ps1 + lamda[2]* ps2 + lamda[3] * ps3 +
        lamda[4] * ps4 + lamda[5]* ps5 + lamda[6] * ps6+ lamda[7] * ps7
      data$weight <- data$Z/ps + (1-data$Z)/(1-ps) #产生治疗组和对照组的权重
      
      ######计算线性平均生存时间#####
      A <- akm_rmst(time = data$Y, 
                    status = data$delta, 
                    group = as.factor(data$Z), 
                    weight=data$weight, 
                    tau= time) #tau
    
    
    
  }else if(method_smd == "max_smd") {
    p = c(rep(0,7))
    par.l = c(rep(0,7)); par.u = c(rep(1,7)); 
    A <- matrix(c(rep(1,7)),1, byrow = T)
    lin.l = 1;lin.u = 1
    
    ret = donlp2(p,s3,par.u = par.u, par.l = par.l,
                 A,lin.l = lin.l,lin.u =lin.u)
    
      lamda <- lamda/(sum(lamda))
      ps <- lamda[1] * ps1 + lamda[2]* ps2 + lamda[3] * ps3 +
        lamda[4] * ps4 + lamda[5]* ps5 + lamda[6] * ps6+ lamda[7] * ps7
      data$weight <- data$Z/ps + (1-data$Z)/(1-ps) #产生治疗组和对照组的权重
      
      ######计算线性平均生存时间#####
      A <- akm_rmst(time = data$Y, 
                    status = data$delta, 
                    group = as.factor(data$Z), 
                    weight=data$weight, 
                    tau= time) #tau
    
  }else{
    p = c(rep(0,7))
    par.l = c(rep(0,7)); par.u = c(rep(1,7)); 
    A <- matrix(c(rep(1,7)),1, byrow = T)
    lin.l = 1;lin.u = 1
    
    ret = donlp2(p,s4,par.u = par.u, par.l = par.l,
                 A,lin.l = lin.l,lin.u =lin.u)
    
    lamda <- ret$par
      ps <- lamda[1] * ps1 + lamda[2]* ps2 + lamda[3] * ps3 +
        lamda[4] * ps4 + lamda[5]* ps5 + lamda[6] * ps6+ lamda[7] * ps7
      data$weight <- data$Z/ps + (1-data$Z)/(1-ps) #产生治疗组和对照组的权重
      
      ######计算线性平均生存时间#####
      A <- akm_rmst(time = data$Y, 
                    status = data$delta, 
                    group = as.factor(data$Z), 
                    weight=data$weight, 
                    tau= time) #tau
  }
  return(A)
}
