########################################################
######PS################################################
########################################################

library(randomForest)
library(Ecume)
ps.model <- function(model_type = model_type, data = data,form = form){
  if(model_type == "lg" ){
    fm <- glm(form, data = data, family = binomial(link = "logit"))
    ps.hat <- as.numeric(predict(fm, newdata = data, type = "response"))
    ps.hat <- pmin(pmax(0.000001, ps.hat), 0.999999)  
    return(list(ps.hat = ps.hat, fm = fm))  
    
  }else if(model_type == "rf" ){ 
    fm <- randomForest(form, data = data)
    ps.hat <- predict(fm, type = "prob")[,2]
    ps.hat <- pmin(pmax(0.000001, ps.hat), 0.999999)
    return(list(ps.hat = ps.hat, fm = fm)) 
    
  }else if(model_type == "svm_l"){
    fm <- svm(form, data= data,kernel = "linear", type = "C-classification", probability = T)
    ps.hat <- predict(fm, data,probability=TRUE)
    r.ps <- as.data.frame(attr(ps.hat, "probabilities"))
    num <- which(colnames(r.ps)== "1")
    ps.hat <-  as.vector(attr(ps.hat, "probabilities")[,num]) 
    ps.hat <- pmin(pmax(0.000001, ps.hat), 0.999999)
    return(list(ps.hat = ps.hat, fm = fm)) 
    
  }else if(model_type == "svm_s"){
    fm <- svm(form, data= data,kernel = "radial", type = "C-classification", probability = T)
    ps.hat <- predict(fm, data,probability=TRUE)
    ps.hat <- as.vector(attr(ps.hat, "probabilities")[,2])
    ps.hat <- pmin(pmax(0.000001, ps.hat), 0.999999)
    return(list(ps.hat = ps.hat, fm = fm)) 
    
  }else if(model_type == "cart"){ 
    fm <- rpart(form, data = data, method = "class")
    ps.hat <- predict(fm, type = "prob")
    ps.hat <- as.data.frame(ps.hat)
    num <- which(colnames(ps.hat) == "1")
    ps.hat <- ps.hat[,num]
    ps.hat <- pmin(pmax(0.000001, ps.hat), 0.999999)
    return(list(ps.hat = ps.hat, fm = fm)) 
    
  }else if(model_type == "pcart"){ 
    fm <- rpart(form, data = data, method = "class")
    fm <- prune(fm,cp = 0.02)
    ps.hat <- predict(fm, type = "prob")[,2]
    ps.hat <- pmin(pmax(0.000001, ps.hat), 0.999999)
    return(list(ps.hat = ps.hat, fm = fm)) 
    
  }else if(model_type == "nn"){ 

    fm <- nnet(form, data = data,size = 3)
    ps.hat <- fm$fitted.values
    ps.hat <- pmin(pmax(0.000001, ps.hat), 0.999999)
    return(list(ps.hat = ps.hat, fm = fm)) 
    
  }else if (model_type == "cbps"){ 
    fm <- weightit(form, data = data, method = "cbps", estimand = "ATE")
    ps.hat <- fm$ps
    ps.hat <- pmin(pmax(0.000001, ps.hat), 0.999999)
    return(list(ps.hat = ps.hat, fm = fm)) 
    
  }else { 
    fm <- weightit(form, data = data, method = "gbm", estimand = "ATE",
                   stop.method = "ks.mean",
                   shrinkage = 0.01 
    )
    ps.hat <- fm$ps
    ps.hat <- pmin(pmax(0.000001, ps.hat), 0.999999)
    return(list(ps.hat = ps.hat, fm = fm)) 
  }
}



#####################################################
########model averaging of two model#################
#####################################################
f1 <- function(lamda){
  # load("ps.1.Rdata",envir = )
  load("ps.1.Rdata",envir=globalenv())
  load("ps.2.Rdata",envir=globalenv())
  
  load("data.Rdata",envir=globalenv())
  
  ps1 <- get("ps1")
  ps2 <- get("ps2")
  
  data <- get("data")
  
  
  ps <- lamda[1] * ps1 + lamda[2]* ps2  
  data$ps <- ps
  treat <- data[data$Z == 1,]
  control <- data[data$Z == 0,]
  

  weight_treat <- 1/treat$ps
  weight_control <- 1/(1-control$ps)
  
  
  treat_1 <- treat[,c(paste0("X",1:10))]
  control_1 <- control[,c(paste0("X",1:10))]
  
  treat.smd.0 <- apply(treat_1,2,function(x){sum(x * weight_treat)/sum(weight_treat)})
  treat.sd.0 <- apply(treat_1,2,sd)
  treat.smd.1 <- treat.smd.0/treat.sd.0
  
  control.smd.0 <- apply(control_1,2,function(x){sum(x * weight_control)/sum(weight_control)})
  control.smd.1 <- control.smd.0/treat.sd.0

  smd_all <- abs(treat.smd.1 - control.smd.1)
  mean_smd <- mean(smd_all)
}

f2 <- function(lamda){
  load("ps.1.Rdata",envir=globalenv())
  load("ps.2.Rdata",envir=globalenv())
  
  load("data.Rdata",envir=globalenv())
  
  ps1 <- get("ps1")
  ps2 <- get("ps2")
  
  data <- get("data")
  
  ps <- lamda[1] * ps1 + lamda[2]* ps2  
  data$ps <- ps
  treat <- data[data$Z == 1,]
  control <- data[data$Z == 0,]
  
  
  weight_treat <- 1/treat$ps
  weight_control <- 1/(1-control$ps)
  
  
  treat_1 <- treat[,c(paste0("X",1:10))]
  control_1 <- control[,c(paste0("X",1:10))]
  
  treat.smd.0 <- apply(treat_1,2,function(x){sum(x * weight_treat)/sum(weight_treat)})
  treat.sd.0 <- apply(treat_1,2,sd)
  treat.smd.1 <- treat.smd.0/treat.sd.0
  
  control.smd.0 <- apply(control_1,2,function(x){sum(x * weight_control)/sum(weight_control)})
  control.smd.1 <- control.smd.0/treat.sd.0
  smd_all <- abs(treat.smd.1 - control.smd.1)
  
  median_smd <- median(smd_all)
}

f3 <- function(lamda){
  load("ps.1.Rdata",envir=globalenv())
  load("ps.2.Rdata",envir=globalenv())
  
  load("data.Rdata",envir=globalenv())
  
  ps1 <- get("ps1")
  ps2 <- get("ps2")
  
  data <- get("data")
  
  ps <- lamda[1] * ps1 + lamda[2]* ps2  
  data$ps <- ps
  treat <- data[data$Z == 1,]
  control <- data[data$Z == 0,]
  
  
  weight_treat <- 1/treat$ps
  weight_control <- 1/(1-control$ps)
  
  
  treat_1 <- treat[,c(paste0("X",1:10))]
  control_1 <- control[,c(paste0("X",1:10))]
  
  treat.smd.0 <- apply(treat_1,2,function(x){sum(x * weight_treat)/sum(weight_treat)})
  treat.sd.0 <- apply(treat_1,2,sd)
  treat.smd.1 <- treat.smd.0/treat.sd.0
  
  control.smd.0 <- apply(control_1,2,function(x){sum(x * weight_control)/sum(weight_control)})
  control.smd.1 <- control.smd.0/treat.sd.0
  smd_all <- abs(treat.smd.1 - control.smd.1)
  
  max_smd <- max(smd_all)
}

f4 <- function(lamda){
  load("ps.1.Rdata",envir=globalenv())
  load("ps.2.Rdata",envir=globalenv())
  
  load("data.Rdata",envir=globalenv())
  
  ps1 <- get("ps1")
  ps2 <- get("ps2")
  
  data <- get("data")
  
  ps <- lamda[1] * ps1 + lamda[2]* ps2  
  data$ps <- ps
  treat <- data[data$Z == 1,]
  control <- data[data$Z == 0,]
  
  
  weight_treat <- 1/treat$ps
  weight_control <- 1/(1-control$ps)
  
  D <- c()
  for(i in 1:10){
    x <- treat[,i,drop =T]
    y <- control[,i,drop =T]
    a <-ks_test(x,y,w_x = weight_treat, w_y = weight_control, thresh = 0)$statistic
    D <- c(D,a)
  }
  mean_D <- mean(D)
}

####################################################
######Seven models##################################
####################################################
s1 <- function(lamda){
  # load("ps.1.Rdata",envir = )
  load("ps.1.Rdata",envir=globalenv())
  load("ps.2.Rdata",envir=globalenv())
  load("ps.3.Rdata",envir=globalenv())
  load("ps.4.Rdata",envir=globalenv())
  load("ps.5.Rdata",envir=globalenv())
  load("ps.6.Rdata",envir=globalenv())
  load("ps.7.Rdata",envir=globalenv())
  load("data.Rdata",envir=globalenv())
  
  ps1 <- get("ps1")
  ps2 <- get("ps2")
  ps3 <- get("ps3")
  ps4 <- get("ps4")
  ps5 <- get("ps5")
  ps6 <- get("ps6")
  ps7 <- get("ps7")
  
  data <- get("data")
  
  
  ps <- lamda[1] * ps1 + lamda[2]* ps2 + lamda[3] * ps3 +
    lamda[4] * ps4 + lamda[5]* ps5 + lamda[6] * ps6+ lamda[7] * ps7
  data$ps <- ps
  treat <- data[data$Z == 1,]
  control <- data[data$Z == 0,]
  
  
  weight_treat <- 1/treat$ps
  weight_control <- 1/(1-control$ps)
  
  
  treat_1 <- treat[,c(paste0("X",1:10))]
  control_1 <- control[,c(paste0("X",1:10))]
  
  treat.smd.0 <- apply(treat_1,2,function(x){sum(x * weight_treat)/sum(weight_treat)})
  treat.sd.0 <- apply(treat_1,2,sd)
  treat.smd.1 <- treat.smd.0/treat.sd.0
  
  control.smd.0 <- apply(control_1,2,function(x){sum(x * weight_control)/sum(weight_control)})
  control.smd.1 <- control.smd.0/treat.sd.0
  smd_all <- abs(treat.smd.1 - control.smd.1)
  
  mean_smd <- mean(smd_all)
}

s2 <- function(lamda){
  # load("ps.1.Rdata",envir = )
  load("ps.1.Rdata",envir=globalenv())
  load("ps.2.Rdata",envir=globalenv())
  load("ps.3.Rdata",envir=globalenv())
  load("ps.4.Rdata",envir=globalenv())
  load("ps.5.Rdata",envir=globalenv())
  load("ps.6.Rdata",envir=globalenv())
  load("ps.7.Rdata",envir=globalenv())
  load("data.Rdata",envir=globalenv())
  
  ps1 <- get("ps1")
  ps2 <- get("ps2")
  ps3 <- get("ps3")
  ps4 <- get("ps4")
  ps5 <- get("ps5")
  ps6 <- get("ps6")
  ps7 <- get("ps7")
  
  data <- get("data")
  
  
  ps <- lamda[1] * ps1 + lamda[2]* ps2 + lamda[3] * ps3 +
    lamda[4] * ps4 + lamda[5]* ps5 + lamda[6] * ps6+ lamda[7] * ps7
  data$ps <- ps
  treat <- data[data$Z == 1,]
  control <- data[data$Z == 0,]
  
  
  weight_treat <- 1/treat$ps
  weight_control <- 1/(1-control$ps)
  
  
  treat_1 <- treat[,c(paste0("X",1:10))]
  control_1 <- control[,c(paste0("X",1:10))]
  
  treat.smd.0 <- apply(treat_1,2,function(x){sum(x * weight_treat)/sum(weight_treat)})
  treat.sd.0 <- apply(treat_1,2,sd)
  treat.smd.1 <- treat.smd.0/treat.sd.0
  
  control.smd.0 <- apply(control_1,2,function(x){sum(x * weight_control)/sum(weight_control)})
  control.smd.1 <- control.smd.0/treat.sd.0
  smd_all <- abs(treat.smd.1 - control.smd.1)
  
  median_smd <- median(smd_all)
}

s3 <- function(lamda){
  # load("ps.1.Rdata",envir = )
  load("ps.1.Rdata",envir=globalenv())
  load("ps.2.Rdata",envir=globalenv())
  load("ps.3.Rdata",envir=globalenv())
  load("ps.4.Rdata",envir=globalenv())
  load("ps.5.Rdata",envir=globalenv())
  load("ps.6.Rdata",envir=globalenv())
  load("ps.7.Rdata",envir=globalenv())
  load("data.Rdata",envir=globalenv())
  
  ps1 <- get("ps1")
  ps2 <- get("ps2")
  ps3 <- get("ps3")
  ps4 <- get("ps4")
  ps5 <- get("ps5")
  ps6 <- get("ps6")
  ps7 <- get("ps7")
  
  data <- get("data")
  
  
  ps <- lamda[1] * ps1 + lamda[2]* ps2 + lamda[3] * ps3 +
    lamda[4] * ps4 + lamda[5]* ps5 + lamda[6] * ps6+ lamda[7] * ps7
  data$ps <- ps
  treat <- data[data$Z == 1,]
  control <- data[data$Z == 0,]
  
  
  weight_treat <- 1/treat$ps
  weight_control <- 1/(1-control$ps)
  
  
  treat_1 <- treat[,c(paste0("X",1:10))]
  control_1 <- control[,c(paste0("X",1:10))]
  
  treat.smd.0 <- apply(treat_1,2,function(x){sum(x * weight_treat)/sum(weight_treat)})
  treat.sd.0 <- apply(treat_1,2,sd)
  treat.smd.1 <- treat.smd.0/treat.sd.0
  
  control.smd.0 <- apply(control_1,2,function(x){sum(x * weight_control)/sum(weight_control)})
  control.smd.1 <- control.smd.0/treat.sd.0
  smd_all <- abs(treat.smd.1 - control.smd.1)
  
  max_smd <- max(smd_all)
}

s4 <- function(lamda){
  load("ps.1.Rdata",envir=globalenv())
  load("ps.2.Rdata",envir=globalenv())
  load("ps.3.Rdata",envir=globalenv())
  load("ps.4.Rdata",envir=globalenv())
  load("ps.5.Rdata",envir=globalenv())
  load("ps.6.Rdata",envir=globalenv())
  load("ps.7.Rdata",envir=globalenv())
  load("data.Rdata",envir=globalenv())
  
  ps1 <- get("ps1")
  ps2 <- get("ps2")
  ps3 <- get("ps3")
  ps4 <- get("ps4")
  ps5 <- get("ps5")
  ps6 <- get("ps6")
  ps7 <- get("ps7")
  
  data <- get("data")
  
  
  ps <- lamda[1] * ps1 + lamda[2]* ps2 + lamda[3] * ps3 +
    lamda[4] * ps4 + lamda[5]* ps5 + lamda[6] * ps6+ lamda[7] * ps7
  data$ps <- ps
  treat <- data[data$Z == 1,]
  control <- data[data$Z == 0,]
  
  
  weight_treat <- 1/treat$ps
  weight_control <- 1/(1-control$ps)
  
  D <- c()
  for(i in 1:10){
    x <- treat[,i,drop =T]
    y <- control[,i,drop =T]
    a <-ks_test(x,y,w_x = weight_treat, w_y = weight_control, thresh = 0)$statistic
    D <- c(D,a)
  }
  mean_D <- mean(D)
}