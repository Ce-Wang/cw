###########################################################
###############one model ##################################
###########################################################
race_sample <- function(
  n.sample = n.sample,
  ua=ua,ub=ub,
  time = time,  
  censor = censor,
  n = n){
  
  ##random forest
  print("rf")
  race_sample.rf <- c()
  for( i in n){
    cat(paste("####",i,"####\n"))
    b <- sample_race(n.sample = n.sample, ua=ua,ub = ub,model_type = "rf",time = time,seed =i,censor = censor) #model_type = rf
    race_sample.rf <- rbind(race_sample.rf,b)
  }
  
  ##logistic
  print("lg")
  race_sample.lg <- c()
  for( i in n){
    cat(paste("####",i,"####\n"))
    b <- sample_race(n.sample = n.sample, ua=ua,ub = ub,model_type = "lg",time = time,seed =i,censor = censor) #model_type = rf
    race_sample.lg <- rbind(race_sample.lg,b)
  }
  
  #svm_liner kernal
  print("svm_l")
  race_sample.svm_l <- c()
  for( i in n){
    cat(paste("####",i,"####\n"))
    b <- sample_race(n.sample = n.sample, ua=ua,ub = ub,model_type = "svm_l",time = time,seed =i,censor = censor)
    race_sample.svm_l <- rbind(race_sample.svm_l,b)
  }
  

  
  ##cart
  print("svm_cart")
  race_sample.cart <- c()
  for( i in n){
    cat(paste("####",i,"####\n"))
    b <- sample_race(n.sample = n.sample, ua=ua,ub = ub,model_type = "cart",time = time,seed =i,censor = censor) 
    race_sample.cart <- rbind(race_sample.cart,b)
  }
  
  
  # ##nn
  print("nn")
  race_sample.nn <- c()
  for( i in n){
    cat(paste("####",i,"####\n"))
    b <- sample_race(n.sample = n.sample, ua=ua,ub = ub,model_type = "nn",time = time,seed =i,censor = censor)
    race_sample.nn <- rbind(race_sample.nn,b)
  }

  ##cbps
  print("cbps")
  race_sample.cbps <- c()
  for( i in n){
    cat(paste("####",i,"####\n"))
    b <- sample_race(n.sample = n.sample, ua=ua,ub = ub,model_type = "cbps",time = time,seed =i,censor = censor) 
    race_sample.cbps <- rbind(race_sample.cbps,b)
  }
  
  # ##gbm
  print("gbm")
  race_sample.gbm <- c()
  for( i in n){
    cat(paste("####",i,"####\n"))
    b <- sample_race(n.sample = n.sample, ua=ua,ub = ub,model_type = "gbm",time = time,seed =i,censor = censor)
    race_sample.gbm <- rbind(race_sample.gbm,b)
  }
  

  race_sample.1 <- cbind(race_sample.lg,race_sample.cbps,
                         race_sample.cart,
                         race_sample.rf, race_sample.gbm,
                         race_sample.svm_l,
                         race_sample.nn)
  colnames(race_sample.1) <- c("lg","lg_se",
                               "cbps","cbps_se",
                               "cart","cart_se",
                               "rf","rf_se",
                               "gbm","gbm_se",
                               "svm_l","svm_l_se",
                               "nn","nn_se")
  return(race_sample.1)
}





###########################################################
###############model averaging function for two model #####
###########################################################
avg_all <- function(
  n.sample = n.sample,
  ua=ua,ub=ub,
  time = time,  
  censor = censor,
  n = n,
  model_type1 = model_type1,
  model_type2 = model_type2
){
  ##mean asmd
  mean_smd_all <- c()
  for( i in n){
    cat(paste("####",i,"####\n"))
    b <- avg_race(n.sample = n.sample, ua=ua,ub = ub,time = time,seed = i,censor =censor,model_type1 = model_type1,
                  model_type2 = model_type2,
                  method_smd = "mean_smd")
    mean_smd_all <- rbind(mean_smd_all,b)
  }
  
  ##median asmd
  median_smd_all <- c()
  for( i in n){
    cat(paste("####",i,"####\n"))
    b <- avg_race(n.sample = n.sample, ua=ua,ub = ub,time = time,seed = i,censor =censor,model_type1 = model_type1,
                  model_type2 = model_type2,
                  method_smd = "median_smd")
    median_smd_all <- rbind(median_smd_all,b)
  }
  
  ##max asmd
  max_smd_all <- c()
  for( i in n){
    cat(paste("####",i,"####\n"))
    b <- avg_race(n.sample = n.sample, ua=ua,ub = ub,time = time,seed = i,censor =censor,model_type1 = model_type1,
                  model_type2 = model_type2,
                  method_smd = "max_smd")
    max_smd_all <- rbind(max_smd_all,b)
  }
  #mean ks
  ks_all <- c()
  for( i in n){
    cat(paste("####",i,"####\n"))
    b <- avg_race(n.sample = n.sample, ua=ua,ub = ub,time = time,seed = i,censor =censor,model_type1 = model_type1,
                  model_type2 = model_type2,
                  method_smd = "ks")
    ks_all <- rbind(ks_all,b)
  }
  
  all_result <- cbind(mean_smd_all,median_smd_all, max_smd_all,ks_all)
  colnames(all_result) <- c("mean","mean_se","median","median_se","max","max_se","ks","ks_se")
  return(all_result)

}

###########################################################
###############7模型平均不同样本量输出结果#################
###########################################################
avg_all_7 <- function(
  n.sample = n.sample,
  ua=ua,ub=ub,
  time = time,  
  censor = censor,
  n = n,
  method_smd = "mean_smd"
){
  mean_smd_all <- c()
  indi <- c()
  for( i in n){
    tryCatch({
      cat(paste("####",i,"####\n"))
      b <- avg_race_m7(n.sample = n.sample, ua=ua,ub = ub,time = time,seed = i,censor =censor,
                       model_type1 = "lg",
                       model_type2 = "nn",
                       model_type3 = "svm_l",
                       model_type4 = "cart",
                       model_type5 = "rf",
                       model_type6 = "cbps",
                       model_type7 = "gbm",
                       method_smd = method_smd)
      mean_smd_all <- rbind(mean_smd_all,b)
    },
    error=function(e){
      indi <- c(indi,i)
      print(indi)
      cat("ERROR :",conditionMessage(e), "\n")
    })
  }
  return(mean_smd_all)
}
