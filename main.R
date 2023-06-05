setwd("~\\scenario19_result")
source("~\\func1.R")
source("~\\func2.R")
source("~\\func3.R")
source("~\\func4.R")


######################################################
#######single PS method estimation####################
######################################################
race_sample.2 <- race_sample(n.sample = 5000, ua=0,ub = 1, time = 30,censor = 50,n=c(1:500))
save(race_sample.2,file = "scenario19_sample_2.Rdata")


######################################################
#######model averaging method#########################
######################################################
################
###two better###
################
avg_sample.2.1 <- avg_all(n.sample = 5000, ua=0,ub = 1, time = 30,censor = 50,n=c(1:500),
                          model_type1 = "lg",model_type2 = "cbps")
save(avg_sample.2.1,file = "scenario19_avg_2_2.Rdata")

#############################
###one better and one bad####
#############################
avg_sample.2.2 <- avg_all(n.sample = 5000, ua=0,ub = 1, time = 30,censor = 50,n=c(1:500),
                          model_type1 = "lg",model_type2 = "cart")
save(avg_sample.2.2,file = "scenario19_avg_2_1.Rdata")

##############
###two bad####
##############
avg_sample.2.3 <- avg_all(n.sample = 5000, ua=0,ub = 1, time = 30,censor = 50,n=c(1:500),
                          model_type1 = "cart",model_type2 = "nn")
save(avg_sample.2.3,file = "scenario19_avg_2_0_1.Rdata")


#########################################
#####model averaging using 7 PS methods##
#########################################
avg_7_mean <- avg_all_7(n.sample = 1000, ua=0,ub = 1, time = 30,censor = 50,n=c(1:500),method_smd = "mean_smd")

avg_7_median <- avg_all_7(n.sample = 1000, ua=0,ub = 1, time = 30,censor = 50,n=c(1:500),method_smd = "median_smd")

avg_7_max <- avg_all_7(n.sample = 1000, ua=0,ub = 1, time = 30,censor = 50,n=c(1:500),method_smd = "max_smd")

avg_7_ks <- avg_all_7(n.sample = 1000, ua=0,ub = 1, time = 30,censor = 50,n=c(1:500),method_smd = "ks_smd")

avg_sample.m7 <- cbind(avg_7_mean,avg_7_median,avg_7_max,avg_7_ks)
save(avg_sample.m7,file = "scenario19_avg_7.Rdata")
apply(avg_sample.m7[,c(T,F)], 2, mean)
apply(avg_sample.m7[,c(T,F)], 2, mean)
