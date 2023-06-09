library(IPWsurvival)
ad_km <- function(times = times, 
                  failures = failures, 
                  variable=variable, 
                  weights=weights,
                  tau = tau 
                  ){
  times[times >= tau] <- tau 
  failures[times >= tau] <- 0 
  res.akm <- adjusted.KM(times=times, failures=failures,
                         variable=variable, weights=weights)
  
  rate1 <- res.akm$survival[res.akm$variable == 1]
  rate2 <- res.akm$survival[res.akm$variable == 0]
  results_all <- rate1[length(rate1)]-rate2[length(rate2)]
  return(results_all)
}




####adjusted rmst####
akm_rmst <- function(time, status, group, weight=NULL, tau=NULL, alpha=.05, 
                     xaxismin=0, xaxismax=max(time)){ 
  if(sum(time<0)>0){print("Error: times must be positive.")
  }else{
    if(sum(weight<=0)>0){print("Error: weights must be greater than 0.")
    }else{
      if(sum(status!=0 & status!=1)>0){print("Error: status must be a vector of 0s and/or 1s.")
      }else{
        
        if(is.null(weight)){weight <- rep(1, length(time))}	
        data <- data.frame(time, status, group, weight)
        data <- data[!is.na(data$group) & !is.na(data$time),]
        data <- data[order(group),] 
        
        #--- If tau not specified, use minimum tau from all groups ---
        j=length(unique(data$group))
        
        if(is.null(tau)){
          taui = rep(999, j)
          for (i in (1:j)){
            groupval <- (levels(data$group)[i])
            dat_group <- data[which(data$group==(groupval)),]
            taui[i] <- max(dat_group$time[dat_group$status==1])
          }
          tau <- min(taui)
        }
        
        #--- Calculate AKM RMST in each group ---
        rmst <- rep(999, length(1:j))
        groupval <- rep(999, length(1:j))
        rmst_var <- rep(999, length(1:j))
        rmst_se <- rep(999, length(1:j))
 
        for (i in 1:j){
          groupval[i] <- (levels(data$group)[i])
          dat_group <- data[which(data$group==(groupval[i])),]
          
          #--- AKM ---
          tj <- c(0,sort(unique(dat_group$time[dat_group$status==1])))
          dj <- sapply(tj, function(x){sum(dat_group$weight[dat_group$time==x & dat_group$status==1])}) 
          yj <- sapply(tj, function(x){sum(dat_group$weight[dat_group$time>=x])}) 
          st <- cumprod(1-(dj/yj))
          m <- sapply(tj, function(x){sum((dat_group$weight[dat_group$time>=x])^2)})
          mj <- ((yj^2)/m)
          ft <- data.frame(tj, yj, dj, st, i, mj)
          
          #--- RMST ---

          rtime <- ft$tj<=tau
          tj_r <- sort(c(ft$tj[rtime],tau))
          st_r <- ft$st[rtime]
          yj_r <- ft$yj[rtime]
          dj_r <- ft$dj[rtime]
          time_diff <- diff(c(0, tj_r))
          areas <- time_diff * c(1, st_r)
          rmst[i] <- sum(areas)
          
          #--- Variance ---
          mj_r <- ft$mj[rtime]
          var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(mj_r *(yj_r - dj_r)))
          var_r <- c(var_r,0)
          rmst_var[i] <- sum(cumsum(rev(areas[-1]))^2 * rev(var_r)[-1])
          rmst_se[i] <- sqrt(rmst_var[i])

        }
      }
    }
  }
  
  results <- data.frame(groupval,rmst,rmst_var,rmst_se,tau)
  pwc <- ((j^2)-j)/2   #number of pairwise comparisons
  
  label_diff <- rep(999,(pwc))
  rmst_diff <- rep(999,(pwc))
  rmst_diff_se <- rep(999,(pwc))
  rmst_diff_low <- rep(999,(pwc))
  rmst_diff_upp <- rep(999,(pwc))
  rmst_diff_pval <- rep(999,(pwc))
  
  output_diff <- data.frame(label_diff,rmst_diff,rmst_diff_se,rmst_diff_low,rmst_diff_upp,rmst_diff_pval)
  l <- 1
  
  for (i in 1:(j-1)){
    for (j in (i+1):j){
      
      #--- RMST Difference ---
      output_diff[l,]$label_diff <- paste('Groups',results[j,]$groupval,'vs.',results[i,]$groupval,' ')
      output_diff[l,]$rmst_diff <- (results[j,]$rmst - results[i,]$rmst)
      output_diff[l,]$rmst_diff_se <- sqrt(results[j,]$rmst_var + results[i,]$rmst_var)
      output_diff[l,]$rmst_diff_low <- output_diff[l,]$rmst_diff - qnorm(1-alpha/2)*output_diff[l,]$rmst_diff_se
      output_diff[l,]$rmst_diff_upp <- output_diff[l,]$rmst_diff + qnorm(1-alpha/2)*output_diff[l,]$rmst_diff_se
      output_diff[l,]$rmst_diff_pval <- 2*(1-pnorm(abs(output_diff[l,]$rmst_diff)/output_diff[l,]$rmst_diff_se))
      
     
      l <- l+1 
    }
  }
  
  colnames(output_diff) <- c("Groups", "Est.", "SE", "CIL", "CIU", "p")
  rownames(output_diff) <- c(output_diff$Groups)
  result_all <- output_diff[c(2,3,4,5,6)]
  return(c(result_all$Est.,result_all$SE))
}

