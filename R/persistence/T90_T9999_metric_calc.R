##paramfile should have five columns: Dataset, Model, k1, k2, k3, BIC
##metofinterest should be T90, T99, T999, or T9999 with quotations
##maxtime is the expected maximum time for the metofinterest (default is 100 days)
metcalc <- function(paramfile, metofinterest, maxtime) {
  if (missing(maxtime)) {
    maxtime <- 100
  }
  param <- paramfile
  DataName <- as.character(param$Dataset)
  
  if (metofinterest=="T90") {
    met <- 1.00
  }
  if (metofinterest=="T99") {
    met <- 2.00
  }
  if (metofinterest=="T999") {
    met <- 3.00
  }
  if (metofinterest=="T9999") {
    met <- 4.00
  }
  if (metofinterest=="T99999") {
    met <- 5.00
  }
  ##Reading in two of the functions for the equations with if statements 
  pred_bi = function(k1,k2,t,LRV){	#k1 must be negative and k3 must be positive. k2 can be pos or neg. An increase in 2nd phase can happen if k2>k1.
    bp=72	#Breakpoint for the biphasic spline
    out=rep(0,length(t));
    if (min(t)<bp) out[t<bp] = exp(-k1*t[t<bp])
    if (max(t)>=bp) out[t>=bp] = exp(-k1*t[t>=bp]+k2*(t[t>=bp]-bp))
    if (LRV==TRUE) out=log10(out)
    return(out)
  }
  pred_bi3 = function(k1,k2,k3,t,LRV){	#k1 must be negative and k3 must be positive. k2 can be pos or neg. An increase in 2nd phase can happen if k2>k1. Reduces to exponential (with k1 as decay rate) if k3=0.
    out=0;
    if (min(t)<k3) out[t<k3] = exp(-k1*t[t<k3])
    if (max(t)>=k3) out[t>=k3] = exp(-k1*t[t>=k3]+k2*(t[t>=k3]-k3))
    if (LRV==TRUE) out=log10(out)
    return(out)
  }
  
  t=seq(0,maxtime, 0.01) #time vector, max currently is 200 days
  out <- array(data=NA, dim=c(length(t), length(param$Dataset)))
  for (i in 1:length(param$Dataset)) {
    if (param[i,2]=='ep') {
      out[,i] <- log10(exp(-param$k1[i]*t))
    }
    if (param[i,2]=='epd') {
      out[,i] <- log10(exp(-param$k1[i]*t*exp(-param$k2[i]*t)))
    }
    if (param[i,2]=='jm1') {
      out[,i] <- log10(1-(1-exp(-param$k1[i]*t))^param$k2[i])
    }
    if (param[i,2]=='jm2') {
      out[,i]  <- log10(1/(1+exp(param$k1[i]+param$k2[i]*log(t))))
    }
    if (param[i,2]=='lg1') {
      out[,i] <- log10(2*1/(1+exp(param$k1[i]*t)))
    }
    if (param[i,2]=='lg2') {
      out[,i] <- log10( 1/(1+exp(param$k1[i]*(t-param$k2[i]))))
    }
    if (param[i,2]=='gz') {
      out[,i]<- log10( exp(-param$k1[i]/param$k2[i]*(exp(param$k2[i]*t)-1))	)
    }
    if (param[i,2]=='gz3') {
      out[,i]<- log10(10^(param$k1[i]*exp(-exp((-param$k2[i]*exp(1)*(param$k3[i]-t)/param$k1[i])+1))))
    }
    if (param[i,2]=='gzm') {
      out[,i] <- log10(exp(-param$k3[i]*t - param$k1[i]/param$k2[i]*(exp(param$k2[i]*t)-1)))
    }
    if (param[i,2]=='wb') {
      out[,i] <- log10(10^(-((t/param$k1[i])^param$k2[i])))
    }
    if (param[i,2]=='dep') {
      out[,i]<- log10(param$k3[i] * exp(-param$k1[i]*t) + (1-param$k3[i]) * exp(-param$k2[i]*t)	)
    }
    if (param[i,2]=='gam') {
      out[,i]<- log10(1-pgamma(t,shape=param$k2[i],scale=param$k1[i])	)
    }
    if (param[i,2]=='ln') {
      out[,i] <- log10( 1-plnorm(t,meanlog=param$k1[i],sdlog=param$k2[i]))
    }
    if (param[i,2]=='sA') {
      out[,i] <- log10(10^-(param$k1[i]*t/((1+param$k2[i]*t)*(param$k3[i]-t))))
    }
    if (param[i,2]=='sB') {
      out[,i] <- log10(10^-(param$k1[i]*t^param$k3[i]/(param$k2[i]+t^param$k3[i])))
    }
    if(param[i,2]=='bi'){
      out[,i] <- pred_bi(param$k1[i], param$k2[i], t, TRUE)
    }
    if(param[i,2]=='bi3') {
      out[,i] <- pred_bi3(param$k1[i], param$k2[i], param$k3[i], t, TRUE)
    }
  }
  
  tt <- 0 ##identify the index of the time closest to a one log reduction
  for (j in 1:length(param$Dataset)) {
    tt[j] <- which.min(abs(out[,j]--met))
  }
  
  metval <- array(0, dim=c(length(tt),1))  ##Turn the index into an actual T90 time
  for (i in 1:length(tt)) {
    metval[i] <- t[tt[i]] 
  }

  metval  ##print out so you can see what they look like
  metval_final <- metval ##Combine back into one table
  
  ##See if we need to expand the time and run any datasets again####
  maxt <- which(metval==max(t) | metval==0) ##Determine if any T90 values were equal to the max time value above or 0. AKA their T90 is out of range.
  param2 <- param[maxt,]  ##Create a new dataframe with the datasets that need a long time frame
  t2 <- seq(0,maxtime*2,1) ##Expand the time frame to whatever you want (change the middle value)
  
  if (length(maxt)>0) {
    out2 <- array(data=NA, dim=c(length(t2), length(param2$Dataset)))
    for (i in 1:length(param2$Dataset)) {
      if (param2[i,2]=='ep') {
        out2[,i] <- log10(exp(-param2$k1[i]*t2))
      }
      if (param2[i,2]=='epd') {
        out2[,i] <- log10(exp(-param2$k1[i]*t2*exp(-param2$k2[i]*t2)))
      }
      if (param2[i,2]=='jm1') {
        out2[,i] <- log10(1-(1-exp(-param2$k1[i]*t2))^param2$k2[i])
      }
      if (param2[i,2]=='jm2') {
        out2[,i]  <- log10(1/(1+exp(param2$k1[i]+param2$k2[i]*log(t2))))
      }
      if (param2[i,2]=='lg1') {
        out2[,i] <- log10(2*1/(1+exp(param2$k1[i]*t2)))
      }
      if (param2[i,2]=='lg2') {
        out2[,i] <- log10( 1/(1+exp(param2$k1[i]*(t2-param2$k2[i]))))
      }
      if (param2[i,2]=='gz') {
        out2[,i]<- log10( exp(-param2$k1[i]/param2$k2[i]*(exp(param2$k2[i]*t2)-1))	)
      }
      if (param2[i,2]=='gz3') {
        out2[,i]<- log10(10^(param2$k1[i]*exp(-exp((-param2$k2[i]*exp(1)*(param2$k3[i]-t2)/param2$k1[i])+1))))
      }
      if (param2[i,2]=='gzm') {
        out2[,i] <- log10(exp(-param2$k3[i]*t2 - param2$k1[i]/param2$k2[i]*(exp(param2$k2[i]*t2)-1)))
      }
      if (param2[i,2]=='wb') {
        out2[,i] <- log10(10^(-((t2/param2$k1[i])^param2$k2[i])))
      }
      if (param2[i,2]=='dep') {
        out2[,i]<- log10(param2$k3[i] * exp(-param2$k1[i]*t2) + (1-param2$k3[i]) * exp(-param2$k2[i]*t2)	)
      }
      if (param2[i,2]=='gam') {
        out2[,i]<- log10(1-pgamma(t2,shape=param2$k2[i],scale=param2$k1[i])	)
      }
      if (param2[i,2]=='ln') {
        out2[,i] <- log10( 1-plnorm(t2,meanlog=param2$k1[i],sdlog=param2$k2[i]))
      }
      if (param2[i,2]=='sA') {
        out2[,i] <- log10(10^-(param2$k1[i]*t2/((1+param2$k2[i]*t2)*(param2$k3[i]-t2))))
      }
      if (param2[i,2]=='sB') {
        out2[,i] <- log10(10^-(param2$k1[i]*t2^param2$k3[i]/(param2$k2[i]+t2^param2$k3[i])))
      }
      if(param2[i,2]=='bi'){
        out2[,i] <- pred_bi(param2$k1[i], param2$k2[i], t2, TRUE)
      }
      if(param2[i,2]=='bi3') {
        out2[,i] <- pred_bi3(param2$k1[i], param2$k2[i], param2$k3[i], t2, TRUE)
      }
    } ##same for loop as above, referencing the new dataframe param2 and t2
    
    tt2 <- 0 ##identify the index of the time closest to a one log reduction
    for (j in 1:length(param2$Dataset)) {
      tt2[j] <- which.min(abs(out2[,j]--met))
    }
    
    metval_2 <- array(0, dim=c(length(tt2),1)) ##Turn the index into an actual T90 time
    for (i in 1:length(tt2)) {
      metval_2[i] <- t2[tt2[i]]
      
    }
    metval_2 ## print out so you can see what they look like 
    
    metval_final[maxt] <- metval_2 
    
  }
  
  ###ONe more time####
  maxt2 <- which(metval_final==max(t2) | metval_final==0) ##Determine if any T90 values were equal to the max time value above or 0. AKA their T90 is out of range.
  param3 <- param[maxt2,]  ##Create a new dataframe with the datasets that need a long time frame
  t3 <- seq(0,maxtime*3,1) ##Expand the time frame to whatever you want (change the middle value)
  
  if (length(maxt2)>0) {
    out3 <- array(data=NA, dim=c(length(t3), length(param3$Dataset)))
    for (i in 1:length(param3$Dataset)) {
      if (param3[i,2]=='ep') {
        out3[,i] <- log10(exp(-param3$k1[i]*t3))
      }
      if (param3[i,2]=='epd') {
        out3[,i] <- log10(exp(-param3$k1[i]*t3*exp(-param3$k2[i]*t3)))
      }
      if (param3[i,2]=='jm1') {
        out3[,i] <- log10(1-(1-exp(-param3$k1[i]*t3))^param3$k2[i])
      }
      if (param3[i,2]=='jm2') {
        out3[,i]  <- log10(1/(1+exp(param3$k1[i]+param3$k2[i]*log(t3))))
      }
      if (param3[i,2]=='lg1') {
        out3[,i] <- log10(2*1/(1+exp(param3$k1[i]*t3)))
      }
      if (param3[i,2]=='lg2') {
        out3[,i] <- log10( 1/(1+exp(param3$k1[i]*(t3-param3$k2[i]))))
      }
      if (param3[i,2]=='gz') {
        out3[,i]<- log10( exp(-param3$k1[i]/param3$k2[i]*(exp(param3$k2[i]*t3)-1))	)
      }
      if (param3[i,2]=='gz3') {
        out3[,i]<- log10(10^(param3$k1[i]*exp(-exp((-param3$k2[i]*exp(1)*(param3$k3[i]-t3)/param3$k1[i])+1))))
      }
      if (param3[i,2]=='gzm') {
        out3[,i] <- log10(exp(-param3$k3[i]*t3 - param3$k1[i]/param3$k2[i]*(exp(param3$k2[i]*t3)-1)))
      }
      if (param3[i,2]=='wb') {
        out3[,i] <- log10(10^(-((t3/param3$k1[i])^param3$k2[i])))
      }
      if (param3[i,2]=='dep') {
        out3[,i]<- log10(param3$k3[i] * exp(-param3$k1[i]*t3) + (1-param3$k3[i]) * exp(-param3$k2[i]*t3)	)
      }
      if (param3[i,2]=='gam') {
        out3[,i]<- log10(1-pgamma(t3,shape=param3$k2[i],scale=param3$k1[i])	)
      }
      if (param3[i,2]=='ln') {
        out3[,i] <- log10( 1-plnorm(t3,meanlog=param3$k1[i],sdlog=param3$k2[i]))
      }
      if (param3[i,2]=='sA') {
        out3[,i] <- log10(10^-(param3$k1[i]*t3/((1+param3$k2[i]*t3)*(param3$k3[i]-t3))))
      }
      if (param3[i,2]=='sB') {
        out3[,i] <- log10(10^-(param3$k1[i]*t3^param3$k3[i]/(param3$k2[i]+t3^param3$k3[i])))
      }
      if(param3[i,2]=='bi'){
        out3[,i] <- pred_bi(param3$k1[i], param3$k2[i], t3, TRUE)
      }
      if(param3[i,2]=='bi3') {
        out3[,i] <- pred_bi3(param3$k1[i], param3$k2[i], param3$k3[i], t3, TRUE)
      }
    } ##same for loop as above, referencing the new dataframe param2 and t2
    
    tt3 <- 0 ##identify the index of the time closest to a one log reduction
    for (j in 1:length(param3$Dataset)) {
      tt3[j] <- which.min(abs(out3[,j]--met))
    }
    
    metval_3 <- array(0, dim=c(length(tt3),1)) ##Turn the index into an actual T90 time
    for (i in 1:length(tt3)) {
      metval_3[i] <- t3[tt3[i]]
      
    }
    metval_3 ## print out so you can see what they look like 
   
    metval_final[maxt2] <- metval_3  
  }
 
  
  zerotrend <- which(metval_final==0 ) ##A T90 of "0" or the highest time value might mean that there was no observed decay in the timeframe, which means that the dataset
  stagnanttrend <- which(metval_final==max(t3))
  
  metval_print <- cbind.data.frame(param$Dataset, param$Model, round(metval_final,2) , param$BIC)
  colnames(metval_print) <- c("Dataset", "Model", metofinterest, "BIC")
  
  
  ##and print that table out as a csv file
  write.csv(metval_print, file.path("./Outputs", paste(metofinterest, " Values_v1.csv", sep='')))
  
  ###MOdel Averaging####
  num.cols <- c(metofinterest)
  metval_print[num.cols] <- sapply(metval_print[num.cols], as.numeric)
  modav_met <- 0; wv_ind <- 0; wv_sum <- 0; met_wv <- 0
  num_dat <- unique(metval_print$Dataset)
  for (i in 1:length(num_dat)) {
    
    xx <- metval_print[metval_print[,1]==num_dat[i],]
    
    for (j in 1:nrow(xx)) {
      wv_ind[j] <- exp(xx[j,4]/2)
    }
    wv_sum <- sum(wv_ind)
    rat <- wv_ind/wv_sum
    
    for (k in 1:nrow(xx)) {
      met_wv[k] <- xx[k,3]*rat[k]
    }
    
    modav_met[i] <- sum(met_wv)
    wv_ind <- 0
    met_wv <- 0
  }
  
  metval_MA <- cbind.data.frame(num_dat, round(modav_met, 2))
  colnames(metval_MA) <- c("Dataset", paste(metofinterest, "_MA", sep=''))
  write.csv(metval_MA, file.path("./Outputs", paste("Model Averaged ", metofinterest, "_v1.csv", sep='')))
  
  return(list("metric"=metval_print, "MA_metric"=metval_MA))
  
}
