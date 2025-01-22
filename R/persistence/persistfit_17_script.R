library(tidyr)
library(dplyr)
library(Metrics)

dataset <- n_trend[1]  # Add name of dataset here- note if the dataset is ABCdatalt in the environment, this line should read: dataset <- 'ABCdatalt'

eval(parse(text=paste("dat=",dataset,"[,c('Time','LR')]",sep=''))) # dat subsets the dataset, pulling the Time and LR columns
datname <- dataset  # Maintains the same dataset name
n <- nrow(dat)  # Number of observations in the dataset

if(n<3) {
  print("STOP- Insufficient observations for fitting all 17 models")
}

trend = lm(-LR ~ 0 + Time, data=dat)	#Forcing intercept through 0
Slope = summary(trend)$coefficients[1,1]
pSlope = summary(trend)$coefficients[1,4]
if (sign(Slope)==-1 & pSlope < 0.05){
  print("Significant negative trend")
} 
if (sign(Slope)==1 & pSlope < 0.05){
  print("STOP- Significant positive trend; the 17 models in this script are for decay")
} 
if (pSlope > 0.05){
  print("STOP- No significant trend")
}

path <- as.character(getwd())  # Name of the working directory
fold <- as.character(dataset)  # Makes the dataset name the name of the new folder
dir.create(file.path(path, 'Outputs', fold), recursive = TRUE, showWarnings = FALSE)  # Creates the dataset folder within "Outputs" in the directory

#### Prediction Functions ####  
## A function for each model predicting concentration over time (in days) ##
## Inputs are the model parameters (k1, k2, k3), time (t), and LRV status (if TRUE function outputs LRV over time) ##

## Exponential Model- (Chick-Watson) ##
pred_ep = function(k1,t,LRV){  
  #k1 = abs(k1)
  out = exp(-k1*t)
  if (LRV==TRUE) out=log10(out)
  return(out)
}
## Biphasic spline, preset breakpoint ##
pred_bi = function(k1,k2,t,LRV){	#k1 must be negative and k3 must be positive. k2 can be pos or neg. An increase in 2nd phase can happen if k2>k1.
  warning("This biphasic spline model ('bi') is not recommended; it has a preset breakpoint of 72. Use 'bi3' instead.")
  bp=72	#Breakpoint for the biphasic spline
  out=rep(0,length(t));
  if (min(t)<bp) out[t<bp] = exp(-k1*t[t<bp])
  if (max(t)>=bp) out[t>=bp] = exp(-k1*t[t>=bp]+k2*(t[t>=bp]-bp))
  if (LRV==TRUE) out=log10(out)
  return(out)
}
## Biphasic spline, with breakpoint as parameter k3 ##
pred_bi3 = function(k1,k2,k3,t,LRV){	#k1 must be negative and k3 must be positive. k2 can be pos or neg. An increase in 2nd phase can happen if k2>k1. Reduces to exponential (with k1 as decay rate) if k3=0.
  out=0;
  if (min(t)<k3) out[t<k3] = exp(-k1*t[t<k3])
  if (max(t)>=k3) out[t>=k3] = exp(-k1*t[t>=k3]+k2*(t[t>=k3]-k3))
  if (LRV==TRUE) out=log10(out)
  return(out)
}
## Exponential Damped ## Cavalli-Sforza 1983
pred_epd = function(k1,k2,t,LRV){
  out = exp(-k1*t*exp(-k2*t))	#k1 must be pos. & k2 neg. to have a monotonic decline. If k1 and k2 are both pos. the curve eventually returns to the t0 level. Reduces to exponential if k2=0. If k1 is neg. then there is no decline from the initial level.
  if (LRV==TRUE) out=log10(out)
  return(out)
}
## Juneja&Marks(1) (i.e. Multi-hit model) ##
pred_jm1 = function(k1,k2,t,LRV){	#k1 increases speed of decline. k2>1 gives shouldering & k2<1 gives tailing. Reduces to exponential if k2=1.
  #k1=abs(k1)
  #k2=abs(k2)
  out = 1-(1-exp(-k1*t))^k2
  if (LRV==TRUE) out=log10(out)
  return(out)
}
## Juneja&Marks(2) (based on log-logistic function) ##
pred_jm2 = function(k1,k2,t,LRV){	#Always tails. Neg. k1 is acceptable. Decreasing k2 leads to more rapid plateauing. Juneja VK 200, eq. 1?
  #k2=abs(k2)
  out = 1/(1+exp(k1+k2*log(t)))
  if (LRV==TRUE) out=log10(out)
  return(out)
}
## General logistic ##See Kamau DN 1990. The name is a misnomer (generalized logistic has 4 parameters).
pred_lg1 = function(k1,t,LRV){	#Slower decay than exponential for same value of k1; weak shouldering.
  out = 2*1/(1+exp(k1*t))	# 2 at beginning forces time of 0 to output 1. NOTE: check graphically if this is a simplified form of JM2.
  if (LRV==TRUE) out=log10(out)
  return(out)
}
## 2 parameter logistic ## Also known as a Fermi curve. See p.24, Peleg 2006, Advanced Quantitative Microbiology for Foods and Biosystems, CRC Press.
pred_lg2 = function(k1,k2,t,LRV){	#Slower decay than exponential for same value of k1; weak shouldering.
  #Only works if k2 >> k1; otherwise, Nt/N0 can be very different from 1 at time 0.
  out = 1/(1+exp(k1*(t-k2)))	#2 at beginning forces time of 0 to output 1. k2 is the time at the inflection point, and k1 is the decay rate at the inflection point.
  if (LRV==TRUE) out=log10(out)
  return(out)	#Only starts at 1 (or 0 on log scale) if k1 and k2 are both 'large', meaning that the product of k1 and k2 should be > ~7.
}
## Gompertz ##
pred_gz = function(k1,k2,t,LRV){	#Gives shouldering if k1 & k2 are both +, or cessation of decay at y=exp(-k1) if k1 & k2 are both -. Increasing k1 or k2 increases the rate of decline.
  #if (sign(k1) != sign(k2)) stop('k1 and k2 must have the same sign: both + for shouldering, both - for tailing.')
  #out = exp(-k1*exp(k2*t)+k1)	#k1 at end forces time of 0 to output 1; not sure where this formulation came from (original Tao Hong's code).
  #out = exp(-k1*(exp(k2*t)-1))	#From Wikipedia; 1 minus the cumulative form of the Gompertz distribution. Appears equivalent to version in Tao's code.
  out = exp(-k1/k2*(exp(k2*t)-1))	#Differs from Wikipedia; see Wu JW 2004, eq. 2. Also El-Gohary A 2013, eq. 1 with theta=1. Also 'gzm' with k3=0.
  if (LRV==TRUE) out=log10(out)
  return(out)
}
## Gompertz (3 parameter) ##	Gil 2011, eq. 14 (see also eq. 7). k1 is 'tail', k2 is max. inact. rate, k3 is 'shoulder'.
pred_gz3 = function(k1,k2,k3,t,LRV){	#Gives shouldering if k1 & k2 are both +, or cessation of decay at y=exp(-k1) if k1 & k2 are both -. Increasing k1 or k2 increases the rate of decline.
  out = 10^(k1*exp(-exp((-k2*exp(1)*(k3-t)/k1)+1)))	#Gil 2011. k1 is 'tail', k2 is max. inact. rate, k3 is 'shoulder'. Might be a convenient reparameterization, but != 0 at t=0.
  if (LRV==TRUE) out=log10(out)
  return(out)
}
## Gompertz-Makeham ##
pred_gzm = function(k1,k2,k3,t,LRV){	#k3 must be positive for decay.
  out = exp(-k3*t - k1/k2*(exp(k2*t)-1))	#Gompertz-Makeham CDF from Wikipedia. If k3 = 0, we get 2-parameter Gompertz.
  if (LRV==TRUE) out=log10(out)
  return(out)
}
## Weibull ##Coroller 2006. Other sources (e.g., Wikipedia) use e as the base instead of 10, but using 10 means that k1 conveniently equals the T90.
pred_wb = function(k1,k2,t,LRV){	#k1 is the time needed for the 1st log10 reduction (delta); k2 is a shape parameter (P). Both must be positive.
  out = 10^(-((t/k1)^k2))			#If k2 < 1, there is tailing; if k2 > 1, there is shouldering.
  if (LRV==TRUE) out=log10(out)	#Equivalent to exponential model if k2==1. k1 for the equivalent exp. model then equals -log(10^(-1/k1)).
  return(out)
}
## Double exponential ##See p.6, Peleg 2006, Advanced Quantitative Microbiology for Foods and Biosystems, CRC Press. See also GInaFiT.
pred_dep = function(k1,k2,k3,t,LRV){	#k1 and k2 are decay rates for each of 2 microbe populations. k3 is the proportion of the total population subject to k1.
  out = k3 * exp(-k1*t) + (1-k3) * exp(-k2*t)	#But see Abraham 1990: allowing k3 to go above 1 gives 'activation shoulder' behavior.
  if (LRV==TRUE) out=log10(out)	#Equivalent to exponential model if k2==1. k1 for the equivalent exp. model then equals -log(10^(-1/k1)).
  return(out)
}
## Gamma ##
pred_gam = function(k1,k2,t,LRV){
  out = 1-pgamma(t,shape=k2,scale=k1)	#If shape=1 (k2=1), reduces to exponential.			
  if (LRV==TRUE) out=log10(out)
  return(out)
}
## Lognormal ## Aragao 2007. Behaves similarly to the Weibull model, but becomes loglinear sooner if there is shouldering.
pred_ln = function(k1,k2,t,LRV){
  out = 1-plnorm(t,meanlog=k1,sdlog=k2)			
  if (LRV==TRUE) out=log10(out)
  return(out)
}
## Sigmoid A (fast-slow-fast) ##See p.35, Peleg 2006, Advanced Quantitative Microbiology for Foods and Biosystems, CRC Press.
pred_sA = function(k1,k2,k3,t,LRV){
  out = 10^-(k1*t/((1+k2*t)*(k3-t)))
  if (LRV==TRUE) out=log10(out)
  return(out)
}
## Sigmoid B (slow-fast-slow) ##See p.36, Peleg 2006, Advanced Quantitative Microbiology for Foods and Biosystems, CRC Press.
pred_sB = function(k1,k2,k3,t,LRV){
  out = 10^-(k1*t^k3/(k2+t^k3))	#k1 appears to be the max LRV.
  if (LRV==TRUE) out=log10(out)
  return(out)
}


#### Likelihood functions ####
## A function for each 1-parameter, 2-parameter, and 3-parameter model ##
## Inputs are the model parameters (k1, k2, k3, model sd), data with Time and LR columns, and abbreviated model name ##
## Abbreviated model name is used to call the respective prediction function within the for loop 

mod_L1 = function(pars,data, mod){ #Common core of all 1-parameter likelihood functions for persistence models
  k1=pars[1]
  sd=pars[2]
  t=data[,1]
  LR=log(10^(-data[,2]))
  model = mod
  lik_temp = matrix(nrow=NROW(LR),ncol=length(k1))
  for (i in 1:length(k1)){
    as.numeric(eval(parse(text=paste("lik_temp[,i]=dnorm(LR,log(pred_", model,"(k1[i],t,LRV=FALSE)),sd,log=TRUE)", sep="" ))))
  }
  lik=-colSums(lik_temp) 
  lik  
}

mod_L2 = function(pars,data, mod){ #Common core of all 2-parameter likelihood functions for persistence models
  k1=pars[1]
  k2=pars[2]
  sd=pars[3]
  t=data[,1]
  LR=log(10^(-data[,2]))
  model = mod
  ncols=1	#Default value corresponding to single values for k1 and k2; gets updated if multiple values for k1 or k2 are being used.
  
  lik_temp = matrix(nrow=NROW(LR),ncol=ncols)
  for (i in 1:ncols){
    as.numeric(eval(parse(text=paste("lik_temp[,i]=dnorm(LR,log(pred_", model,"(k1[i], k2[i], t,LRV=FALSE)),sd,log=TRUE)", sep=""))))
  }
  lik=-colSums(lik_temp)
  lik
}

mod_L3 = function(pars,data, mod){ #Common core of all 3-parameter likelihood functions for persistence models
  k1=pars[1]
  k2=pars[2]
  k3=pars[3]
  sd=pars[4]
  t=data[,1]
  LR=log(10^(-data[,2]))
  model=mod
  ncols=1	#Default value corresponding to single values for k1 and k2; gets updated if multiple values for k1 or k2 are being used.
  lik_temp = matrix(nrow=NROW(LR),ncol=ncols)
  for (i in 1:ncols){
    as.numeric(eval(parse(text=paste("lik_temp[,i]=dnorm(LR,log(pred_", model,"(k1[i], k2[i], k3[i],t,LRV=FALSE)),sd,log=TRUE)", sep=""))))
  }
  lik=-colSums(lik_temp)
  lik
}

#### Identification and Transformations ####
## Assign unique color+line-type combination to each model for graphing purposes
GSs = data.frame('model'=c('ep','lg1','lg2','epd','jm1','jm2','gz','wb',
                           'ln','gam','bi','bi3','dep','gz3','gzm','sA','sB'),
                 'color'= c('black','red', 'orange', 'yellow', 'green', 'blue', 'purple',
                                    'red', 'orange', 'yellow', 'green', 'blue', 'purple',
                                    'red', 'orange', 'yellow', 'green'),
                 'lty' = c('solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid',
                           'dashed', 'dashed', 'dashed', 'dashed', 'dashed', 'dashed',
                           'dotted', 'dotted', 'dotted', 'dotted'))
row.names(GSs) = c('ep','lg1','lg2','epd','jm1','jm2','gz','wb','ln','gam','bi','bi3','dep','gz3','gzm','sA','sB')
GSs[,c('k1','k2','k3','k4','k5')] = '-'	#Making fields for transformation info

## Transforming some parameters to improve model fitting efficiency
GSs[,'suffix'] = ''	#Suffix for likelihood function, denoting transformation. Used later on in the gutsMLEx() functions.
GSs['lg2',c('k1','k2','suffix')] = c('real2pos','real2pos','.T1')
GSs['jm1',c('k1','k2','suffix')] = c('real2pos','real2pos','.T1')
GSs['gz',c('k1','k2','suffix')] = c('real2extreme','real2extreme','.T1')
GSs['gam',c('k1','suffix')] = c('real2extreme','.T1')
GSs['bi3',c('k3','suffix')] = c('real2pos','.T1')
GSs['dep',c('k3','suffix')] = c('real2extreme','.T2')
GSs['gz3',c('k2','k3','suffix')] = c('real2pos','real2pos','.T1')
GSs['gzm',c('k1','k2','k3','suffix')] = c('real2extreme','real2extreme','real2extreme','.T1')

real2pos = function(x,reverse){  #Seems obvious, but implementing it this way for consistency.
  if (reverse==1){
    y=log(x)
  }else{
    y=exp(x)
  }
  return(y)
}
real2extreme = function(x,reverse){	#In reverse, makes extremely huge- or tiny-magnitude numbers closer to 1, preserving the sign.
  if (reverse==1){
    y = nthroot(x,7)	#Real 7th root. Could change 7 to some other odd integer if desired.
  }else{
    y = x^7
  }
}

#### Starting Parameters - Creating a dataframe of starting guessses for parameters in each model  ####
DesiredModels = c('ep','lg1','lg2','epd','jm1','jm2','gz','wb','ln','gam','bi','bi3','dep','gz3','gzm','sA','sB')
start = data.frame('model'=DesiredModels,'k1'=NA,'k2'=NA,'k3'=NA,'k4'=NA,'k5'=NA,'sd'=NA)	#Starting values for parameters. These are always untransformed.
row.names(start) = DesiredModels
start = start[,-1]	#Removing the now-useless start$model column
start$sd = 1
start$MLEmethod = 'Nelder-Mead'	#The optimization method chosen at the top of this code. Applied to all models.
start['wb','MLEmethod'] = 'SANN'
if (is.na(start['ep','k1'])) start['ep','k1']   		= 2e-4
if (is.na(start['lg1','k1'])) start['lg1','k1'] 			= 1e-4
if (is.na(start['lg2','k1'])) start['lg2',c('k1','k2')] 	= c(1e-4,1e4)
if (is.na(start['epd','k1'])) start['epd',c('k1','k2')] 	= c(1e-3,1e-3)	#Gets screwy if k2 < 1e-3
if (is.na(start['jm1','k1'])) start['jm1',c('k1','k2')] 	= c(1e-5,1)
if (is.na(start['jm2','k1'])) start['jm2',c('k1','k2')]	= c(-0.9,0.01)
if (is.na(start['gz','k1'])) start['gz',c('k1','k2')] 	= c(-1,-1.38)
if (is.na(start['wb','k1'])) start['wb',c('k1','k2')] = c(15,0.5)
if (is.na(start['ln','k1'])) start['ln',c('k1','k2')] 	= c(10,5)
if (is.na(start['gam','k1'])) start['gam',c('k1','k2')] 	= c(5,0.1)
if (is.na(start['bi','k1'])) start['bi',c('k1','k2')] 	= c(7e-3,7e-3)
if (is.na(start['bi3','k1'])) start['bi3',c('k1','k2','k3')] 	= c(1e-4,1e-4,96); #startT['bi3','k3'] = 'real2pos'
if (is.na(start['dep','k1'])) start['dep',c('k1','k2','k3')] 	= c(1e-3,1e-5,1.06); #startT['dep','k3'] = 'real2prob'	#NOTE: startT must agree with the transformations used in the nll.XXX.T() functions.
if (is.na(start['gz3','k1'])) start['gz3',c('k1','k2','k3')] 	= c(-0.1,1e-4,1e1)
if (is.na(start['gzm','k1'])) start['gzm',c('k1','k2','k3')] 	= c(0.25, -0.01, 0.10)
if (is.na(start['sA','k1'])) start['sA',c('k1','k2','k3')] 	= c(-1e3,0.0001,500)
if (is.na(start['sB','k1'])) start['sB',c('k1','k2','k3')] 	= c(0.1,1e6,2)
startT = start
outputT = startT	# Output dataframe looks like starting parameter dataframe
outputT$nll = NA  # Add a column from negative log-likelihood value
outputT$BIC = NA  # Add a column for BIC value
outputT$kSource = 'start'	 # Indicate that the parameter estimates inititated from the starting guesses

####First Optimization ####
conv1 <- array(data=NA, dim=c(length(DesiredModels), 1)) # Setting up an empty dataframe for the convergence diagnostic
rownames(conv1) <- DesiredModels
for (m in 1:length(DesiredModels)) {
  mod <- DesiredModels[m]
  if (mod %in% c(c("ep", "lg1"))) {
    tryCatch(eval(parse(text=paste(mod, "_results <- optim(par=c(k1=start['", mod,"','k1'], sd=start['", mod,
                                   "','sd']), mod='", mod, "', fn=mod_L1, data=dat, hessian=T)", sep=""))), 
             error=function(cond) {paste("optim with hessian=T failed")})
    
    if(exists(paste(mod, "_results", sep='')) & is.list(eval(parse(text=paste(mod, "_results", sep=''))))) {
      eval(parse(text=paste(mod, "_results =", mod, "_results", sep="")))
      eval(parse(text=paste("conv1['", mod, "', 1] <- ", mod,"_results$convergence", sep="")))
      eval(parse(text=paste("outputT['", mod, "', 'k1'] <- ", mod,"_results$par[1]", sep="")))
      eval(parse(text=paste("outputT['", mod, "', 'sd'] <- ", mod,"_results$par[2]", sep="")))
      eval(parse(text=paste("outputT['", mod, "', 'nll'] <- ", mod,"_results$value", sep="")))
      
    } else {
      eval(parse(text=paste(mod, "_results <- optim(par=c(k1=start['", mod,"','k1'], sd=start['", mod,
                            "','sd']), mod='", mod, "', fn=mod_L1, data=dat, hessian=F)", sep="")))
      eval(parse(text=paste("conv1['", mod, "', 1] <- ", mod,"_results$convergence", sep="")))
      eval(parse(text=paste("outputT['", mod, "', 'k1'] <- ", mod,"_results$par[1]", sep="")))
      eval(parse(text=paste("outputT['", mod, "', 'sd'] <- ", mod,"_results$par[2]", sep="")))
      eval(parse(text=paste("outputT['", mod, "', 'nll'] <- ", mod,"_results$value", sep="")))
     }
    
  }
  if(mod %in% c('lg2','epd','jm1','jm2','gz','wb','ln','gam','bi')){
    tryCatch(eval(parse(text=paste(mod, "_results <- optim(par=c(k1=start['", mod,"','k1'], k2=start['", mod, "', 'k2'],
                               sd=start['", mod, "','sd']), mod='", mod, "', fn=mod_L2, data=dat, hessian=T)", sep=""))), 
             error=function(cond) {paste("optim with hessian=T failed")})
    
    
    if(exists(paste(mod, "_results", sep=''))){
      if(is.list(eval(parse(text=paste(mod, "_results", sep=''))))) { 
        eval(parse(text=paste(mod, "_results =", mod, "_results", sep="")))
        eval(parse(text=paste("conv1['", mod, "', 1] <- ", mod,"_results$convergence", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'k1'] <- ", mod,"_results$par[1]", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'k2'] <- ", mod,"_results$par[2]", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'sd'] <- ", mod,"_results$par[3]", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'nll'] <- ", mod,"_results$value", sep="")))
       }} else {
        eval(parse(text=paste(mod, "_results <- optim(par=c(k1=start['", mod,"','k1'], k2=start['", mod, "', 'k2'],
                               sd=start['", mod, "','sd']), mod='", mod, "', fn=mod_L2, data=dat, hessian=F)", sep="")))
        eval(parse(text=paste("conv1['", mod, "', 1] <- ", mod,"_results$convergence", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'k1'] <- ", mod,"_results$par[1]", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'k2'] <- ", mod,"_results$par[2]", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'sd'] <- ", mod,"_results$par[3]", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'nll'] <- ", mod,"_results$value", sep="")))
       } 
  }
  if (mod %in% c('bi3', 'dep', 'gz3', 'gzm', 'sA', 'sB')) {
    tryCatch(eval(parse(text=paste(mod, "_results <- optim(par=c(k1=start['", mod,"','k1'], k2=start['", mod, "', 'k2'], 
  k3=start['", mod, "', 'k3'], sd=start['", mod, "','sd']), mod='", mod, "', fn=mod_L3, method='L-BFGS-B',lower=lower_", mod, ", 
                                 upper=upper_", mod, ", data=dat, hessian=T)", sep=""))), 
             error=function(cond) {paste("optim with hessian=T failed")})
    
    
    if(exists(paste(mod, "_results", sep=''))){
      if(is.list(eval(parse(text=paste(mod, "_results", sep=''))))) { 
        eval(parse(text=paste(mod, "_results =", mod, "_results", sep="")))
        eval(parse(text=paste("conv1['", mod, "', 1] <- ", mod,"_results$convergence", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'k1'] <- ", mod,"_results$par[1]", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'k2'] <- ", mod,"_results$par[2]", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'k3'] <- ", mod,"_results$par[3]", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'sd'] <- ", mod,"_results$par[4]", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'nll'] <- ", mod,"_results$value", sep="")))
      }} else {
        eval(parse(text=paste(mod, "_results <- optim(par=c(k1=start['", mod,"','k1'], k2=start['", mod, "', 'k2'], k3=start['", mod, "', 'k3'],
                               sd=start['", mod, "','sd']), mod='", mod, "', fn=mod_L3, data=dat, hessian=F)", sep="")))
        eval(parse(text=paste("conv1['", mod, "', 1] <- ", mod,"_results$convergence", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'k1'] <- ", mod,"_results$par[1]", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'k2'] <- ", mod,"_results$par[2]", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'k3'] <- ", mod,"_results$par[3]", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'sd'] <- ", mod,"_results$par[4]", sep="")))
        eval(parse(text=paste("outputT['", mod, "', 'nll'] <- ", mod,"_results$value", sep="")))
      } 
  }
}

## The following lines add all of the gathered data (parameters, standard errors, convergence and nll values into the results dataframe)
results <- as.data.frame(matrix(NA, nrow = 17, ncol = 8))
results[,1:3] <- round(outputT[,1:3],4) # Rounded parameter estimates
results[,4:7] <- outputT[, 6:9]  # sd, method, nll, BIC
results[,4] <- round(results[,4],4) # Rounded model sd estimate
results[,7] <- round(results[,7],2) # Rounded BIC estimate
results[,8] <- conv1  # Convergence indicator

## Calculate BIC values for each model - 2*nll + k*ln(n) where nll is the negative log-likelihood, k is the number of parameters in the likelihood function, and n is the number of data-points
BIC <- c(2*results[1,6]+2*log(nrow(dat)), 2*results[2,6]+2*log(nrow(dat)), 2*results[3,6]+3*log(nrow(dat)),
         2*results[4,6]+3*log(nrow(dat)), 2*results[5,6]+3*log(nrow(dat)), 2*results[6,6]+3*log(nrow(dat)),
         2*results[7,6]+3*log(nrow(dat)), 2*results[8,6]+3*log(nrow(dat)), 2*results[9,6]+3*log(nrow(dat)),
         2*results[10,6]+3*log(nrow(dat)), 2*results[11,6]+3*log(nrow(dat)),
         2*results[12,6]+4*log(nrow(dat)), 2*results[13,6]+4*log(nrow(dat)), 2*results[14,6]+4*log(nrow(dat)),
         2*results[15,6]+4*log(nrow(dat)), 2*results[16,6]+4*log(nrow(dat)), 2*results[17,6]+4*log(nrow(dat)))

results[,7] <- round(BIC,2)  # Rounded BIC values


colnames(results) <- c("k1", "k2", "k3", "sd", "Method", "nll", "BIC",  "Converged?")
row.names(results) <- c('ep','lg1','lg2','epd','jm1','jm2','gz','wb','ln','gam','bi','bi3','dep','gz3','gzm','sA','sB')

write.csv(results, file.path(path, 'Outputs', fold, paste(datname, '_Output1', '.csv', sep=""))) # Saves first optimization in Outputs folder

#### Second Optimization - Optimized parameters from first optimization are now the starting parameters ####
convmid <- array(data=NA, dim=c(length(DesiredModels), 1)) ##setting up an empty dataframe for the convergence diagnostic
rownames(convmid) <- DesiredModels
outputmid <- startT

for (m in 1:length(DesiredModels)) {
  mod <- DesiredModels[m]
  if (mod %in% c("ep", "lg1")) {
    tryCatch(eval(parse(text=paste(mod, "_results_mid <- optim(par=c(k1=outputT['", mod,"','k1'], sd=outputT['", mod,
                                   "','sd']), mod='", mod, "', fn=mod_L1, data=dat, hessian=T)", sep=""))), 
             error=function(cond) {paste("optim with hessian=T failed")})
    
    if(exists(paste(mod, "_results_mid", sep='')) & is.list(eval(parse(text=paste(mod, "_results_mid", sep=''))))) {
      eval(parse(text=paste(mod, "_results_mid  =", mod, "_results_mid ", sep="")))
      eval(parse(text=paste("convmid['", mod, "', 1] <- ", mod,"_results_mid$convergence", sep="")))
      eval(parse(text=paste("outputmid['", mod, "', 'k1'] <- ", mod,"_results_mid$par[1]", sep="")))
      eval(parse(text=paste("outputmid['", mod, "', 'sd'] <- ", mod,"_results_mid$par[2]", sep="")))
      eval(parse(text=paste("outputmid['", mod, "', 'nll'] <- ", mod,"_results_mid$value", sep="")))
    } else {
      eval(parse(text=paste(mod, "_results_mid <- optim(par=c(k1=outputT['", mod,"','k1'], sd=outputT['", mod,
                            "','sd']), mod='", mod, "', fn=mod_L1, data=dat, hessian=F)", sep="")))
      eval(parse(text=paste("convmid['", mod, "', 1] <- ", mod,"_results_mid$convergence", sep="")))
      eval(parse(text=paste("outputmid['", mod, "', 'k1'] <- ", mod,"_results_mid$par[1]", sep="")))
      eval(parse(text=paste("outputmid['", mod, "', 'sd'] <- ", mod,"_results_mid$par[2]", sep="")))
      eval(parse(text=paste("outputmid['", mod, "', 'nll'] <- ", mod,"_results_mid$value", sep="")))
    }  
    
  }
  if (mod %in% c('lg2','epd','jm1','jm2','gz','wb','ln','gam','bi')) {
    tryCatch(eval(parse(text=paste(mod, "_results_mid <- optim(par=c(k1=outputT['", mod,"','k1'], k2=outputT['", mod, "', 'k2'],
                               sd=outputT['", mod, "','sd']), mod='", mod, "', fn=mod_L2, data=dat, hessian=T)", sep=""))), 
             error=function(cond) {paste("optim with hessian=T failed")})
    
    
    if(exists(paste(mod, "_results_mid", sep=''))){
      if(is.list(eval(parse(text=paste(mod, "_results_mid", sep=''))))) {
        eval(parse(text=paste(mod, "_results_mid =", mod, "_results_mid", sep="")))
        eval(parse(text=paste("convmid['", mod, "', 1] <- ", mod,"_results_mid$convergence", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'k1'] <- ", mod,"_results_mid$par[1]", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'k2'] <- ", mod,"_results_mid$par[2]", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'sd'] <- ", mod,"_results_mid$par[3]", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'nll'] <- ", mod,"_results_mid$value", sep="")))
       }} else {
        eval(parse(text=paste(mod, "_results_mid <- optim(par=c(k1=outputT['", mod,"','k1'], k2=outputT['", mod, "', 'k2'],
                               sd=outputT['", mod, "','sd']), mod='", mod, "', fn=mod_L2, data=dat, hessian=F)", sep="")))
        
        eval(parse(text=paste("convmid['", mod, "', 1] <- ", mod,"_results_mid$convergence", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'k1'] <- ", mod,"_results_mid$par[1]", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'k2'] <- ", mod,"_results_mid$par[2]", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'sd'] <- ", mod,"_results_mid$par[3]", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'nll'] <- ", mod,"_results_mid$value", sep="")))
      } 
  }
  if (mod %in% c('bi3', 'dep', 'gz3', 'gzm', 'sA', 'sB')) {
    eval(parse(text=paste("tryCatch(", mod, "_results_mid <- optim(par=c(k1=outputT['", mod,"','k1'], k2=outputT['", mod, "', 'k2'], 
  k3=outputT['", mod, "', 'k3'], sd=outputT['", mod, "','sd']), mod='", mod, "', fn=mod_L3, method='L-BFGS-B', lower=lower_", mod, ", 
                                 upper=upper_", mod, ", data=dat, hessian=T), 
                        error=function(cond) {paste('optim with hessian=T failed')})", sep=""))) 
    
    
    
    if(exists(paste(mod, "_results_mid", sep=''))){
      if(is.list(eval(parse(text=paste(mod, "_results_mid", sep=''))))) { 
        eval(parse(text=paste(mod, "_results_mid =", mod, "_results_mid", sep="")))
        eval(parse(text=paste("convmid['", mod, "', 1] <- ", mod,"_results_mid$convergence", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'k1'] <- ", mod,"_results_mid$par[1]", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'k2'] <- ", mod,"_results_mid$par[2]", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'k3'] <- ", mod,"_results_mid$par[3]", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'sd'] <- ", mod,"_results_mid$par[4]", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'nll'] <- ", mod,"_results_mid$value", sep="")))
       }} else {
        eval(parse(text=paste(mod, "_results_mid <- optim(par=c(k1=outputT['", mod,"','k1'], k2=outputT['", mod, "', 'k2'], k3=outputT['", mod, "', 'k3'],
                               sd=outputT['", mod, "','sd']), mod='", mod, "', fn=mod_L3, data=dat, hessian=F)", sep="")))
        eval(parse(text=paste("convmid['", mod, "', 1] <- ", mod,"_results_mid$convergence", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'k1'] <- ", mod,"_results_mid$par[1]", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'k2'] <- ", mod,"_results_mid$par[2]", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'k3'] <- ", mod,"_results_mid$par[3]", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'sd'] <- ", mod,"_results_mid$par[4]", sep="")))
        eval(parse(text=paste("outputmid['", mod, "', 'nll'] <- ", mod,"_results_mid$value", sep="")))
       } 
  }
  
}

###Third Optimization- Optimized parameters from second optimization are now the starting parameters ####
conv2 <- array(data=NA, dim=c(length(DesiredModels), 1)) ##setting up an empty dataframe for the convergence diagnostic
rownames(conv2) <- DesiredModels
output2 <- outputT
for (m in 1:length(DesiredModels)) {
  mod <- DesiredModels[m]
  if (mod %in% c("ep", "lg1")) {
    eval(parse(text=paste("tryCatch(", mod, "_results2 <- optim(par=c(k1=outputmid['", mod,"','k1'], sd=outputmid['", mod,
                          "','sd']), mod='", mod, "', fn=mod_L1, data=dat, hessian=T), 
                                 error=function(cond) {paste('optim with hessian=T failed')})", sep=""))) 
    
    if(exists(paste(mod, "_results2", sep='')) & is.list(eval(parse(text=paste(mod, "_results2", sep=''))))) {
      eval(parse(text=paste(mod, "_results2  =", mod, "_results2 ", sep="")))
      eval(parse(text=paste("conv2['", mod, "', 1] <- ", mod,"_results2$convergence", sep="")))
      eval(parse(text=paste("output2['", mod, "', 'k1'] <- ", mod,"_results2$par[1]", sep="")))
      eval(parse(text=paste("output2['", mod, "', 'sd'] <- ", mod,"_results2$par[2]", sep="")))
      eval(parse(text=paste("output2['", mod, "', 'nll'] <- ", mod,"_results2$value", sep="")))
    } else {
      eval(parse(text=paste(mod, "_results2 <- optim(par=c(k1=outputmid['", mod,"','k1'], sd=outputmid['", mod,
                            "','sd']), mod='", mod, "', fn=mod_L1, data=dat, hessian=F)", sep="")))
      eval(parse(text=paste("conv2['", mod, "', 1] <- ", mod,"_results2$convergence", sep="")))
      eval(parse(text=paste("output2['", mod, "', 'k1'] <- ", mod,"_results2$par[1]", sep="")))
      eval(parse(text=paste("output2['", mod, "', 'sd'] <- ", mod,"_results2$par[2]", sep="")))
      eval(parse(text=paste("output2['", mod, "', 'nll'] <- ", mod,"_results2$value", sep="")))
      }
  }
  if (mod %in% c('lg2','epd','jm1','jm2','gz','wb','ln','gam','bi')) {
    eval(parse(text=paste("tryCatch(", mod, "_results2 <- optim(par=c(k1=outputmid['", mod,"','k1'], k2=outputmid['", mod, "', 'k2'],
                               sd=outputmid['", mod, "','sd']), mod='", mod, "', fn=mod_L2, data=dat, hessian=T), 
                                 error=function(cond) {paste('optim with hessian=T failed')})", sep=""))) 
    
    
    if(exists(paste(mod, "_results2", sep=''))){
      if(is.list(eval(parse(text=paste(mod, "_results2", sep=''))))) {
        eval(parse(text=paste(mod, "_results2 =", mod, "_results2", sep="")))
        eval(parse(text=paste("conv2['", mod, "', 1] <- ", mod,"_results2$convergence", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'k1'] <- ", mod,"_results2$par[1]", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'k2'] <- ", mod,"_results2$par[2]", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'sd'] <- ", mod,"_results2$par[3]", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'nll'] <- ", mod,"_results2$value", sep="")))
       
      }} else {
        eval(parse(text=paste(mod, "_results2 <- optim(par=c(k1=outputmid['", mod,"','k1'], k2=outputmid['", mod, "', 'k2'],
                               sd=outputmid['", mod, "','sd']), mod='", mod, "', fn=mod_L2, data=dat, hessian=F)", sep="")))
        
        eval(parse(text=paste("conv2['", mod, "', 1] <- ", mod,"_results2$convergence", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'k1'] <- ", mod,"_results2$par[1]", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'k2'] <- ", mod,"_results2$par[2]", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'sd'] <- ", mod,"_results2$par[3]", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'nll'] <- ", mod,"_results2$value", sep="")))
        
      } 
  }
  if (mod %in% c('bi3', 'dep', 'gz3', 'gzm', 'sA', 'sB')) {
    eval(parse(text=paste("tryCatch(", mod, "_results2 <- optim(par=c(k1=outputmid['", mod,"','k1'], k2=outputmid['", mod, "', 'k2'], 
  k3=outputmid['", mod, "', 'k3'], sd=outputmid['", mod, "','sd']), mod='", mod, "', fn=mod_L3, method='L-BFGS-B', lower=lower_", mod, ", 
                                 upper=upper_", mod, ", data=dat, hessian=T), 
                        error=function(cond) {paste('optim with hessian=T failed')})", sep=""))) 
    
    if(exists(paste(mod, "_results2", sep=''))){
      if(is.list(eval(parse(text=paste(mod, "_results2", sep=''))))) { 
        eval(parse(text=paste(mod, "_results2 =", mod, "_results2", sep="")))
        eval(parse(text=paste("conv2['", mod, "', 1] <- ", mod,"_results2$convergence", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'k1'] <- ", mod,"_results2$par[1]", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'k2'] <- ", mod,"_results2$par[2]", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'k3'] <- ", mod,"_results2$par[3]", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'sd'] <- ", mod,"_results2$par[4]", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'nll'] <- ", mod,"_results2$value", sep="")))
       
      }} else {
        eval(parse(text=paste(mod, "_results2 <- optim(par=c(k1=outputmid['", mod,"','k1'], k2=outputmid['", mod, "', 'k2'], k3=outputmid['", mod, "', 'k3'],
                               sd=outputmid['", mod, "','sd']), mod='", mod, "', fn=mod_L3, data=dat, hessian=F)", sep="")))
        eval(parse(text=paste("conv2['", mod, "', 1] <- ", mod,"_results2$convergence", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'k1'] <- ", mod,"_results2$par[1]", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'k2'] <- ", mod,"_results2$par[2]", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'k3'] <- ", mod,"_results2$par[3]", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'sd'] <- ", mod,"_results2$par[4]", sep="")))
        eval(parse(text=paste("output2['", mod, "', 'nll'] <- ", mod,"_results2$value", sep="")))
        
      } 
  }
}

results2 <- as.data.frame(matrix(NA, nrow = length(DesiredModels), ncol = 10))
row.names(results2) <- DesiredModels
results2[,1:3] <- round(output2[,1:3],4) # Rounded parameter estimates
results2[,4] <- round(output2[, 6],4) # Rounded sd estimate
results2[,5:7] <- output2[,7:9]  # Method, nll, BIC
results2[,10] <- conv2 # Convergence indicator
row.names(results2) <- DesiredModels
colnames(results2) <- c("k1", "k2", "k3", "sd","Method", "nll", "BIC","RMSE_Normalized","AdRSquare", "Converged?")

check_p <- array(NA, dim=c(length(DesiredModels), 1)) # Empty dataframe that will house the difference between second and third optimization parameter estimates (instability)
row.names(check_p) <- DesiredModels;
output_mid <- outputmid
for (m in 1:length(DesiredModels)) {
  mod <- DesiredModels[m]
  if (mod %in% c('ep', 'lg1')) {
    eval(parse(text=paste(mod, " <- pred_", mod, "(results2['", mod, "', 1], dat$Time, TRUE)", sep='')))  # Uses optimized model to predict LRV over time
    eval(parse(text=paste("results2['", mod, "', 'BIC'] <- round(2*results2['", mod, "', 6]+2*log(nrow(dat)), 2)", sep=''))) # Calculates BIC values
    r_ind <- 1 # Single k parameter (k1)
    eval(parse(text=paste("check_p['", mod, "',1] <- abs(output_mid['", mod, "',1]-output2['", mod, "',1])/output_mid['", mod, "',1]", sep='')))
  }
  if (mod %in% c('lg2','epd','jm1','jm2','gz','wb','ln','gam','bi')) {
    eval(parse(text=paste(mod, " <- pred_", mod, "(results2['", mod, "', 1], results2['", mod, "', 2], dat$Time, TRUE)", sep='')))
    eval(parse(text=paste("results2['", mod, "', 'BIC'] <- round(2*results2['", mod, "', 6]+3*log(nrow(dat)), 2)", sep='')))
    r_ind <- 2 # Model with two k parameters (k1, k2)
    eval(parse(text=paste("check_p['", mod, "',1] <- max(abs(output_mid['", mod, "',1]-output2['", mod, "',1])/output_mid['", mod, "',1], 
                          abs(output_mid['", mod, "',2]-output2['", mod, "',2])/output_mid['", mod, "',2])", sep='')))
  }
  if (mod %in% c('bi3', 'dep', 'gz3', 'gzm', 'sA', 'sB')) {
    eval(parse(text=paste(mod, " <- pred_", mod, "(results2['", mod, "', 1], results2['", mod, "', 2], results2['", mod, "', 3], dat$Time, TRUE)", sep='')))
    eval(parse(text=paste("results2['", mod, "', 'BIC'] <- round(2*results2['", mod, "', 6]+4*log(nrow(dat)), 2)", sep='')))
    r_ind <- 3 # Model with three k parameters (k1, k2, k3)
    eval(parse(text=paste("check_p['", mod, "',1] <- max(abs(output_mid['", mod, "',1]-output2['", mod, "',1])/output_mid['", mod, "',1], 
                          abs(output_mid['", mod, "',2]-output2['", mod, "',2])/output_mid['", mod, "',2], 
                          abs(output_mid['", mod, "',3]-output2['", mod, "',3])/output_mid['", mod, "',3])", sep='')))
  }
  
  eval(parse(text=paste("results2['", mod, "', 'RMSE_Normalized'] <- rmse(-dat$LR,", mod, ")/(max(dat$LR)-min(dat$LR))", sep='')))  # Calculates nRMSE values with observed and predicted LRV over time, normalized by range of LRV observed
  eval(parse(text=paste("results2['", mod, "', 'AdRSquare']  <- round((1-((1-(1-(sum((-dat$LR-", mod, ")^2)/sum((dat$LR-mean(dat$LR))^2))))*(n-1))/(n-r_ind-1)), 2)", sep=''))) # Calculates Adjusted R2 value
  
}

# Labels models with parameter estimates changing by greater than 1% after two optimization routines as unstable (>1% arbitrarily selected)
for (p in 1:nrow(check_p)) {
  if (check_p[p , 1]>0.01 | is.na(check_p[p, 1])) {
    results2[p, 'nll'] <- "Unstable"
    
  } 
}

## Remove models that have unstable parameter estimates, impossible nRMSE or AdR2 values, and that failed to converge
for (i in 1:nrow(results2)) {
  if (results2[i,"Method"]=="Unstable" | results2[i, "RMSE_Normalized"]>1 | results2[i, "AdRSquare"]< 0 | results2[i, "Converged?"] != 0 | 
      is.na(results2[i, "RMSE_Normalized"]) | is.na(results2[i, "AdRSquare"])) {
    results2[i, "BIC"] <- NA
    results2[i, "RMSE_Normalized"] <- NA
    results2[i, "AdRSquare"] <- NA
    results2[i, "Converged?"] <- NA
    
  } else {
    results2[i, "BIC"] <- results2[i, "BIC"]
    results2[i, "RMSE_Normalized"] <- results2[i, "RMSE_Normalized"]
    results2[i, "AdRSquare"] <- results2[i, "AdRSquare"]
    results2[i, "Converged?"] <- results2[i, "Converged?"]
  }
  
}

####Check for U-shaped EPD or DEP ####
## This will just indicate in the output that it was U-shaped, it will not prevent it from being selected as a best fitting model ##
if (which.min(pred_epd(results2["epd",1],results2["epd",2], dat$Time, TRUE))<length(dat$Time)) {
  
  results2['epd',5] <- "U-Shaped"
} else{
  results2['epd',5] <- results2['epd',5]
}

if (which.min(pred_dep(results2["dep",1],results2["dep",2], results["dep",3],dat$Time, TRUE))<length(dat$Time)) {
  results2["dep", 5] <- "U-Shaped"
} else{
  results2["dep",5] <- results2["dep",5]
  
}

# Function for calculating se
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

###Bootstrap the residuals####
minb <- results2$BIC[which.min(results2$BIC)]  # Model with lowest BIC
mod <- DesiredModels[which(abs(results2$BIC-minb)<=2)] # All models within 2 units of the lowest BIC
modb <- mod[which(results2[mod, "Converged?"]==0)]  # All models within 2 units of the lowest BIC that also converged
se_boot <- array(NA, dim=c(length(DesiredModels), 9))  # Empty dataframe for se of each parameter estimates
row.names(se_boot) <- DesiredModels

# The following for loop bootstraps the residuals for each model in modb (this part of the code takes the most time- future versions 
# should look at making this part of the code optional if a user wants to quickly look at optimized parameter values only)
for (i in 1:length(modb)) {
  mod <- modb[i]
  if (any(mod %in% c("ep", "lg1"))) {
    eval(parse(text=paste("boot_", mod, " <- pred_", mod, "(results2['",mod,"',1], dat$Time, TRUE)", sep="")))
    eval(parse(text=paste("res1 <- -dat$LR-boot_", mod, sep="")))
    eval(parse(text=paste("boot", mod, " <- array(NA, dim=c(1000,1))", sep="")))
    for (i in 1:1000) {
      bootid<-sample(1:n, n, replace = T)
      LRnew <- dat$LR+res1[bootid]
      LRnew[which(dat$Time==0)] <- 0
      datnew <- cbind(dat$Time, LRnew)
      colnames(datnew) <- c("Time", "LR")
      tryCatch(eval(parse(text=paste("optimnew <- optim(par=c(k1=results2['", mod,"','k1'], sd=results2['", mod,
                                     "','sd']), mod='", mod, "', fn=mod_L1, data=datnew, hessian=T)", sep=""))), 
               error=function(cond) {paste("optim for 1-parameter model with hessian=T failed")})
      
      if(exists("optimnew") & is.list(optimnew)) {
        eval(parse(text=paste("boot", mod, "[i, 1]=c(optimnew$par[1])", sep="")))
        
      } else {
        eval(parse(text=paste("optimnew <- optim(par=c(k1=results2['", mod,"','k1'], sd=results2['", mod,
                              "','sd']), mod='", mod, "', fn=mod_L1, data=datnew, hessian=F)", sep="")))
        
        
        if(exists("optimnew") & is.list(optimnew)) {
          eval(parse(text=paste("boot", mod, "[i, 1]=c(optimnew$par[1])", sep="")))
        } else {
          eval(parse(text=paste("boot", mod, "[i, 1]=NA", sep="")))
        }
      }
    }
    if (exists(paste("boot_", mod, sep=''))) {
      eval(parse(text=paste("se_boot['", mod, "',1:3] <- c(results2['", mod,"',1], quantile(boot", mod,"[,1], 0.025, na.rm=TRUE), quantile(boot", mod, "[,1], 0.975, na.rm=TRUE))", sep="")))
      
    } else { paste("se_boot['",mod, "', 1:3] <- c(NA, NA, NA)", sep="")}
    png(file.path(path, 'Outputs', fold, paste(datname, '_', mod, 'bootfig.png', sep="")), width=600, height=600)
    title=paste(mod, " & bootstrapped residuals", sep="")
    plot(dat$Time, -dat$LR, xlab="Time (Days)", ylab="log10(Nt/No)", main=title)
    eval(parse(text=paste("lines(dat$Time, pred_", mod, "(results2['", mod, "',1], dat$Time, TRUE), col='red')", sep="")))
    eval(parse(text=paste("lines(dat$Time, pred_", mod, "(quantile(boot", mod,", 0.025), dat$Time, TRUE), col='red', lty='dashed')", sep="")))
    eval(parse(text=paste("lines(dat$Time, pred_", mod, "(quantile(boot", mod,", 0.975), dat$Time, TRUE), col='red', lty='dashed')", sep="")))
    dev.off()
  }
  
  
  if (any(mod %in% c('lg2','epd','jm1','jm2','gz','wb','ln','gam','bi'))) {
    eval(parse(text=paste("boot_", mod, " <- pred_", mod, "(results2['",mod,"',1], results2['",mod,"',2], dat$Time, TRUE)", sep="")))
    eval(parse(text=paste("res1 <- -dat$LR-boot_", mod, sep="")))
    eval(parse(text=paste("boot", mod, " <- array(NA, dim=c(1000,2))", sep="")))
    for (i in 1:1000) {
      bootid<-sample(1:n, n, replace = T)
      LRnew <- dat$LR+res1[bootid]
      LRnew[which(dat$Time==0)] <- 0
      datnew <- cbind(dat$Time, LRnew)
      colnames(datnew) <- c("Time", "LR")
      eval(parse(text=paste("optimnew <- tryCatch(optim(par=c(k1=results2['", mod,"','k1'], k2=results2['", mod,"','k2'], 
                                    sd=results2['", mod,"','sd']), mod='", mod, "', 
                                    fn=mod_L2, data=datnew, hessian=T),
               error=function(cond) {paste('optim for 2-parameter model with hessian=T failed')})", sep="")))
      
      if(exists("optimnew") & is.list(optimnew)) {
        eval(parse(text=paste("boot", mod, "[i, 1:2]=c(optimnew$par[1], optimnew$par[2])", sep="")))
        
      } else {
        eval(parse(text=paste("optimnew <- tryCatch(optim(par=c(k1=results2['", mod,"','k1'], k2=results2['", mod,"','k2'],
                             sd=results2['", mod,"','sd']), mod='", mod, "', fn=mod_L2, data=datnew, hessian=F),
               error=function(cond) {paste('bootstrap failed for 2-parameter model')})", sep="")))
        
        if(exists("optimnew") & is.list(optimnew)) {
          eval(parse(text=paste("boot", mod, "[i, 1:2]=c(optimnew$par[1], optimnew$par[2])", sep="")))
        } else {
          eval(parse(text=paste("boot", mod, "[i, 1:2]=c(NA,NA)", sep="")))
        }
      }
    }
    if (exists(paste("boot_", mod, sep=''))) {
      eval(parse(text=paste("se_boot['", mod, "',1:6] <- c(results2['", mod,"',1], quantile(boot", mod,"[,1], 0.025, na.rm=TRUE), quantile(boot", mod, "[,1], 0.975, na.rm=TRUE),
                           results2['", mod,"',2], quantile(boot", mod,"[,2], 0.025, na.rm=TRUE), quantile(boot", mod, "[,2], 0.975, na.rm=TRUE))", sep="")))
      
    } else { paste("se_boot['",mod, "', 1:6] <- c(NA, NA, NA, NA, NA, NA)", sep="") }
    png(file.path(path, 'Outputs', fold, paste(datname, '_', mod, 'bootfig.png', sep="")), width=600, height=600)
    title=paste(mod, " & bootstrapped residuals", sep="")
    plot(dat$Time, -dat$LR, xlab="Time (Days)", ylab="log10(Nt/No)", main=title)
    eval(parse(text=paste("lines(dat$Time, pred_", mod, "(results2['", mod, "',1], results2['", mod, "',2], dat$Time, TRUE), col='red')", sep="")))
    for (i in 1:1000) {
      eval(parse(text=paste("lines(dat$Time, pred_", mod, "(boot", mod, "[i,1], boot", mod, "[i,2], dat$Time, TRUE), col='red', lty='dashed')", sep="")))
    }
    points(dat$Time, -dat$LR, cex=2, pch=16)
    dev.off()
  }
  
  
  
  if (any(mod %in% c('bi3', 'dep', 'gzm', 'gz3', 'sA', 'sB'))) {
    eval(parse(text=paste("boot_", mod, " <- pred_", mod, "(results2['",mod,"',1], results2['",mod,"',2],
                         results2['",mod,"',3], dat$Time, TRUE)", sep="")))
    eval(parse(text=paste("res1 <- -dat$LR-boot_", mod, sep="")))
    eval(parse(text=paste("boot", mod, " <- array(NA, dim=c(1000,3))", sep="")))
    for (i in 1:1000) {
      bootid<-sample(1:n, n, replace = T)
      LRnew <- dat$LR+res1[bootid]
      LRnew[which(dat$Time==0)] <- 0
      datnew <- cbind(dat$Time, LRnew)
      colnames(datnew) <- c("Time", "LR")
      eval(parse(text=paste("optimnew <- tryCatch(optim(par=c(k1=results2['", mod,"','k1'], k2=results2['", mod,"','k2'], 
                                    k3=results2['", mod,"','k3'], sd=results2['", mod,"','sd']), mod='", mod, "', 
                                    fn=mod_L3, data=datnew, hessian=T),
               error=function(cond) {paste('optim for 3-parameter model with hessian=T failed')})", sep="")))
      
      if(exists("optimnew") & is.list(optimnew)) {
        eval(parse(text=paste("boot", mod, "[i, 1:3]=c(optimnew$par[1], optimnew$par[2], optimnew$par[3])", sep="")))
        
      } else {
        eval(parse(text=paste("optimnew <- tryCatch(optim(par=c(k1=results2['", mod,"','k1'], k2=results2['", mod,"','k2'],
                             k3=results2['", mod,"','k3'], sd=results2['", mod,"','sd']), mod='", mod, "', 
                             fn=mod_L3, data=datnew, hessian=F),
               error=function(cond) {paste('bootstrap for 3-parameter model failed')})", sep="")))
        
        if(exists("optimnew") & is.list(optimnew)) {
          eval(parse(text=paste("boot", mod, "[i, 1:3]=c(optimnew$par[1], optimnew$par[2], optimnew$par[3])", sep="")))
        } else {
          eval(parse(text=paste("boot", mod, "[i, 1:3]=c(NA,NA, NA)", sep="")))
        }
      }
    }
    if (exists(paste("boot_", mod, sep=''))) {
      eval(parse(text=paste("se_boot['", mod, "', 1:9] <- c(results2['", mod,"',1], quantile(boot", mod,"[,1], 0.025, na.rm=TRUE), quantile(boot", mod, "[,1], 0.975, na.rm=TRUE),
                           results2['", mod,"',2], quantile(boot", mod,"[,2], 0.025, na.rm=TRUE), quantile(boot", mod, "[,2], 0.975, na.rm=TRUE),
                           results2['", mod,"',3], quantile(boot", mod,"[,3], 0.025, na.rm=TRUE), quantile(boot", mod, "[,3], 0.975, na.rm=TRUE))", sep="")))
      
    } else { paste("se_boot['",mod, "', 1:9] <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA)", sep="") }
    png(file.path(path, 'Outputs', fold, paste(datname, '_', mod, 'bootfig.png', sep="")), width=600, height=600)
    title=paste(mod, " & bootstrapped residuals", sep="")
    plot(dat$Time, -dat$LR, xlab="Time (Days)", ylab="log10(Nt/No)", main=title)
    eval(parse(text=paste("lines(dat$Time, pred_", mod, "(results2['", mod, "',1], results2['", mod, "',2], results2['", mod, "',3],
                         dat$Time, TRUE), col='red')", sep="")))
    for (i in 1:1000) {
      eval(parse(text=paste("lines(dat$Time, pred_", mod, "(boot", mod, "[i,1], boot", mod, "[i,2], boot", mod, "[i,3],
                           dat$Time, TRUE), col='red', lty='dashed')", sep="")))
    }
    points(dat$Time, -dat$LR, cex=2, pch=16)
    dev.off()
  }
}


row.names(se_boot) <- DesiredModels
colnames(se_boot) <- c("k1", "k1 2.5%", "k1 97.5%", "k2", "k2 2.5%", "k2 97.5%","k3", "k3 2.5%", "k3 97.5%")

se_boot_subset <- subset(se_boot, rownames(se_boot) %in% modb)

# Plot of all best fitting models together #
leg_lab <- 0; leg_col <- 0; leg_lty <- 0
png(file.path(path, 'Outputs', fold, paste(datname, '_Summary.png', sep="")), width=600, height=600)
plot(dat$Time, -dat$LR, main=paste(datname, "Best Fitting Models", sep=" "), xlab="Time (Days)", ylab="LRV")
if (max(dat$Time) < 1) {
  tt <- seq(0, max(dat$Time), 0.01)
} else {
  tt <- seq(0, max(dat$Time), 0.1)
}
for (plot in 1:length(modb)){
  plot_mod <- modb[plot]
  if (plot_mod %in% c('ep', 'lg1')) {
    eval(parse(text=paste("lines(tt, pred_", plot_mod, "(results2['", plot_mod, "',1],tt, TRUE), col=GSs['", plot_mod, "', 'color'], 
                          lty=GSs[", plot_mod, ", 'lty'])", sep="")))
  }
  if (plot_mod %in% c('lg2','epd','jm1','jm2','gz','wb','ln','gam','bi')) {
    eval(parse(text=paste("lines(tt, pred_", plot_mod, "(results2['", plot_mod, "',1], results2['", plot_mod, "',2],
                          tt, TRUE), col=GSs['", plot_mod, "', 'color'], lty=GSs[", plot_mod, ", 'lty'])", sep="")))
  }
  if (plot_mod %in% c('bi3', 'dep', 'gzm', 'gz3', 'sA', 'sB')) {
    eval(parse(text=paste("lines(tt, pred_", plot_mod, "(results2['", plot_mod, "',1], results2['", plot_mod, "',2],
    results2['", plot_mod, "',3], tt, TRUE), col=GSs['", plot_mod, "', 'color'], lty=GSs[", plot_mod, ", 'lty'])", sep="")))
  }
  leg_lab[plot] <- eval(parse(text=paste("paste(GSs['",plot_mod,"','model'],': ',results2['",plot_mod,"','BIC'], sep='')",sep='' )))
  leg_col[plot] <-  eval(parse(text=paste("GSs['", plot_mod, "', 'color']", sep="")))
  leg_lty[plot] <-  eval(parse(text=paste("GSs['", plot_mod, "', 'lty']", sep="")))
}
legend("topright", legend=leg_lab, col=leg_col, lty=leg_lty)
dev.off()

###Hessian Check####
## The mle function does produce Fischer matrix that can be transformed to calculate parameter se
se_tot2 <- array(NA, dim=c(length(DesiredModels), 4))
rownames(se_tot2) <- DesiredModels
for (j in 1:length(DesiredModels)) {
  mod <- DesiredModels[j]
  eval(parse(text=paste("hess <- ", mod, "_results2$hessian", sep="")))
  if (is.null(hess)) {
    se=NA
  } else{
    se <- tryCatch(sqrt(diag(solve(hess))), error=function(cond) {paste(NA)})
  }
  
  if (mod %in% c('ep', 'lg1')) {
    se_tot2[mod, 1] <- round(as.numeric(se[1]),3); se_tot2[mod,4] <- round(as.numeric(se[2]),3);
  }
  if (mod %in% c('lg2','epd','jm1','jm2','gz','wb','ln','gam','bi')) {
    se_tot2[mod, 1:2] <- round(as.numeric(se[1:2]),3); se_tot2[mod,4] <- round(as.numeric(se[3]),3);
  }
  if (mod %in% c('bi3', 'dep', 'gzm', 'gz3', 'sA', 'sB')) {
    se_tot2[mod, 1:4] <- round(as.numeric(se),3);
  }
}

colnames(se_tot2) <- c("k1 se", "k2 se", "k3 se", "sd se")
se_tot_subset <- subset(se_tot2, rownames(se_tot2) %in% modb)


write.csv(results2, file.path(path, 'Outputs', fold, paste(datname, '_Output2', '.csv', sep="")))
write.csv(se_boot, file.path(path, 'Outputs', fold, paste(datname, '_bootCI', '.csv', sep="")))

####Output- Create a succinct dataframe for quick comparisons ####
bfparam <- results2[which(rownames(results2) %in% modb),]
param <- cbind.data.frame(rep(datname, nrow(bfparam)), rownames(bfparam), bfparam)
colnames(param) <- c("Dataset", "Model", "k1", "k2", "k3", "sd","Method", "nll", "BIC","RMSE_Normalized","AdRSquare", "Converged?")
param_se <- cbind.data.frame(rep(datname, nrow(se_boot_subset)), rownames(se_boot_subset), round(se_boot_subset,2))
colnames(param_se) <- c("Dataset", "Model", "k1", "k1 2.5%", "k1 97.5%", "k2", "k2 2.5%", "k2 97.5%","k3", "k3 2.5%", "k3 97.5%")
outputinfo <- list("optim2"=param, "boot"=param_se)

outputinfo
