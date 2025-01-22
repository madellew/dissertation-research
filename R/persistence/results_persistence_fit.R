## Script version of the R Markdown Results file (RMarkdown_Results_Persistence_Fit)
## Running the full script will have the same outputs of the R Markdown file but 
## will allow the user to work with the intermediate variables in the environment. 

##Created by Jade Mitchell, PhD and Ombaka-Owade, PhD; Michigan State University
##Modified on 20 July 2024 by Madeline Lewis - madellew@umich.edu

# Needed packages: 
library(tidyr)
library(dplyr)
library(Metrics)
library(readr)
library(knitr)


# Read in the data: 
source("DataRead_Folder.R") # Calls the function that reads all the datasets in the "Data" folder into the environment

alldatasets <- names(mget(ls(pattern="lt", all.names = TRUE))) # Stores the names of each dataset in the variable "alldatasets"
length(alldatasets) # The number of datasets being analyzed

# Confirm datasets have at least 4 observations:
min_obs=function(dataset){
  eval(parse(text=paste("dat=",dataset,"[,c('Time','LR')]",sep='')))
  m <- length(dat$Time)
}
min_obs_check <- tryCatch(lapply(alldatasets, min_obs), 
                          error=function(cond) {paste("One dataset could not be checked")})
names(min_obs_check) <- alldatasets
datasets_4up <- names(min_obs_check[which(min_obs_check>3)]) # New list of dataset names with 4 or more observations

# Evaluate for signficant trend
toft <- array(NA, dim=c(length(datasets_4up), 2)); 
for (i in 1:length(datasets_4up)) {
  dataset =  datasets_4up[i] #The name of the dataset that you wish to analyze.
  toft[i,1] <- dataset
  eval(parse(text=paste("dat=",dataset,"[,c('Time','LR')]",sep='')))
  trend = lm(-LR ~ 0 + Time, data=dat)	#Forcing intercept through 0
  summary(trend)
  Slope = summary(trend)$coefficients[1,1]
  pSlope = summary(trend)$coefficients[1,4]
  if (sign(Slope)==-1 & pSlope < 0.05){
    toft[i,2] <- "N" # significant negative trend
  } 
  if (sign(Slope)==1 & pSlope < 0.05){
    toft[i,2] <- "P" # significant positive trend
  } 
  if (pSlope > 0.05){
    toft[i,2] <- "NA" # no significant trend
  }
}
row.names(toft) <- datasets_4up
notrend <- subset(toft, toft[,2]=="NA")[,1]
n_trend <- subset(toft, toft[,2]=="N")[,1]
p_trend <- subset(toft, toft[,2]=="P")[,1]

write.csv(toft, "toft.csv")

# Fit the models to all of the datasets with a significant negative trend at once: 
source('persistfit_17_function.R') # Calls the function that fits 17 models to each dataset
start_time <- Sys.time()  # Notes the starting time 
fulloutput <- tryCatch(lapply(p_trend, persistfit_17_function),  #specify trend
                       error=function(cond) {paste("One dataset could not be optimized")})  # Fits the 17 models to each dataset in "alldatasets"
end_time <- Sys.time() # Notes the ending time
end_time-start_time  # The time required to run the number of datasets in "alldatasets"

optimparam <- data.table::rbindlist(sapply(fulloutput, '[', 1))  # Collates the best fitting models for all the datasets in one dataframe
optimparam
write.csv(optimparam, "optimparam.csv")

optimboot <- data.table::rbindlist(sapply(fulloutput, '[', 2))  # Collates the bootstrapped parameter estimates for each dataset's best fitting model in one dataframe
optimboot
write.csv(optimboot, "optimboot.csv")

# Fit the models to one dataset at a time: 
persistfit_17_function(n_trend)  # Calls the function for the first dataset in n_trend


# Use model averaging techniques to calculate metrics of interest: 
source('T90_T9999_metric_calc.R') # Calls the function that calculates T90s, T99s, T999s, or T9999s for each dataset using the best fitting models (model averaging)

param4pred <- cbind.data.frame(optimparam[,1:5], optimparam[,"BIC"]) # Creates a dataframe composed of dataset name, model name, k1, k2, k3, and BIC
colnames(param4pred) <- c("Dataset", "Model", "k1", "k2", "k3", "BIC")

T90pred <- metcalc(param4pred, 'T90', 50) # Calculates the T90 for each dataset, with an assumed maximum time of 50 days (function searches for the T90 between 0 and 30 days)
T90pred
write.csv(T90pred$metric, "T90pred_metric.csv")
write.csv(T90pred$MA_metric, "T90pred_MAmetric.csv")

T99pred <- metcalc(param4pred, 'T99', 100) # Calculates the T99 for each dataset, with an assumed maximum time of 100 days (function searches for the T99 between 0 and 60 days)
T99pred
write.csv(T99pred$metric, "T99pred_metric.csv")
write.csv(T99pred$MA_metric, "T99pred_MAmetric.csv")

T999pred <- metcalc(param4pred, 'T999', 150) # Calculates the T90 for each dataset, with an assumed maximum time of 150 days (function searches for the T90 between 0 and 30 days)
T999pred
write.csv(T999pred$metric, "T999pred_metric.csv")
write.csv(T999pred$MA_metric, "T999pred_MAmetric.csv")

T9999pred <- metcalc(param4pred, 'T9999', 200) # Calculates the T99 for each dataset, with an assumed maximum time of 200 days (function searches for the T99 between 0 and 60 days)
T9999pred
write.csv(T9999pred$metric, "T9999pred_metric.csv")
write.csv(T9999pred$MA_metric, "T9999pred_MAmetric.csv")

T99999pred <- metcalc(param4pred, 'T99999', 300) # Calculates the T99 for each dataset, with an assumed maximum time of 300 days (function searches for the T99 between 0 and 60 days)
T99999pred
write.csv(T99999pred$metric, "T99999pred_metric.csv")
write.csv(T99999pred$MA_metric, "T99999pred_MAmetric.csv")
