#Reading in data#

library(dplyr)
library(readr)

files <- list.files(path = "./Data_Cdiff",
                    pattern = "*.csv", full.names = T)  # Pulls alls the csv files from the Data folder
tbl <- sapply(files, read_csv, simplify=FALSE) %>% 
  bind_rows(.id = "id")       # Reads in the data from each csv file and creates a giant table

data = tbl[,1:16] #Limits the columns being considered to those 16 shown in the template 

data$Species <-  as.character(data$Species)   
data$Fomite = as.character(data$Fomite)
data$Experiment = mapply(paste,data$Species,data$Fomite, sep="") # Experiment name is the Species + Fomite
data$Fomite = as.factor(data$Fomite)
data$Species = as.factor(data$Species)
data$Source= as.character(data$Source)

table(data$Species,data$Fomite)	


dLTm = aggregate(Count ~ Species  + Fomite + log10 + Time + Conc + Experiment+ Source, data=data, FUN=mean) # Core information needed from each dataset
dLTmZeros = dLTm
dLTm$Count[which(dLTm$Count==0)] = 10  # Opportunity to assign a limit of detection to observations of 0
dLTm$Time[which(dLTm$Time < 0)] = 0 # Remove any negative time points from digitizing errors
dLTm$pLeft = dLTm$Count / (dLTm$Conc)  
dLTm$LR = -log10(dLTm$pLeft)	# Calculate log10 reduction.

exps <- unique(data$Experiment) # Number of experiments based on unique combinations of Species + Fomite

for (i in 1:length(exps)) {
  datasetName =  paste(exps[i], "lt", sep="")  # Adding "lt" to each experiment name
  assign(datasetName,  dLTm %>% dplyr::filter(Experiment== exps[i]))
}


