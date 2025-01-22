#Load HCUP data and prepare for data analysis

#Created: 12 Nov 2024
#Last modified: 13 Nov 2024

#Load required packages
library(tidyverse)
library(psych)
library(Hmisc)
library(gtsummary)
library(hablar)
library(dplyr)
library(tidycomm)
library(reshape2)
library(fitdistrplus)
library(VGAM)
library(mixtools)
library(mclust)
library(paletteer)
library(ggeffects)

#Import data
hcup <- read.csv("hcup_new.csv", header=TRUE)
hcup_hosp <- read.csv("hcup_hosp.csv", header=TRUE)

#Let's look at hospital data before joining
hcup_hosp <- hcup_hosp %>%
  convert(fct(HOSP_BEDSIZE, HOSP_LOCTEACH, HOSP_REGION))

hcup_hosp %>%
  dplyr::select(HOSP_BEDSIZE, HOSP_LOCTEACH, HOSP_REGION, TOTAL_DISC) %>%
  Hmisc::describe()

#Joining datasets
data_hcup <- left_join(hcup, hcup_hosp, by = "HOSP_NIS")

#Meaningful categories
data <- as_tibble(data_hcup)
data <- data %>%
  mutate(DISPUNIFORM=ifelse(DISPUNIFORM==99, NA, DISPUNIFORM)) %>%
  mutate(HCUP_ED=if_else(HCUP_ED>0,1,0)) %>%
  mutate(URCATH1=if_else(URCATH>0,1,0)) %>%
  mutate(CVC1=if_else(CVC>0,1,0)) %>%
  mutate(CAN=if_else(CAN>0,1,0))

#Create variables for outcomes when LOS is greater than or equal to 2 days
data <- data %>%
  mutate(CDI_HA=if_else(CDI_np==1 & LOS>=2,1,0)) %>%
  mutate(CAN_HA=if_else(CAN_np==1 & LOS>=2,1,0)) %>%
  mutate(CLABSI_HA=if_else(CLABSI_np==1 & LOS>=2,1,0)) %>%
  mutate(CAUTI_HA=if_else(CAUTI_np==1 & LOS>=2,1,0))

#Convert categorical variables to factor variable
str(data) #check variable type
data <- data %>%
  convert(fct(AMONTH, DIED, DISPUNIFORM, DQTR, ELECTIVE, FEMALE, HCUP_ED, 
              RACE, TRAN_IN, TRAN_OUT, CDI, CDI1, CDI_np, CDI_HA, CAN_UTI,
              CAN_UTI_np, CAN, CAN_HA, CAN_np, CLABSI, CLABSI_HA, CLABSI_np,
              CAUTI, CAUTI_HA, CAUTI_np, HAI, ABX, MRSA, URCATH, URCATH1, 
              CVC, CVC1)) %>%
  mutate(LOS.num=as.numeric(LOS))
  
#Univariate summary statistics
Hmisc::describe(data$AGE)
Hmisc::describe(data$LOS)


data %>%
  dplyr::select(AMONTH, DIED, DISPUNIFORM, DQTR, ELECTIVE, FEMALE, HCUP_ED, 
                RACE, TRAN_IN, TRAN_OUT, CDI, CDI1, CDI_np, CDI_HA, CAN_UTI,
                CAN_UTI_np, CAN, CAN_HA, CAN_np, CLABSI, CLABSI_HA, CLABSI_np,
                CAUTI, CAUTI_HA, CAUTI_np, HAI, ABX, MRSA, URCATH, URCATH1, 
                CVC, CVC1) %>%
  Hmisc::describe()

#Devices and related infections
tabyl(data, CLABSI_HA, CVC1)
stats::chisq.test(data$CLABSI_HA, data$CVC1)
tabyl(data, CAUTI_HA, URCATH1)
stats::chisq.test(data$CAUTI_HA, data$URCATH1)

#Candidiasis
tabyl(data, CAUTI_HA, CAN_HA)
stats::chisq.test(data$CVC1, data$CAN_HA)

