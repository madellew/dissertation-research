
#QMRA model for Candida auris catheter-associated urinary tract infections

#Created: 9 Dec 2023
#Last modified: 21 Jan 2025

#Calculate the probability of infection per HCW patient visit
#Estimates for three measures of environmental contamination: (1) start patient shed low, (2) start patient shed high, (3) existing contam
#Can factor in disinfection log-reduction (yes/no), whether the contact is guaranteed to be with a dirty fomite (yes/no), 
#and if a guaranteed break in sterility occurred during HCW-catheter contact (yes/no)

#Require necessary packages
library(triangle)
library(stats)
library(ggplot2)
library(reshape2)
library(paletteer)


set.seed(1234)
runtime = 'full' # 'test'
if(runtime == 'full'){iters = 10000}else{iters = 1000}

#Read CSV for bootstrapped dose-response parameters
k_bp_cauti <- read.csv("exp_bp_AR0384.csv")

#####
#Parameterize environmental contamination

contam <- c("Shed Low", "Shed High", "All")
N.contam <- length(contam)

p_col_conc1 <- rnorm(iters,1.983,1.765) #low shed
p_col_conc2 <- runif(iters,-1.183,5.957) #high shed
shed_effc <- 0.536
shed_freq <- 0.5
t_shed <-runif(iters,0,1)
f_dirty_conc_shed1 <- ((10^p_col_conc1*shed_effc)/shed_freq)*t_shed
f_dirty_conc_shed2 <- ((10^p_col_conc2*shed_effc)/shed_freq)*t_shed

f_dirty_conc_all <- rnorm(iters, 0.9886, 1.063)
f_dirty_conc_all <- 10^f_dirty_conc_all

contam.all <- data.frame(f_dirty_conc_shed1, f_dirty_conc_shed2, f_dirty_conc_all)
colnames(contam.all) <- contam

summary(contam.all)
summary(log10(contam.all))
write.csv(contam.all,file = "Env_Conc.csv")

#####
#Estimate fomite concentration with or without surface disinfection

dis <- vector() #cleaning efficacy

f_conc <- matrix(0, iters, ncol = N.contam)  # fomite concentration 
colnames(f_conc) <- contam

#yes, include disinfection; no, no disinfection
is.dis = 'no'
if(is.dis == 'yes'){dis.yes = 1}else{dis.yes = 0}

for(i in 1:iters)
{
  for(j in 1:ncol(contam.all))
  {
    dis[i] <- runif(1,0,7)*dis.yes
    f_conc[i,j] <- contam.all[i,j]*10^(-dis[i])
  }
}

f_conc <- as.data.frame(f_conc)

rm(j,i)  

#####
#Estimate concentration on HCW hand following pick up from fomite

f_h_N <- vector() #number of fomite-hand contacts during patient visit
f_h_effc <- vector() #fomite-hand transfer efficiency 
f_decay <- vector() #decay rate - fomite
f_t <- vector() #time between fomite was contaminated to hand contact

#Juneja and Marks 2
decay_jm2 = function(K1,K2,t){	
  1/(1+exp(K1+K2*log(t)))
}

k1 <- -4.1147
k2 <- 1.6605

h_conc <- matrix(0, iters, ncol = N.contam)  # fingertip concentration
colnames(h_conc) <- contam

f_h_N <- rnbinom(iters,1.04,0.229)

is.dirt = 'yes' #yes, assume all dirty; no, probability contacted fomite is dirty
if(is.dirt == 'yes'){dirt.yes = f_h_N+1}else{dirt.yes = rbinom(iters,f_h_N,0.322)} 


for(i in 1:iters)
{
  for(j in 1:ncol(contam.all))
  {
    f_h_effc[i] <- rtriangle(1,0.00007143,0.007143,0.001757) #min, max, mode
    #f_h_effc[i] <- runif(1, 0.396,0.619) #high RH (sens)
    f_t[i] <- runif(1,0,1)
    f_decay[i] <- decay_jm2(k1,k2,f_t[i])
    #f_decay[i] <- decay_epd(c,d,f_t[i]) #plastic decay (sens)
    h_conc[i,j] <- dirt.yes[i]*f_conc[i,j]*f_h_effc[i]*f_decay[i]
  }
}

h_conc <- as.data.frame(h_conc)

summary(h_conc)
summary(log10(h_conc))
write.csv(h_conc,file = "Hand_Conc.csv")

rm(j,i)

#####
#Estimate concentration transferred to urinary catheter (dose)

h_area <- vector() #fingertip surface area
h_c_N <- vector() #hand-catheter transfer frequency
h_c_effc <- vector() #hand-catheter transfer efficiency
h_decay <- vector() #decay rate - hands
h_t <- vector() #time between fomite contacted and catheter contact

dose <- matrix(0, iters, ncol = N.contam) #Dose
colnames(dose) <- contam

is.sterile = 'no' #yes, probability break in sterile barriers; no, assume all contacts
if(is.sterile == 'yes'){sterile.yes = rbinom(iters,1,0.59)}else{sterile.yes = runif(iters,1,1)} 


for(i in 1:iters)
{
  for(j in 1:ncol(contam.all))
  {
    h_area[i] <- 10
    h_c_N[i] <- sample(c(1,2), 1, replace = TRUE)
    h_c_effc[i] <- rtriangle(1,0.00007143,0.06429,0.008131)
    #h_c_effc[i] <- rtriangle(1,0,0.3447,0.1026) #(sens)
    h_t[i] <- rtriangle(1,3.33,120,13.33)
    h_decay[i] <- exp(-0.07296*h_t[i])
    dose[i,j] <- sterile.yes[i]*h_conc[i,j]*h_area[i]*h_c_N[i]*h_c_effc[i]*h_decay[i]
  }
}

dose <- as.data.frame(dose)

summary(dose)
summary(log10(dose))
write.csv(dose,file = "Dose.csv")

rm(j,i)

#####
#Estimate risk of infection (CAUTI) using exponential dose-response model

Risk.contam <- matrix(0, iters, ncol = N.contam)  #risk of infection
colnames(Risk.contam) <- contam

for(i in 1:iters)
{
  for(j in 1:ncol(contam.all))
  {
  k <- k_bp_cauti$k
  Risk.contam[i,j] <- 1-exp(-k[i]*dose[i,j]) #per HCW visit to patient
  }
}

Risk.contam <- as.data.frame(Risk.contam)

summary(Risk.contam)
summary(log10(Risk.contam))
write.csv(Risk.contam,file = "Risk.csv")
colSums(Risk.contam==0) #frequency of "zero" risk

rm(j,i)


#### Summary statistics and plotting ####

#Summary statistics 

CAUTI.descrips <- matrix(0,7,3)

for(i in 1:ncol(Risk.contam))
{
  CAUTI.descrips[1,i] <- min(Risk.contam[,i])
  CAUTI.descrips[2,i] <- quantile(Risk.contam[,i], probs=0.05)
  CAUTI.descrips[3,i] <- median(Risk.contam[,i])
  CAUTI.descrips[4,i] <- quantile(Risk.contam[,i], probs=0.95)
  CAUTI.descrips[5,i] <- mean(Risk.contam[,i])
  CAUTI.descrips[6,i] <- sd(Risk.contam[,i])
  CAUTI.descrips[7,i] <- max(Risk.contam[,i])
}
colnames(CAUTI.descrips) <- colnames(Risk.contam)
rownames(CAUTI.descrips) <- c("Minimum","Lower 95th","Median","Upper 95th","Mean","Standard Deviation","Maximum")
write.csv(CAUTI.descrips,file = "CAUTI_Risk_Descriptive_Statistics.csv")

#Distribution of risk estimates 

Risk.contam.plot <- as.matrix(Risk.contam) 
Risk.contam.ggplot <- melt(Risk.contam.plot)
colnames(Risk.contam.ggplot) <- c("dummy","Source","Risk")

#Histogram
CAUTI.hist <- ggplot(Risk.contam.ggplot, aes(x=Risk, fill=Source, color=Source)) + 
  geom_histogram(position = "identity", alpha=0.5, bins=40) +
  scale_x_log10() +
  labs(x = c(expression(~log[10]~Probability~of~Infection)),
       y = "Frequency", fill="Contamination Source", color="Contamination Source") +  
  theme_bw() +
  theme(text=element_text(family="serif",size=13)) +
  scale_color_paletteer_d("rcartocolor::Bold", labels=c("Low Shed", "High Shed", "Pre-seeded")) +
  scale_fill_paletteer_d("rcartocolor::Bold", labels=c("Low Shed", "High Shed", "Pre-seeded"))
ggsave("CAUTI_Risk_hist.png",CAUTI.hist,dpi=600)

#Box/Violin plot
CAUTI.box <- ggplot(Risk.contam.ggplot, aes(Source, Risk, fill=Source)) + 
  #geom_boxplot(alpha=0.75) + 
  geom_violin(alpha=0.75, trim=FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(text=element_text(family="serif",size=13), axis.text.x=element_text(size = 12)) + 
  labs(x = "Environmental Contamination Source", 
       y = expression(~log[10]~Probability~of~Infection),fill="Source") + 
  #stat_summary(fun=median, geom="point", size=2.5, shape=5) + 
  #stat_summary(fun=mean, geom="point", size=2, color="black",shape=16) + 
  stat_summary(fun.data="mean_sdl", fun.args = list(mult = 1), geom="pointrange") +
  scale_fill_paletteer_d("rcartocolor::Bold")
ggsave("CAUTI_Risk_box.png",CAUTI.box)

#Jitterplot
CAUTI.jitter <- ggplot(Risk.contam.ggplot, aes(Source,Risk)) +
  geom_jitter(aes(colour=Source)) +
  theme_bw() +
  scale_y_log10() +
  labs(x = "Environmental Contamination Source",
       y = c(expression(~log[10]~Probability~of~Infection))) 
ggsave("CAUTI_Risk_jitter.png",CAUTI.jitter)

#####
#Additional R scripts for Chapter 7 QMRA
source('cauti_qmra_sens.R') #CAUTI sensitivity analysis
source('clabsi_qmra.R') #CLABSI risk estimation
source('clabsi_qmra_sens.R') #CLABSI sensitivity analysis
source('cath_qmra_growth.R') #Risk estimation with pathogen growth on/in catheter
source('cath_qmra_cum.R') #Cumulative risk estimation

#End

