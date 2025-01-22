
#QMRA model for Clostridioides difficile infection (CDI)

#Created: 2 May 2022
#Last modified: 21 Jan 2025

#Calculate the probability of infection per exposure - single hand-to-surface contact, followed by single hand-to-mouth contact (dose).
#Or extend to estiamte daily infection risks.
#Estimates for four measures of environmental contamination: (1) colonized patient shedding, (2) infected patient shedding, 
#(3) existing contamination (pre-seeded), and (4) HCW pick up from direct patient care.
#Can factor in disinfection log-reduction (yes/no), whether the contact is always with a dirty fomite (yes/no), whether contacted patient
#is under contact precautions (yes/no), and hand hygiene (yes/no).

#Require necessary packages
library(EnvStats)
library(mc2d)
library(gsl)
library(stats)
library(ggplot2)
library(reshape2)


set.seed(1234)
runtime = 'full' # 'full'
if(runtime == 'full'){iters = 10000}else{iters = 1000}

#####
#Vectors and constants

contam <- c("Shed-Col", "Shed-Inf", "Env") #environmental contamination where "shed-col" and "shed-inf" start with an infected or 
                                            #colonized patient shedding onto a fomite; or a fomite is pre-seeded ("env")
contam2 <- c("Shed-Col", "Shed-Inf", "Env", "P-HCW") #where "p-hcw" is HCW hand contamination from direct patient contact 

N.contam <- length(contam)
N.contam2 <- length(contam2)
HealthEffects <- c("Infection","Illness")

#Concentrations
p_conc_col <- vector() #concentration of colonized patient
p_conc_inf <- vector() #concentration of infected patient

#Concentration on "dirty" fomite
f_dirty_conc_col <- vector()  #concentration on fomite shed from colonized patient
f_dirty_conc_inf <- vector()  #concentration on fomite shed from infected patient
f_dirty_conc_env <- vector() #concentration on pre-seeded fomite

f_conc <- matrix(0, iters, ncol = N.contam) #concentration on fomite after disinfection, if any

#HCW hand concentration
hcw_p_conc <- vector() #concentration on HCW hands after direct patient contact
hcw_f_conc <- matrix(0, iters, ncol=N.contam) #concentration on HCW hand after fomite contact
hcw_conc <- matrix (0, iters, ncol=N.contam2) # concentration on HCW after hand hygiene, if any

f_conc_cc <-  matrix (0, iters, ncol=N.contam2) #concentration on fomite after cross-contamination

h_conc <- matrix(0, iters, ncol = N.contam2)  #concentration on patient hand
dose <- matrix(0, iters, ncol = N.contam2) #Dose

colnames(f_conc) <- contam
colnames(hcw_f_conc) <- contam
colnames(hcw_conc) <- contam2
colnames(f_conc_cc) <- contam2
colnames(h_conc) <- contam2
colnames(dose) <- contam2

#Times
f_t <- vector() #time between fomite was contaminated/disinfected to hand contact

#Transfer efficiencies and number of contact events
shed_effc <- 4.85E-07
hcw_p_pick <- vector() #pick-up rate
hcw_p_N <- vector() #number of HCW contacts with patient
hcw_f_N <- vector() #number of HCW contacts with fomite
hcw_f_effc <- vector() #hand-to-fomite transfer efficiency
f_h_N <- vector() #number of patient contacts with fomites
f_h_effc <- vector() #fomite-hand transfer efficiency 
h_m_N <- vector() #number of hand-to-mouth contacts
h_m_effc <- vector() #hand-mouth transfer efficiency
h_area <- vector() #fingertip surface area

#Disinfection, growth and decay
dis <- vector() #cleaning efficacy
hh <- vector() #hand hygiene
f_decay <- vector() #decay rate from fomite

#Dose-response and risk characterization
p.inf.exact <- function(dose, a, b) 
  1 - hyperg_1F1(a,a+b,-dose)

alpha.inf <- 0.07403
beta.inf <- 1.120

alpha.ill <- 0.01202
beta.ill <- 0.03711

Risk.inf <- matrix(0, iters, ncol = N.contam2) #probability of infection
Risk.ill <- matrix(0, iters, ncol = N.contam2) #probability of illness

colnames(Risk.inf) <- c(contam2)
colnames(Risk.ill) <- c(contam2)


#####
#Exposure model

#Model patient shedding onto fomite

p_conc_col <- rnorm(iters,3.6,1.3)
p_conc_inf <- rnorm(iters,5.6,1.4)
f_dirty_conc_col <- (10^p_conc_col)*shed_effc
f_dirty_conc_inf <- (10^p_conc_inf)*shed_effc
f_dirty_conc_env <- rnormTrunc(iters, 11, 14, min=0)   
f_dirty_conc_env <- f_dirty_conc_env/100

contam.env <- data.frame(f_dirty_conc_col, f_dirty_conc_inf, f_dirty_conc_env)
colnames(contam.env) <- contam

summary(contam.env)
summary(log10(contam.env))
#write.csv(contam.env,file = "Env_Conc.csv")

#####
#Model disinfection of fomite

is.dis = 'no' #yes, include disinfection; no, no disinfection
if(is.dis == 'yes'){dis.yes = 1}else{dis.yes = 0}

for(i in 1:iters)
{
  for(j in 1:ncol(contam.env))
  {
    dis[i] <- runif(1,0,6)*dis.yes
    f_conc[i,j] <- contam.env[i,j]*10^(-dis[i])
  }
}

f_conc <- as.data.frame(f_conc)

rm(j,i)

#####
#Model HCW pick-up from fomite

f_h_N <- rnbinom(iters,1.83,0.167)

decay_epd = function(K1,K2,t){
  out = exp(-K1*t*exp(-K2*t))
}

k1 <- 0.117
k2 <- 0.0257

#k1 <- -10.09 #growth (sens)
#k2 <- 0.9572 #growth (sens)

is.dirt = 'yes' #yes, assume all dirty; no, probability contacted fomite is dirty
if(is.dirt == 'yes'){dirt.yes = f_h_N+1}else{dirt.yes = rbinom(iters,f_h_N,0.5)} 


for(i in 1:iters)
{
  for(j in 1:ncol(contam.env))
  {
    f_h_effc[i] <- rtriang(1,min=0.475,max=0.714,mean=0.57) 
    f_t[i] <- runif(1,0,1)
    f_decay[i] <- decay_epd(k1,k2,f_t[i])
    hcw_f_conc[i,j] <- dirt.yes[i]*f_conc[i,j]*f_h_effc[i]*f_decay[i]
  }
}

rm(j,i)

hcw_f_conc <- as.data.frame(hcw_f_conc)
summary(hcw_f_conc)
summary(log10(hcw_f_conc))
colSums(hcw_f_conc==0)

#####
#Model HCW pick-up from patient

is.isolation = 'no' #yes, probability of contacting asymptomatic colonized patient or efficiency of contact precautions; 
#no, assume all contacts
if(is.isolation == 'yes'){isolation.yes = rbinom(iters,1,0.60)}else{isolation.yes = runif(iters,1,1)} 
#if(is.isolation == 'yes'){isolation.yes = rbinom(iters,1,0.33)}else{isolation.yes = runif(iters,1,1)} #corrected for effectiveness

for(i in 1:iters)
{
    hcw_p_pick[i] <- rtriang(1,min=0.00204,max=0.408,mode=0.0286)
    #hcw_p_pick[i] <- rnormTrunc(1, 0.247, 0.0429, min=0) #alt pick up rate (sens)
    hcw_p_N[i] <- rnbinom(1,3.47,0.542)+1
    hcw_p_conc[i] <- isolation.yes[i]*hcw_p_pick[i]*hcw_p_N[i]
}

contam.all <- data.frame(hcw_f_conc, hcw_p_conc)
colnames(contam.all) <- contam2

summary(contam.all)
summary(log10(contam.all))
colSums(contam.all==0)
#write.csv(contam.all,file = "HCW_Conc.csv")                

#####

#Model hand hygiene

is.hh = 'no' #yes, include hand hygiene; no, no hand hygiene
if(is.hh == 'yes'){hh.yes = rbinom(iters,1,0.60)}else{hh.yes = 0}

for(i in 1:iters)
{
  for(j in 1:ncol(contam.all))
  {
    hh[i] <- 2.1*hh.yes
    hcw_conc[i,j] <- contam.all[i,j]*10^(-hh[i])
  }
}

hcw_conc <- as.data.frame(hcw_conc)

summary(hcw_conc)
summary(log10(hcw_conc))
colSums(hcw_conc==0)

rm(j,i)

#####
#Model HCW deposition onto fomite (cross-contamination)

for(i in 1:iters)
{
  for(j in 1:ncol(contam.all))
  {
    hcw_f_N[i] <- sample(c(1,2), 1, replace = TRUE) #restricted
    #hcw_f_N[i] <- rnbinom(1,1.83,0.167) #for risk characterization (daily risk)
    hcw_f_effc[i] <- rtriang(1,min=0.0067,max=0.1670,mode=0.0559)
    f_conc_cc[i,j] <- hcw_conc[i,j]*hcw_f_N[i]*hcw_f_effc[i]
  }
}

f_conc_cc <- as.data.frame(f_conc_cc)  

summary(f_conc_cc)
summary(log10(f_conc_cc))
#write.csv(f_conc_cc,file = "Fomite_CC_Conc.csv")

rm(j,i)

#####
#Model patient pick-up from fomite
f_p_N <-  rtriang(iters,min=5,max=31,mean=14.2)

is.dirt.cc = 'yes' #yes, assume single contact; no, probability contacted fomite is dirty (for daily risk)
if(is.dirt.cc == 'yes'){dirt.yes.cc = runif(iters,1,1)}else{dirt.yes.cc = rbinom(iters,floor(f_p_N),0.08)} 

for(i in 1:iters)
{
  for(j in 1:ncol(contam.all))
  {
    f_h_effc[i] <- rtriang(1,min=0.475,max=0.714,mean=0.57)
    f_t[i] <- runif(1,0,1)
    f_decay[i] <- decay_epd(k1,k2,f_t[i])
    h_conc[i,j] <- dirt.yes.cc[i]*f_conc_cc[i,j]*f_h_effc[i]*f_decay[i]
  }
}

h_conc <- as.data.frame(h_conc)

summary(h_conc)
summary(log10(h_conc))
colSums(h_conc==0)
#write.csv(h_conc,file = "Hand_Conc.csv")

rm(j,i)

#####
#Model patient pick-up and transfer to mouth

for(i in 1:iters)
{
  for(j in 1:ncol(contam.all))
  {
    h_area[i] <- runif(1,2,6)
    h_m_N[i] <- 1 #restricted
    #h_m_N[i] <- 0.59*12*rnormTrunc(1,2.9,2.5,min=0) #risk characterization (daily risk)
    h_m_effc[i] <- 0.4099
    dose[i,j] <- h_conc[i,j]*h_area[i]*h_m_N[i]*h_m_effc[i]
  }
}

dose <- as.data.frame(dose)

summary(dose)
summary(log10(dose))
#write.csv(dose,file = "Dose.csv")

rm(j,i)

#####
#Estimate risk of infection or illness (CDI)

for(i in 1:iters)
{
  for(j in 1:ncol(contam.all))
  {
    Risk.inf[i,j] <- p.inf.exact(dose[i,j],alpha.inf,beta.inf)
    Risk.ill[i,j] <- p.inf.exact(dose[i,j],alpha.ill,beta.ill)
  }
}

Risk.inf <- data.frame(Risk.inf)
Risk.ill <- data.frame(Risk.ill)

summary(Risk.inf)
summary(log10(Risk.inf))
summary(Risk.ill)
summary(log10(Risk.ill))

write.csv(Risk.inf,file = "Risk.inf.csv")
write.csv(Risk.ill,file = "Risk.ill.csv")
colSums(Risk.inf==0) #frequency of "zero" risk
colSums(Risk.ill==0)

rm(j,i)

#### Summary statistics and plotting ####

#Summary statistics 

CDI.descrips <- matrix(0,7,4)
CDI.ill.descrips <- matrix(0,7,4)

for(i in 1:ncol(Risk.inf))
{
  CDI.descrips[1,i] <- min(Risk.inf[,i])
  CDI.descrips[2,i] <- quantile(Risk.inf[,i], probs=0.05)
  CDI.descrips[3,i] <- median(Risk.inf[,i])
  CDI.descrips[4,i] <- quantile(Risk.inf[,i], probs=0.95)
  CDI.descrips[5,i] <- mean(Risk.inf[,i])
  CDI.descrips[6,i] <- sd(Risk.inf[,i])
  CDI.descrips[7,i] <- max(Risk.inf[,i])
}
colnames(CDI.descrips) <- colnames(Risk.inf)
rownames(CDI.descrips) <- c("Minimum","Lower 95th","Median","Upper 95th","Mean","Standard Deviation","Maximum")
write.csv(CDI.descrips,file = "CDI_Risk_Inf_Descriptive_Statistics.csv")

for(i in 1:ncol(Risk.ill))
{
  CDI.ill.descrips[1,i] <- min(Risk.ill[,i])
  CDI.ill.descrips[2,i] <- quantile(Risk.ill[,i], probs=0.05)
  CDI.ill.descrips[3,i] <- median(Risk.ill[,i])
  CDI.ill.descrips[4,i] <- quantile(Risk.ill[,i], probs=0.95)
  CDI.ill.descrips[5,i] <- mean(Risk.ill[,i])
  CDI.ill.descrips[6,i] <- sd(Risk.ill[,i])
  CDI.ill.descrips[7,i] <- max(Risk.ill[,i])
}
colnames(CDI.ill.descrips) <- colnames(Risk.ill)
rownames(CDI.ill.descrips) <- c("Minimum","Lower 95th","Median","Upper 95th","Mean","Standard Deviation","Maximum")
write.csv(CDI.ill.descrips,file = "CDI_Risk_Ill_Descriptive_Statistics.csv")

#Distribution of risk estimates 

Risk.inf.plot <- as.matrix(Risk.inf); 
Risk.inf.ggplot <- melt(Risk.inf.plot)
colnames(Risk.inf.ggplot) <- c("dummy","Source","Risk")

Risk.ill.plot <- as.matrix(Risk.ill); 
Risk.ill.ggplot <- melt(Risk.ill.plot)
colnames(Risk.ill.ggplot) <- c("dummy","Source","Risk")

#Histogram
CDI.hist <- ggplot(Risk.inf.ggplot, aes(x=Risk, fill=Source, color=Source)) + 
  geom_histogram(position = "identity", alpha=0.5, bins=40) +
  scale_x_log10() +
  labs(x = c(expression(~log[10]~Probability~of~Infection)),
       y = "Frequency", fill="Contamination Source", color="Contamination Source") +  
  theme_bw() +
  theme(text=element_text(family="serif",size=13)) +
  scale_color_manual(values=c("#7F3C8DFF","#11A579FF","#3969ACFF","#E73F74FF"), 
                     labels=c("Colonized Shed", "Infected Shed", "Pre-seeded", "HCW Pick up")) +
  scale_fill_manual(values=c("#7F3C8DFF","#11A579FF","#3969ACFF","#E73F74FF"), 
                    labels=c("Colonized Shed", "Infected Shed", "Pre-seeded", "HCW Pick up"))
ggsave("CDI_Risk_Inf_hist.png",CDI.hist,dpi=600)

CDI.ill.hist <- ggplot(Risk.ill.ggplot, aes(x=Risk, fill=Source, color=Source)) + 
  geom_histogram(position = "identity", alpha=0.5, bins=40) +
  scale_x_log10() +
  labs(x = c(expression(~log[10]~Probability~of~Illness)),
       y = "Frequency", fill="Contamination Source", color="Contamination Source") +  
  theme_bw() +
  theme(text=element_text(family="serif",size=13)) +
  scale_color_manual(values=c("#7F3C8DFF","#11A579FF","#3969ACFF","#E73F74FF"), 
                     labels=c("Colonized Shed", "Infected Shed", "Pre-seeded", "HCW Pick up")) +
  scale_fill_manual(values=c("#7F3C8DFF","#11A579FF","#3969ACFF","#E73F74FF"), 
                    labels=c("Colonized Shed", "Infected Shed", "Pre-seeded", "HCW Pick up"))
ggsave("CDI_Risk_Ill_hist.png",CDI.ill.hist,dpi=600)

#Box/Violin plot
CDI.box <- ggplot(Risk.inf.ggplot, aes(Source, Risk, fill=Source)) + 
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
ggsave("CDI_Risk_Inf_box.png",CDI.box)

CDI.ill.box <- ggplot(Risk.ill.ggplot, aes(Source, Risk, fill=Source)) + 
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
ggsave("CDI_Risk_Ill_box.png",CDI.ill.box)

#Jitterplot
CDI.jitter <- ggplot(Risk.inf.ggplot, aes(Source,Risk)) +
  geom_jitter(aes(colour=Source), alpha=0.75) +
  scale_y_log10() +
  labs(x = "Environmental Contamination Source",
       y = c(expression(~log[10]~Probability~of~Infection))) +
  theme_bw() +
  theme(text=element_text(family="serif",size=13)) +
  scale_color_manual(values=c("#7F3C8DFF","#11A579FF","#3969ACFF","#E73F74FF"), 
                     labels=c("Colonized Shed", "Infected Shed", "Pre-seeded", "HCW Pick up"))
ggsave("CDI_Risk_Inf_jitter-HH.png",CDI.jitter,dpi=600)

CDI.ill.jitter <- ggplot(Risk.ill.ggplot, aes(Source,Risk)) +
  geom_jitter(aes(colour=Source), alpha=0.75) +
  scale_y_log10() +
  labs(x = "Environmental Contamination Source",
       y = c(expression(~log[10]~Probability~of~Infection))) +
  theme_bw() +
  theme(text=element_text(family="serif",size=13)) +
  scale_color_manual(values=c("#7F3C8DFF","#11A579FF","#3969ACFF","#E73F74FF"), 
                     labels=c("Colonized Shed", "Infected Shed", "Pre-seeded", "HCW Pick up"))
ggsave("CDI_Risk_Ill_jitter-HH.png",CDI.ill.jitter)


#END

