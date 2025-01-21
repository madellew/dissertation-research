
#QMRA model for Candida auris central line-associated bloodstream infections

#Created: 20 May 2024
#Last modified: 21 Jan 2025

#Follows cauti_qmra.R up to dose estimation (line 123)

#Read CSV for bootstrapped dose-response parameters
k_bp_clabsi <- read.csv("exp_bp_b311.csv")

#####
#Estimate concentration transferred to central venous catheter (dose)

dose.clabsi <- matrix(0, iters, ncol = N.contam)
colnames(dose.clabsi) <- contam

#Can skip to dose-response if assume break in sterility always occurs, if so...
dose.clabsi <- dose

h_area <- vector() #fingertip surface area
h_c_N <- vector() #hand-catheter transfer frequency
h_c_effc <- vector() #hand-catheter transfer efficiency
h_decay <- vector() #decay rate - hands
h_t <- vector() #time between fomite contacted and catheter contact

is.sterile = 'yes' #yes, probability break in sterile barriers; no, assume all contacts
if(is.sterile == 'yes'){sterile.yes = rbinom(iters,1,0.62)}else{sterile.yes = runif(iters,1,1)} 

for(i in 1:iters)
{
  for(j in 1:ncol(contam.all))
  {
    h_area[i] <- 10
    h_c_N[i] <- sample(c(1,2), 1, replace = TRUE)
    h_c_effc[i] <- rtriangle(1,0.00007143,0.06429,0.008131)
    h_t[i] <- rtriangle(1,3.33,120,13.33)
    h_decay[i] <- exp(-0.07296*h_t[i])
    dose.clabsi[i,j] <- sterile.yes[i]*h_conc[i,j]*h_area[i]*h_c_N[i]*h_c_effc[i]*h_decay[i]
  }
}

dose.clabsi <- as.data.frame(dose.clabsi)

summary(dose.clabsi)
summary(log10(dose.clabsi))
write.csv(dose,file = "Dose_CLABSI.csv")

rm(j,i)

#####
#Estimate risk of infection (CLABSI) using exponential dose-response model

Risk.contam.clabsi <- matrix(0, iters, ncol = N.contam)  #risk of infection
colnames(Risk.contam.clabsi) <- contam

for(i in 1:iters)
{
  for(j in 1:ncol(contam.all))
  {
    k_bsi <- k_bp_clabsi$k
    Risk.contam.clabsi[i,j] <- 1-exp(-k_bsi[i]*dose.clabsi[i,j]) #per HCW visit to patient
  }
}

Risk.contam.clabsi <- as.data.frame(Risk.contam.clabsi)

rm(j,i)

summary(Risk.contam.clabsi)
summary(log10(Risk.contam.clabsi))
colSums(Risk.contam.clabsi==0) 
write.csv(Risk.contam.clabsi,file = "Risk_CLASBI.csv")


#### Summary statistics and plotting ####

#Summary statistics 

CLABSI.descrips <- matrix(0,7,3)

for(i in 1:ncol(Risk.contam.clabsi))
{
  CLABSI.descrips[1,i] <- min(Risk.contam.clabsi[,i])
  CLABSI.descrips[2,i] <- quantile(Risk.contam.clabsi[,i], probs=0.05)
  CLABSI.descrips[3,i] <- median(Risk.contam.clabsi[,i])
  CLABSI.descrips[4,i] <- quantile(Risk.contam.clabsi[,i], probs=0.95)
  CLABSI.descrips[5,i] <- mean(Risk.contam.clabsi[,i])
  CLABSI.descrips[6,i] <- sd(Risk.contam.clabsi[,i])
  CLABSI.descrips[7,i] <- max(Risk.contam.clabsi[,i])
}
colnames(CLABSI.descrips) <- colnames(Risk.contam.clabsi)
rownames(CLABSI.descrips) <- c("Minimum","Lower 95th","Median","Upper 95th","Mean","Standard Deviation","Maximum")
write.csv(CLABSI.descrips,file = "CLABSI_Risk_Descriptive_Statistics.csv")

#Distribution of risk estimates 

Risk.contam.clabsi.plot <- as.matrix(Risk.contam.clabsi); 
Risk.contam.clabsi.ggplot <- melt(Risk.contam.clabsi.plot)
colnames(Risk.contam.clabsi.ggplot) <- c("dummy","Source","Risk")

#Histogram
CLABSI.hist <- ggplot(Risk.contam.clabsi.ggplot, aes(x=Risk, fill=Source, color=Source)) + 
  geom_histogram(position = "identity", alpha=0.5, bins=40) +
  scale_x_log10() +
  labs(x = c(expression(~log[10]~Probability~of~Infection)),
       y = "Frequency", fill="Contamination Source", color="Contamination Source") +  
  theme_bw() +
  theme(text=element_text(family="serif",size=13)) +
  scale_color_paletteer_d("rcartocolor::Bold", labels=c("Low Shed", "High Shed", "Pre-seeded")) +
  scale_fill_paletteer_d("rcartocolor::Bold", labels=c("Low Shed", "High Shed", "Pre-seeded"))
ggsave("CLABSI_Risk_hist.png",CLABSI.hist,dpi=600)

#Boxplot
CLABSI.box <- ggplot(Risk.contam.clabsi.ggplot, aes(Source, Risk, fill=Source)) + 
  geom_boxplot() + 
  scale_y_log10() +
  theme(axis.text.x = element_text(face="bold", color = "black",size = 10, angle=45)) + 
  labs(x = c(expression("Environmental Contamination Source")), 
       y = expression(~log[10]~Probability~of~Infection),fill="Source") + 
  stat_summary(fun=median, geom="point", size=2.5, color="white") + 
  stat_summary(fun=mean, geom="point", size=2, color="black",shape=15) + 
  guides(fill = guide_legend(override.aes = list(shape = NA)))
ggsave("CLABSI_Risk_box.png",CLABSI.box)

#Jitter plot
CLABSI.jitter <- ggplot(Risk.contam.clabsi.ggplot, aes(Source,Risk)) +
  geom_jitter(aes(colour=Source)) +
  scale_y_log10() +
  labs(x = c(expression("Environmental Contamination Source")),
       y = c(expression(~log[10]~Probability~of~Infection))) 
ggsave("CLABSI_Risk_jitter.png",CLABSI.jitter)


#End