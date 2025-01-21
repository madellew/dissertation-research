
#QMRA model for Candida auris device-related infections - Microbial growth

#Created: 15 Aug 2024
#Last modified: 21 Jan 2025

#Following the initial exposure event and subsequent dose, the dose resulting through the duration of catheterization
#is a function of the initial exposure dose, microbial growth rate, and the length of catheterization.
#Risk following microbial growth must be estimated for each source of contamination individually.

#####
#Specify growth model (function and parameters)

#Dampened exponential
decay_epd = function(K1,K2,t){
  out = exp(-K1*t*exp(-K2*t))
}

k_1 <- -0.2125
k_2 <- 0.0497

#####
#Set up matrices to estimate dose and risk everyday for 40 days

t <- vector()
t <- c(0:40)
k <- k_bp_cauti$k

dose_grow <- matrix(0,iters,ncol=41)
colnames(dose_grow) <- paste0(0:40)
Risk.grow <- matrix(0,iters,ncol=41)
colnames(Risk.grow) <- paste0(0:40)

#Estimate daily risk

dose_vector <- vector()

for(i in 1:iters)
{
  dose_vector = dose$All[i]*decay_epd(k_1,k_2,t)
  #dose_vector = dose$All[i]*decay_gz(k1,k2,t) #alt growth model
  for(j in 1:ncol(Risk.grow))
  {
    dose_grow[i,j] <- dose_vector[j]
    Risk.grow[i,j] <- 1-exp(-k[i]*dose_grow[i,j]) #per HCW visit to patient
  }
}

Risk.grow <- as.data.frame(Risk.grow)
rm(j,i)

#####
#Calculate daily median and mean risk and define for each contamination source

colMed_Shed1 <- apply(Risk.grow,2,median)
colMed_Shed2 <- apply(Risk.grow,2,median)
colMed_All <- apply(Risk.grow,2,median)

colMeans_Shed1 <- colMeans(Risk.grow)
colMeans_Shed2 <- colMeans(Risk.grow)
colMeans_All <- colMeans(Risk.grow)

grow.mean <- data.frame(colMeans_Shed1,colMeans_Shed2,colMeans_All)
colnames(grow.mean) <- contam
write.csv(grow.mean,file = "CAUTI.grow.mean.csv")

grow.med <- data.frame(colMed_Shed1,colMed_Shed2,colMed_All)
colnames(grow.med) <- contam
write.csv(grow.med,file = "CAUTI.grow.med.csv")

#####
#Plot daily mean/median risk of infection by source of contamination

#Mean
grow.plot <- as.matrix(grow.mean) 
grow.ggplot <- melt(grow.plot)
t_vec <-rep(t, times=3)
grow.ggplot$t <- t_vec
colnames(grow.ggplot) <- c("dummy","Source","Mean","t")

CAUTI.growth <- ggplot(grow.ggplot, aes(x=t, y=Mean, group=Source, color=Source)) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  labs(x = "Environmental Contamination Source",
       y = c(expression(~log[10]~Mean~Probability~of~Infection)), color="Contamination Source") +
  theme_bw() +
  theme(text=element_text(family="serif",size=13)) +
  scale_color_paletteer_d("rcartocolor::Bold", labels=c("Low Shed", "High Shed", "Pre-seeded"))
ggsave("CAUTI_Growth_logMean.png",CAUTI.growth,dpi=600)

#Median
grow.plot.med <- as.matrix(grow.med) 
grow.ggplot.med <- melt(grow.plot.med)
grow.ggplot.med$t <- t_vec
colnames(grow.ggplot.med) <- c("dummy","Source","Median","t")

CAUTI.growth.Med <- ggplot(grow.ggplot.med, aes(x=t, y=Median, group=Source, color=Source)) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  labs(x = "Environmental Contamination Source",
       y = c(expression(~log[10]~Median~Probability~of~Infection)), color="Contamination Source") +
  theme_bw() +
  theme(text=element_text(family="serif",size=13)) +
  scale_color_paletteer_d("rcartocolor::Bold", labels=c("Low Shed", "High Shed", "Pre-seeded"))
ggsave("CAUTI_Growth_logMed.png",CAUTI.growth.Med,dpi=600)


#####
#Growth by exposure duration group

#Define exposure groups
exp_dur <- c("Initial", "Short","Moderate", "Extended", "Prolonged")
N.exp_dur <- length(exp_dur)

exp_dur_1 <- 1
exp_dur_s <- rtriangle(iters, 1,7,3)
exp_dur_m <- rtriangle(iters, 8,14,11.5)
exp_dur_l <- rtriangle(iters, 15,28,25)
exp_dur_xl <- rtriangle(iters,30,43,36)
exp_dur.all <- data.frame(exp_dur_1, exp_dur_s, exp_dur_m, exp_dur_l, exp_dur_xl)
colnames(exp_dur.all) <- exp_dur

grow <- matrix(0, iters, ncol = N.exp_dur)
dose_grow_LOS <- matrix(0, iters, ncol = N.exp_dur)
Risk_grow_LOS <- matrix(0, iters, ncol = N.exp_dur) 

for(i in 1:iters)
{
  for(j in 1:ncol(exp_dur.all))
  {
    k <- k_bp_cauti$k
    grow[i,j] <- decay_epd(k_1,k_2,exp_dur.all[i,j])
    dose_grow_LOS[i,j] <- dose$All[i]*grow[i,j]
    Risk_grow_LOS[i,j] <- 1-exp(-k[i]*dose_grow_LOS[i,j]) 
  }
}

Risk_grow_LOS <- as.data.frame(Risk_grow_LOS)
colnames(Risk_grow_LOS) <- exp_dur

summary(Risk_grow_LOS)
summary(log10(Risk_grow_LOS))
colSums(Risk_grow_LOS==0)

rm(j,i)

Risk.grow.LOS.plot <- as.matrix(Risk_grow_LOS); 
Risk.grow.LOS.ggplot <- melt(Risk.grow.LOS.plot)
colnames(Risk.grow.LOS.ggplot) <- c("dummy","Duration","Risk")

Risk.LOS.box.out <- ggplot(Risk.grow.LOS.ggplot, aes(Duration, Risk, fill=Duration)) + 
  #geom_boxplot(alpha=0.75) + 
  geom_violin(alpha=0.75, trim=FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(text=element_text(family="serif",size=13), axis.text.x=element_text(size = 12)) +
  labs(x = "Length of catheterization", 
       y = expression(~log[10]~Probability~of~Infection), fill="Duration") +
  #stat_summary(fun=median, geom="point", size=2.5, shape=5) + 
  #stat_summary(fun=mean, geom="point", size=2, color="black",shape=16) + 
  stat_summary(fun.data="mean_sdl", fun.args = list(mult = 1), geom="pointrange") +
  scale_fill_paletteer_d("rcartocolor::Bold")
ggsave("CAUTI_Risk_LOS_All.png",Risk.LOS.box.out,dpi=600)


#END
