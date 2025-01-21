
#QMRA model for Candida auris device-related infections - Cumulative Risk

#Created: 1 Aug 2024
#Last modified: 21 Jan 2025

#Multiple exposure days are modeled to estimate cumulative risk, the probability that infection will occur over n days.
#Exposure days are group based on duration of catheterization.
#Cumulative risk is estimated for a single source of contamination at a time.

#####
#Define CAUTI Exposure Groups

exp_dur <- c("Initial", "Short","Moderate", "Extended", "Prolonged")
N.exp_dur <- length(exp_dur)

exp_dur_1 <- 1 # this should replicate initial risk
exp_dur_s <- rtriangle(iters, 1,7,3) #min, max, mean
exp_dur_m <- rtriangle(iters, 8,14,11.5)
exp_dur_l <- rtriangle(iters, 15,28,25)
exp_dur_xl <- rtriangle(iters,30,43,36)
exp_dur.all <- data.frame(exp_dur_1, exp_dur_s, exp_dur_m, exp_dur_l, exp_dur_xl)
colnames(exp_dur.all) <- exp_dur

#####
#Estimate cumulative risk for specified contamination source

Risk.cum <- matrix(0, iters, ncol = N.exp_dur)
colnames(Risk.cum) <- exp_dur

for(i in 1:iters)
{
  for(j in 1:ncol(exp_dur.all))
  {
    Risk.cum[i,j] <- 1-(1-Risk.contam$`All`[i])^exp_dur.all[i,j]
  }
}

Risk.cum <- as.data.frame(Risk.cum)

summary(Risk.cum)
summary(log10(Risk.cum))
colSums(Risk.contam==0)

#####
#Box/violin plot of cumulative risk over exposure duration

Risk.cum.plot <- as.matrix(Risk.cum) 
Risk.cum.ggplot <- melt(Risk.cum.plot)
colnames(Risk.cum.ggplot) <- c("dummy","Duration","Risk")

Risk.cum.box.out <- ggplot(Risk.cum.ggplot, aes(Duration, Risk, fill=Duration)) + 
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
ggsave("CAUTI_Cum_Risk_All.png",Risk.cum.box.out,dpi=600)


#####
#Define CLABSI exposure groups

exp_dur <- c("Initial", "Short","Moderate", "Extended")
N.exp_dur <- length(exp_dur)

exp_dur_1 <- 1 # this should replicate initial risk
exp_dur_s <- rtriangle(iters, 2,10,6.4) #min, max, mean
exp_dur_m <- rtriangle(iters, 11,20,12.7)
exp_dur_l <- rtriangle(iters, 20,39,21.7)
exp_dur.all <- data.frame(exp_dur_1, exp_dur_s, exp_dur_m, exp_dur_l)
colnames(exp_dur.all) <- exp_dur

#####
#Estimate cumulative risk for specified contamination source

Risk.cum <- matrix(0, iters, ncol = N.exp_dur)
colnames(Risk.cum) <- exp_dur

for(i in 1:iters)
{
  for(j in 1:ncol(exp_dur.all))
  {
    Risk.cum[i,j] <- 1-(1-Risk.contam.clabsi$All[i])^exp_dur.all[i,j]
  }
}

Risk.cum <- as.data.frame(Risk.cum)

summary(Risk.cum)
summary(log10(Risk.cum))
colSums(Risk.cum==0)

#####
#Box/violin plot of cumulative risk over exposure duration

Risk.cum.plot <- as.matrix(Risk.cum)
Risk.cum.ggplot <- melt(Risk.cum.plot)
colnames(Risk.cum.ggplot) <- c("dummy","Duration","Risk")

Risk.cum.box.out <- ggplot(Risk.cum.ggplot, aes(Duration, Risk, fill=Duration)) + 
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
ggsave("CLABSI_Cum_Risk_All.png",Risk.cum.box.out,dpi=600)


#END

