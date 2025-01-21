
#QMRA model for Candida auris catheter-associated urinary tract infections - Sensitivity analysis

#Created: 1 Aug 2024
#Last modified: 21 Jan 2025

#####
#Plot parameter values against risk estimates (by source of contamination)

#Shedding Low
par(mfrow = c(1,4), cex=1, oma=c(0,0,0,0))
plot(p_col_conc1,Risk.contam$`Shed Low`, log="y", xlab="Skin Colonization Concentration", ylab="CAUTI Risk")
lines(lowess(p_col_conc1,Risk.contam$`Shed Low`), lwd = 4, col = "#de712b")
plot(t_shed,Risk.contam$`Shed Low`, log="y", xlab="Duration of shedding", ylab="CAUTI Risk")
lines(lowess(t_shed,Risk.contam$`Shed Low`), lwd = 4, col = "#de712b")
plot(f_decay,Risk.contam$`Shed Low`, log="y", xlab="Fomite decay", ylab="CAUTI Risk")
lines(lowess(f_decay,Risk.contam$`Shed Low`), lwd = 4, col = "#de712b")
plot(f_h_effc,Risk.contam$`Shed Low`, log="y", xlab="Fomite-to-hand transfer efficiency", ylab="CAUTI Risk")
lines(lowess(f_h_effc,Risk.contam$`Shed Low`), lwd = 4, col = "#de712b")

par(mfrow = c(1,4), cex=1, oma=c(0,0,0,0))
plot(f_h_N,Risk.contam$`Shed Low`, log="y", xlab="No. of fomite contacts", ylab="CAUTI Risk")
lines(lowess(f_h_N,Risk.contam$`Shed Low`), lwd = 4, col = "#de712b")
plot(f_t,Risk.contam$`Shed Low`, log="y", xlab="Time to fomite contact", ylab="CAUTI Risk")
lines(lowess(f_t,Risk.contam$`Shed Low`), lwd = 4, col = "#de712b")
plot(h_c_effc,Risk.contam$`Shed Low`, log="y", xlab="Hand-to-catheter transfer efficiency", ylab="CAUTI Risk")
lines(lowess(h_c_effc,Risk.contam$`Shed Low`), lwd = 4, col = "#de712b")
plot(h_decay,Risk.contam$`Shed Low`, log="y", xlab="Hand decay", ylab="CAUTI Risk")
lines(lowess(h_decay,Risk.contam$`Shed Low`), lwd = 4, col = "#de712b")

par(mfrow = c(1,2), cex=1, oma=c(0,0,0,0))
plot(h_t,Risk.contam$`Shed Low`, log="y", xlab="Time to catheter contact", ylab="CAUTI Risk")
lines(lowess(h_t,Risk.contam$`Shed Low`), lwd = 4, col = "#de712b")
plot(k,Risk.contam$`Shed Low`, log="y", xlab="Dose-response parameter", ylab="CAUTI Risk")
lines(lowess(k,Risk.contam$`Shed Low`), lwd = 4, col = "#de712b")

#Shedding High
par(mfrow = c(1,4), cex=1, oma=c(0,0,0,0))
plot((10^p_col_conc2),Risk.contam$`Shed High`, log="y", xlab="Skin Colonization Concentration", ylab="CAUTI Risk")
lines(lowess((10^p_col_conc2),Risk.contam$`Shed High`), lwd = 4, col = "#de712b")
plot(t_shed,Risk.contam$`Shed High`, log="y", xlab="Duration of shedding", ylab="CAUTI Risk")
lines(lowess(t_shed,Risk.contam$`Shed High`), lwd = 4, col = "#de712b")
plot(f_decay,Risk.contam$`Shed High`, log="y", xlab="Fomite decay", ylab="CAUTI Risk")
lines(lowess(f_decay,Risk.contam$`Shed High`), lwd = 4, col = "#de712b")
plot(f_h_effc,Risk.contam$`Shed High`, log="y", xlab="Fomite-to-hand transfer efficiency", ylab="CAUTI Risk")
lines(lowess(f_h_effc,Risk.contam$`Shed High`), lwd = 4, col = "#de712b")

par(mfrow = c(1,4), cex=1, oma=c(0,0,0,0))
plot(f_h_N,Risk.contam$`Shed High`, log="y", xlab="No. of fomite contacts", ylab="CAUTI Risk")
lines(lowess(f_h_N,Risk.contam$`Shed High`), lwd = 4, col = "#de712b")
plot(f_t,Risk.contam$`Shed High`, log="y", xlab="Time to fomite contact", ylab="CAUTI Risk")
lines(lowess(f_t,Risk.contam$`Shed High`), lwd = 4, col = "#de712b")
plot(h_c_effc,Risk.contam$`Shed High`, log="y", xlab="Hand-to-catheter transfer efficiency", ylab="CAUTI Risk")
lines(lowess(h_c_effc,Risk.contam$`Shed High`), lwd = 4, col = "#de712b")
plot(h_decay,Risk.contam$`Shed High`, log="y", xlab="Hand decay", ylab="CAUTI Risk")
lines(lowess(h_decay,Risk.contam$`Shed High`), lwd = 4, col = "#de712b")

par(mfrow = c(1,2), cex=1, oma=c(0,0,0,0))
plot(h_t,Risk.contam$`Shed High`, log="y", xlab="Time to catheter contact", ylab="CAUTI Risk")
lines(lowess(h_t,Risk.contam$`Shed High`), lwd = 4, col = "#de712b")
plot(k,Risk.contam$`Shed High`, log="y", xlab="Dose-response parameter", ylab="CAUTI Risk")
lines(lowess(k,Risk.contam$`Shed High`), lwd = 4, col = "#de712b")

#Pre-seeded fomite
par(mfrow = c(1,3), cex=1, oma=c(0,0,0,0))
plot(f_dirty_conc_all,Risk.contam$All, log="y", xlab="Environment Concentration", ylab="CAUTI Risk")
lines(lowess(f_dirty_conc_all,Risk.contam$All), lwd = 4, col = "#de712b")
plot(f_decay,Risk.contam$All, log="y", xlab="Fomite decay", ylab="CAUTI Risk")
lines(lowess(f_decay,Risk.contam$All), lwd = 4, col = "#de712b")
plot(f_h_effc,Risk.contam$All, log="y", xlab="Fomite-to-hand transfer efficiency", ylab="CAUTI Risk")
lines(lowess(f_h_effc,Risk.contam$All), lwd = 4, col = "#de712b")

par(mfrow = c(1,4), cex=1, oma=c(0,0,0,0))
plot(f_h_N,Risk.contam$All, log="y", xlab="No. of fomite contacts", ylab="CAUTI Risk")
lines(lowess(f_h_N,Risk.contam$All), lwd = 4, col = "#de712b")
plot(f_t,Risk.contam$All, log="y", xlab="Time to fomite contact", ylab="CAUTI Risk")
lines(lowess(f_t,Risk.contam$All), lwd = 4, col = "#de712b")
plot(h_c_effc,Risk.contam$All, log="y", xlab="Hand-to-catheter transfer efficiency", ylab="CAUTI Risk")
lines(lowess(h_c_effc,Risk.contam$All), lwd = 4, col = "#de712b")
plot(h_decay,Risk.contam$All, log="y", xlab="Hand decay", ylab="CAUTI Risk")
lines(lowess(h_decay,Risk.contam$All), lwd = 4, col = "#de712b")

par(mfrow = c(1,2), cex=1, oma=c(0,0,0,0))
plot(h_t,Risk.contam$All, log="y", xlab="Time to catheter contact", ylab="CAUTI Risk")
lines(lowess(h_t,Risk.contam$All), lwd = 4, col = "#de712b")
plot(k,Risk.contam$All, log="y", xlab="Dose-response parameter", ylab="CAUTI Risk")
lines(lowess(k,Risk.contam$All), lwd = 4, col = "#de712b")

#####
#Tornado plots by source of contamination

#Shedding Low
Risk.Vars <- data.frame(p_col_conc1, t_shed, f_decay, f_h_effc, f_h_N, f_t, h_c_effc, h_decay, h_t, k) 
colnames(Risk.Vars) <- c("c_col", "t_shed", "decay_f","eff_f_h", "N_f_h","t_f", "eff_h_c", "decay_h", "t_h", "k")
Sens.Risk <- matrix(NA, 1, ncol(Risk.Vars)); colnames(Sens.Risk) <- colnames(Risk.Vars)
for(j in 1:ncol(Sens.Risk))
{
  Sens.Risk[,j] <- cor(log(Risk.Vars[,j],10),log(Risk.contam$`Shed Low`,10), method = "spearman", use="complete.obs")
}
Risk.Sens <- melt(Sens.Risk); colnames(Risk.Sens) <- c("dummy","Variable","Rho")
Risk.Sens.PLOT <- ggplot(Risk.Sens,aes(x=Variable,y=Rho, fill=Variable)) + 
  geom_bar(stat = "identity", color = "black", fill = "grey") + 
  coord_flip() + 
  theme_bw() +
  theme(text=element_text(family="serif",size=13), legend.position = "none") + 
  labs(y = expression(~Spearman~rho), x = "Model Variable", title = " ")
ggsave("Risk.Sens.PLOT_Low.png",Risk.Sens.PLOT,dpi=600)

#Shedding High
Risk.Vars <- data.frame(p_col_conc2, t_shed, f_decay, f_h_effc, f_h_N, f_t, h_c_effc, h_decay, h_t, k) 
colnames(Risk.Vars) <- c("c_col", "t_shed", "decay_f","eff_f_h", "N_f_h","t_f", "eff_h_c", "decay_h", "t_h", "k")
Sens.Risk <- matrix(NA, 1, ncol(Risk.Vars)); colnames(Sens.Risk) <- colnames(Risk.Vars)
for(j in 1:ncol(Sens.Risk))
{
  Sens.Risk[,j] <- cor(log(Risk.Vars[,j],10),log(Risk.contam$`Shed High`,10), method = "spearman", use="complete.obs")
}
Risk.Sens <- melt(Sens.Risk); colnames(Risk.Sens) <- c("dummy","Variable","Rho")
Risk.Sens.PLOT <- ggplot(Risk.Sens,aes(x=Variable,y=Rho, fill=Variable)) + 
  geom_bar(stat = "identity", color = "black", fill = "grey") + 
  coord_flip() + 
  theme_bw() +
  theme(text=element_text(family="serif",size=13), legend.position = "none") + 
  labs(y = expression(~Spearman~rho), x = "Model Variable", title = " ")
ggsave("Risk.Sens.PLOT_High.png",Risk.Sens.PLOT,dpi=600)

#Pre-seeded fomite
Risk.Vars <- data.frame(f_dirty_conc_all, f_h_effc, f_h_N, h_c_effc, h_t, k) 
colnames(Risk.Vars) <- c("Fomite concentration","Fomite-hand transfer effc", 
                         "HCW-fomite contact events", "Hand-catheter transfer effc", "Time to catheter contact", "k")
Sens.Risk <- matrix(NA, 1, ncol(Risk.Vars)); colnames(Sens.Risk) <- colnames(Risk.Vars)
for(j in 1:ncol(Sens.Risk))
{
  Sens.Risk[,j] <- cor(log(Risk.Vars[,j],10),log(Risk.contam$All,10), method = "spearman", use="complete.obs")
}
Risk.Sens <- melt(Sens.Risk); colnames(Risk.Sens) <- c("dummy","Variable","Rho")
Risk.Sens.PLOT <- ggplot(Risk.Sens,aes(x=Variable,y=Rho, fill=Variable)) + 
  geom_bar(stat = "identity", color = "black", fill = "grey", width = 0.5) + 
  coord_flip() + 
  theme_bw() +
  theme(text=element_text(family="serif",size=13), legend.position = "none", aspect.ratio = 3/2, 
        axis.title.y = element_text(margin = margin(r = 20))) + 
  labs(y = expression(~Spearman~rho), x = "Model Variable", title = " ") +
  geom_hline(yintercept = 0)
ggsave("Risk.Sens.PLOT_All.png",Risk.Sens.PLOT,dpi=600)


#####
#Alternative parameterization

#Contamination source
shed_rate <- 0.3 #MRSA shedding rate
t_shed <-runif(iters,0,24)
f_dirty_conc_shed <- vector()

for(i in 1:iters)
{
f_dirty_conc_shed[i] <- shed_effc*t_shed[i]*shed_rate
}

f_dirty_conc_all2 <- rtriangle(iters, 0.024,9.70,0.92) #alt fomite concentration

contam2 <- c("Shed MRSA", "Fomite")
N.contam2 <- length(contam2)
contam.all <- data.frame(f_dirty_conc_shed, f_dirty_conc_all2)
colnames(contam.all) <- contam2

#Transfer efficiency
f_h_effc[i] <- runif(1, 0.034,0.203) #alt TE under low RH
f_h_effc[i] <- runif(1, 0.396,0.619) #alt TE under high RH

h_c_effc[i] <- rtriangle(1,0,0.3447,0.1026) #alt TE (hands)

#Persistence

#Dampened exponential
decay_epd = function(K1,K2,t){
  out = exp(-K1*t*exp(-K2*t))
}

#Gompertz
decay_gz = function(K1,K2,t){	
  out = exp(-k1/k2*(exp(k2*t)-1))
}

#Plastic decay
c <- 1.2533
d <- 0.0407
f_decay[i] <- decay_epd(k1,k2,f_t[i])

#Steel growth
a <- -0.4002 
b <- 0.1102
f_growth[i] <- decay_epd(k1,k2,f_t[i])

#####
#Tornado plots by alt source of contamination

#MRSA shedding parameter
Risk.Vars <- data.frame(f_dirty_conc_shed, t_shed, f_decay, f_h_effc, f_h_N, f_t, h_c_effc, h_decay, h_t, k) 
colnames(Risk.Vars) <- c("f_conc_shed", "t_shed", "decay_f","eff_f_h", "N_f_h","t_f", "eff_h_c", "decay_h", "t_h", "k")
Sens.Risk <- matrix(NA, 1, ncol(Risk.Vars)); colnames(Sens.Risk) <- colnames(Risk.Vars)
for(j in 1:ncol(Sens.Risk))
{
  Sens.Risk[,j] <- cor(log(Risk.Vars[,j],10),log(Risk.contam$`Shed MRSA`,10), method = "spearman", use="complete.obs")
}
Risk.Sens <- melt(Sens.Risk); colnames(Risk.Sens) <- c("dummy","Variable","Rho")
Risk.Sens.PLOT <- ggplot(Risk.Sens,aes(x=Variable,y=Rho, fill=Variable)) + 
  geom_bar(stat = "identity", color = "black", fill = "grey") + 
  coord_flip() + 
  theme_bw() +
  theme(text=element_text(family="serif",size=13), legend.position = "none") + 
  labs(y = expression(~Spearman~rho), x = "Model Variable", title = " ")
Risk.Sens.PLOT

#Pre-seeded fomite from environmental sampling data
Risk.Vars <- data.frame(f_dirty_conc_all2,f_decay, f_h_effc, f_h_N, f_t, h_c_effc, h_decay, h_t, k) 
colnames(Risk.Vars) <- c("f_conc", "decay_f","eff_f_h", "N_f_h","t_f", "eff_h_c", "decay_h", "t_h", "k")
Sens.Risk <- matrix(NA, 1, ncol(Risk.Vars)); colnames(Sens.Risk) <- colnames(Risk.Vars)
for(j in 1:ncol(Sens.Risk))
{
  Sens.Risk[,j] <- cor(log(Risk.Vars[,j],10),log(Risk.contam$`Fomite`,10), method = "spearman", use="complete.obs")
}
Risk.Sens <- melt(Sens.Risk); colnames(Risk.Sens) <- c("dummy","Variable","Rho")
Risk.Sens.PLOT <- ggplot(Risk.Sens,aes(x=Variable,y=Rho, fill=Variable)) + 
  geom_bar(stat = "identity", color = "black", fill = "grey") + 
  coord_flip() + 
  theme_bw() +
  theme(text=element_text(family="serif",size=13), legend.position = "none") + 
  labs(y = expression(~Spearman~rho), x = "Model Variable", title = " ")
Risk.Sens.PLOT

#END
