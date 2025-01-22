
#QMRA model for Clostridioides difficile infection (CDI) - Sensitivity analysis

#Created: 1 Aug 2024
#Last modified: 21 Jan 2025

#####
#Plot parameter values against risk estimates (by source of contamination)

#Shedding - colonized
par(mfrow = c(1,4), cex=1, oma=c(0,0,0,0))
plot(p_conc_col,Risk.inf$Shed.Col, log="y", xlab="Colonized Shedding Concentration", ylab="CDI Risk")
lines(lowess(p_conc_col,Risk.inf$Shed.Col), lwd = 4, col = "#de712b")
plot(f_decay,Risk.inf$Shed.Col, log="y", xlab="Fomite decay", ylab="CDI Risk")
lines(lowess(f_decay,Risk.inf$Shed.Col), lwd = 4, col = "#de712b")
plot(f_t,Risk.inf$Shed.Col, log="y", xlab="Time to fomite contact", ylab="CDI Risk")
lines(lowess(f_t,Risk.inf$Shed.Col), lwd = 4, col = "#de712b")
plot(f_h_effc,Risk.inf$Shed.Col, log="y", xlab="Fomite-to-hand transfer efficiency", ylab="CDI Risk")
lines(lowess(f_h_effc,Risk.inf$Shed.Col), lwd = 4, col = "#de712b")

par(mfrow = c(1,3), cex=1, oma=c(0,0,0,0))
plot(f_h_N,Risk.inf$Shed.Col, log="y", xlab="No. of HCW fomite contacts", ylab="CDI Risk")
lines(lowess(f_h_N,Risk.inf$Shed.Col), lwd = 4, col = "#de712b")
plot(hcw_f_effc,Risk.inf$Shed.Col, log="y", xlab="Hand-to-fomite transfer efficiency", ylab="CDI Risk")
lines(lowess(hcw_f_effc,Risk.inf$Shed.Col), lwd = 4, col = "#de712b")
plot(h_area,Risk.inf$Shed.Col, log="y", xlab="Fingertip surface area", ylab="CDI Risk")
lines(lowess(h_area,Risk.inf$Shed.Col), lwd = 4, col = "#de712b")

par(mfrow = c(1,2), cex=1, oma=c(0,0,0,0))
plot(f_p_N,Risk.inf$Shed.Col, log="y", xlab="No. of patient fomite contacts", ylab="CDI Risk")
lines(lowess(f_p_N,Risk.inf$Shed.Col), lwd = 4, col = "#de712b")
plot(h_m_N,Risk.inf$Shed.Col, log="y", xlab="No. of hand-mouth contacts", ylab="CDI Risk")
lines(lowess(h_m_N,Risk.inf$Shed.Col), lwd = 4, col = "#de712b")

#Shedding -infected
par(mfrow = c(1,4), cex=1, oma=c(0,0,0,0))
plot(p_conc_inf,Risk.inf$Shed.Inf, log="y", xlab="Infected Shedding Concentration", ylab="CDI Risk")
lines(lowess(p_conc_inf,Risk.inf$Shed.Inf), lwd = 4, col = "#de712b")
plot(f_decay,Risk.inf$Shed.Inf, log="y", xlab="Fomite decay", ylab="CDI Risk")
lines(lowess(f_decay,Risk.inf$Shed.Inf), lwd = 4, col = "#de712b")
plot(f_t,Risk.inf$Shed.Inf, log="y", xlab="Time to fomite contact", ylab="CDI Risk")
lines(lowess(f_t,Risk.inf$Shed.Inf), lwd = 4, col = "#de712b")
plot(f_h_effc,Risk.inf$Shed.Inf, log="y", xlab="Fomite-to-hand transfer efficiency", ylab="CDI Risk")
lines(lowess(f_h_effc,Risk.inf$Shed.Inf), lwd = 4, col = "#de712b")

par(mfrow = c(1,3), cex=1, oma=c(0,0,0,0))
plot(f_h_N,Risk.inf$Shed.Inf, log="y", xlab="No. of fomite contacts", ylab="CDI Risk")
lines(lowess(f_h_N,Risk.inf$Shed.Inf), lwd = 4, col = "#de712b")
plot(hcw_f_effc,Risk.inf$Shed.Inf, log="y", xlab="Hand-to-fomite transfer efficiency", ylab="CDI Risk")
lines(lowess(hcw_f_effc,Risk.inf$Shed.Inf), lwd = 4, col = "#de712b")
plot(h_area,Risk.inf$Shed.Inf, log="y", xlab="Fingertip surface area", ylab="CDI Risk")
lines(lowess(h_area,Risk.inf$Shed.Inf), lwd = 4, col = "#de712b")

#Pre-seeded fomite
par(mfrow = c(1,4), cex=1, oma=c(0,0,0,0))
plot(f_dirty_conc_env,Risk.inf$Env, log="y", xlab="Fomite Concentration", ylab="CDI Risk")
lines(lowess(f_dirty_conc_env,Risk.inf$Env), lwd = 4, col = "#de712b")
plot(f_decay,Risk.inf$Env, log="y", xlab="Fomite decay", ylab="CDI Risk")
lines(lowess(f_decay,Risk.inf$Env), lwd = 4, col = "#de712b")
plot(f_t,Risk.inf$Env, log="y", xlab="Time to fomite contact", ylab="CDI Risk")
lines(lowess(f_t,Risk.inf$Env), lwd = 4, col = "#de712b")
plot(f_h_effc,Risk.inf$Env, log="y", xlab="Fomite-to-hand transfer efficiency", ylab="CDI Risk")
lines(lowess(f_h_effc,Risk.inf$Env), lwd = 4, col = "#de712b")

par(mfrow = c(1,3), cex=1, oma=c(0,0,0,0))
plot(f_h_N,Risk.inf$Env, log="y", xlab="No. of fomite contacts", ylab="CDI Risk")
lines(lowess(f_h_N,Risk.inf$Env), lwd = 4, col = "#de712b")
plot(hcw_f_effc,Risk.inf$Env, log="y", xlab="Hand-to-fomite transfer efficiency", ylab="CDI Risk")
lines(lowess(hcw_f_effc,Risk.inf$Env), lwd = 4, col = "#de712b")
plot(h_area,Risk.inf$Env, log="y", xlab="Fingertip surface area", ylab="CDI Risk")
lines(lowess(h_area,Risk.inf$Env), lwd = 4, col = "#de712b")

#HCW-Patient
par(mfrow = c(1,4), cex=1, oma=c(0,0,0,0))
plot(hcw_p_pick,Risk.inf$P.HCW, log="y", xlab="Pick-up rate", ylab="CDI Risk")
lines(lowess(hcw_p_pick,Risk.inf$P.HCW), lwd = 4, col = "#de712b")
plot(hcw_p_N,Risk.inf$P.HCW, log="y", xlab="No. of HCW-patient contacts", ylab="CDI Risk")
lines(lowess(hcw_p_N,Risk.inf$P.HCW), lwd = 4, col = "#de712b")
plot(hcw_f_effc,Risk.inf$P.HCW, log="y", xlab="Hand-to-fomite transfer efficiency", ylab="CDI Risk")
lines(lowess(hcw_f_effc,Risk.inf$P.HCW), lwd = 4, col = "#de712b")
plot(h_area,Risk.inf$P.HCW, log="y", xlab="Fingertip surface area", ylab="CDI Risk")
lines(lowess(h_area,Risk.inf$P.HCW), lwd = 4, col = "#de712b")

#####
#Tornado plots by source of contamination

#Shedding - colonized
Risk.Vars <- data.frame(p_conc_col, f_decay, f_h_effc, f_h_N, hcw_f_effc, h_area) 
colnames(Risk.Vars) <- c("p_col", "decay_f","eff_f_h", "N_f_h", "eff_h_f", "sa_h")
#Risk.Vars <- data.frame(p_conc_col, f_decay, f_h_effc, f_h_N, hcw_f_effc, h_area, f_p_N, h_m_N) #daily risk
#colnames(Risk.Vars) <- c("p_col", "decay_f","eff_f_h", "N_f_h","t_f", "eff_h_f", "sa_h", "N_f_p", "N_h_m") #daily risk
Sens.Risk <- matrix(NA, 1, ncol(Risk.Vars)); colnames(Sens.Risk) <- colnames(Risk.Vars)
for(j in 1:ncol(Sens.Risk))
{
  Sens.Risk[,j] <- cor(log(Risk.Vars[,j],10),log(Risk.inf$Shed.Col,10), method = "spearman", use="complete.obs")
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
Risk.Sens.PLOT

#Shedding - infected
Risk.Vars <- data.frame(p_conc_inf, f_decay, f_h_effc, f_h_N, hcw_f_effc, h_area) 
colnames(Risk.Vars) <- c("p_inf", "decay_f","eff_f_h", "N_f_h", "eff_h_f", "sa_h")
#Risk.Vars <- data.frame(p_conc_inf, f_decay, f_h_effc, f_h_N, hcw_f_effc, h_area, f_p_N, h_m_N) #daily
#colnames(Risk.Vars) <- c("p_inf", "decay_f","eff_f_h", "N_f_h", "eff_h_f", "sa_h", "N_f_p", "N_h_m") #daily
Sens.Risk <- matrix(NA, 1, ncol(Risk.Vars)); colnames(Sens.Risk) <- colnames(Risk.Vars)
for(j in 1:ncol(Sens.Risk))
{
  Sens.Risk[,j] <- cor(log(Risk.Vars[,j],10),log(Risk.inf$Shed.Inf,10), method = "spearman", use="complete.obs")
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
Risk.Sens.PLOT

#Pre-seeded fomite
Risk.Vars <- data.frame(f_dirty_conc_env, f_decay, f_h_effc, f_h_N, hcw_f_effc, h_area) 
colnames(Risk.Vars) <- c("Fomite concentration", "Decay","Fomite-hand transfer effc", "HCW-fomite contact events", 
                         "Hand-fomite transfer effc", "Fingertip surface area")
#Risk.Vars <- data.frame(f_dirty_conc_env, f_decay, f_h_effc, f_h_N, hcw_f_effc, h_area, f_p_N, h_m_N) #daily
#colnames(Risk.Vars) <- c("f_conc", "decay_f","eff_f_h", "N_f_h", "eff_h_f", "sa_h", "N_f_p", "N_h_m") #daily
Sens.Risk <- matrix(NA, 1, ncol(Risk.Vars)); colnames(Sens.Risk) <- colnames(Risk.Vars)
for(j in 1:ncol(Sens.Risk))
{
  Sens.Risk[,j] <- cor(log(Risk.Vars[,j],10),log(Risk.inf$Env,10), method = "spearman", use="complete.obs")
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
ggsave("Risk.Sens.PLOT_Seed.png",Risk.Sens.PLOT,dpi=600)

#HCW-Patient
Risk.Vars <- data.frame(hcw_p_pick, f_h_N, hcw_f_effc, h_area) 
colnames(Risk.Vars) <- c("pickup", "N_f_h","eff_h_f", "sa_h")
#Risk.Vars <- data.frame(hcw_p_pick, f_h_N, hcw_f_effc, h_area, f_p_N, h_m_N) #daily
#colnames(Risk.Vars) <- c("pickup", "N_f_h","eff_h_f", "sa_h", "N_f_p", "N_h_m") #daily
Sens.Risk <- matrix(NA, 1, ncol(Risk.Vars)); colnames(Sens.Risk) <- colnames(Risk.Vars)
for(j in 1:ncol(Sens.Risk))
{
  Sens.Risk[,j] <- cor(log(Risk.Vars[,j],10),log(Risk.inf$Env,10), method = "spearman", use="complete.obs")
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
Risk.Sens.PLOT


#####
#Alternative parameterization

#Contamination source
p_conc_treat <- vector()
f_dirty_conc_treat <- vector()

p_conc_treat <- rtriang(iters,min=0,max=4.4,mean=2.4) #shedding patient under treatment
f_dirty_conc_treat <- (10^p_conc_treat)*shed_effc

f_dirty_conc_env2 <- rnormTrunc(iters, 1.33, 4.69, min=0) #alt pre-seeded fomite

hcw_p_pick[i] <- rnormTrunc(iters, 0.247, 0.0429, min=0) #alt HCW pick up rate

contam <- c("Shed-Treat", "Env2") 
contam2 <- c("Shed-Treat", "Env2", "P-HCW2")
N.contam <- length(contam)
N.contam2 <- length(contam2)
contam.env <- data.frame(f_dirty_conc_treat, f_dirty_conc_env2)
colnames(contam.env) <- contam

summary(contam.env)
summary(log10(contam.env))

#Persistence
#Growth on fomite (replace decay)
k1 <- -10.09
k2 <- 0.9572

f_growth[i] <- decay_epd(k1,k2,f_t[i])


#END
