
#QMRA model for Candida auris central line-associated bloodstream infections - Sensitivity analysis

#Created: 1 Aug 2024
#Last modified: 21 Jan 2025

#Tornado plots by alt source of contamination

#Shedding Low
Risk.Vars <- data.frame(p_col_conc1, t_shed, f_decay, f_h_effc, f_h_N, f_t, h_c_effc, h_decay, h_t, k_bsi) 
colnames(Risk.Vars) <- c("c_col", "t_shed", "decay_f","eff_f_h", "N_f_h","t_f", "eff_h_c", "decay_h", "t_h", "k")
Sens.Risk <- matrix(NA, 1, ncol(Risk.Vars)); colnames(Sens.Risk) <- colnames(Risk.Vars)
for(j in 1:ncol(Sens.Risk))
{
  Sens.Risk[,j] <- cor(log(Risk.Vars[,j],10),log(Risk.contam.clabsi$`Shed Low`,10), method = "spearman", use="complete.obs")
}
Risk.Sens <- melt(Sens.Risk); colnames(Risk.Sens) <- c("dummy","Variable","Rho")
Risk.Sens.PLOT <- ggplot(Risk.Sens,aes(x=Variable,y=Rho, fill=Variable)) + 
  geom_bar(stat = "identity", color = "black", fill = "grey") + 
  coord_flip() + 
  theme_bw() +
  theme(text=element_text(family="serif",size=13), legend.position = "none") + 
  labs(y = expression(~Spearman~rho), x = "Model Variable", title = " ")
Risk.Sens.PLOT

#Shedding High
Risk.Vars <- data.frame(p_col_conc2, t_shed, f_decay, f_h_effc, f_h_N, f_t, h_c_effc, h_decay, h_t, k_bsi) 
colnames(Risk.Vars) <- c("c_col", "t_shed", "decay_f","eff_f_h", "N_f_h","t_f", "eff_h_c", "decay_h", "t_h", "k")
Sens.Risk <- matrix(NA, 1, ncol(Risk.Vars)); colnames(Sens.Risk) <- colnames(Risk.Vars)
for(j in 1:ncol(Sens.Risk))
{
  Sens.Risk[,j] <- cor(log(Risk.Vars[,j],10),log(Risk.contam.clabsi$`Shed High`,10), method = "spearman", use="complete.obs")
}
Risk.Sens <- melt(Sens.Risk); colnames(Risk.Sens) <- c("dummy","Variable","Rho")
Risk.Sens.PLOT <- ggplot(Risk.Sens,aes(x=Variable,y=Rho, fill=Variable)) + 
  geom_bar(stat = "identity", color = "black", fill = "grey") + 
  coord_flip() + 
  theme_bw() +
  theme(text=element_text(family="serif",size=13), legend.position = "none") + 
  labs(y = expression(~Spearman~rho), x = "Model Variable", title = " ")
Risk.Sens.PLOT

#Pre-Seeded
Risk.Vars <- data.frame(f_dirty_conc_all, f_decay, f_h_effc, f_h_N, f_t, h_c_effc, h_decay, h_t, k_bsi) 
colnames(Risk.Vars) <- c("c_f", "decay_f","eff_f_h", "N_f_h","t_f", "eff_h_c", "decay_h", "t_h", "k")
Sens.Risk <- matrix(NA, 1, ncol(Risk.Vars)); colnames(Sens.Risk) <- colnames(Risk.Vars)
for(j in 1:ncol(Sens.Risk))
{
  Sens.Risk[,j] <- cor(log(Risk.Vars[,j],10),log(Risk.contam.clabsi$All,10), method = "spearman", use="complete.obs")
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
