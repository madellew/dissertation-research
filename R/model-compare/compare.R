#Compare mean and median estimates from EITS and QMRA models, HCUP NIS 2021 data
#analysis, and HCUP SID analysis by Miller et al

#Created: 12 Nov 2024
#Last modified: 20 Jan 2025

#Load estimates
mean_est <- read.csv("mean.csv")
med_est <- read.csv("median.csv")

#Plot mean estimates

plt_Mean <- ggplot(mean_est, aes(Outcome, INC, group=Outcome, shape=Method, color=Method)) +
  geom_jitter(alpha=.75 ,size=4, width=0.25) +
  scale_y_log10() +
  scale_shape_manual(values=c(18, 15, 17, 16)) +
  theme_bw() +
  theme(text=element_text(family="serif", size=13), axis.title.y = element_text(margin = margin(r = 20))) +
  labs(x="Outcome", y="Mean Cumulative Incidence Estimates") +
  scale_color_paletteer_d("calecopal::superbloom3")
ragg::agg_png("EST-mean.png", res=500, height=4, width=5, units="in", scaling=0.6)
plt_Mean
dev.off()

#Plot median estimates

plot_Med <- ggplot(med_est, aes(x=Outcome, y=INC, group=Outcome, shape=Method, color=Method)) +
  geom_jitter(alpha=.75, size=4, width=0.25) +
  scale_y_log10() +
  scale_shape_manual(values=c(18, 15, 17, 16)) +
  theme_bw() +
  theme(text=element_text(family="serif", size=13), axis.title.y = element_text(margin = margin(r = 20))) +
  labs(x="Outcome", y="Median Cumulative Incidence Estimates") +
  scale_color_paletteer_d("calecopal::superbloom3")
ragg::agg_png("EST-med.png", res=500, height=4, width=5, units="in", scaling=0.6)
plot_Med
dev.off()
