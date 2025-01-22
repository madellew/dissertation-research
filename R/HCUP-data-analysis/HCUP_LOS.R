#Analysis of length of stay by outcome variable and device

#Created: 12 Nov 2024
#Last modified: 21 Jan 2025

#LOS summarized by variables
describeBy(data$LOS, data$CDI_HA)
describeBy(data$LOS, data$CAN_HA)
describeBy(data$LOS, data$CAN_UTI_np)
describeBy(data$LOS, data$CAUTI_HA)
describeBy(data$LOS, data$CLABSI_HA)
describeBy(data$LOS, data$URCATH1)
describeBy(data$LOS, data$CVC1)

wilcox.test(LOS~CDI_HA, data=data)
wilcox.test(LOS~CAN_HA, data=data)
wilcox.test(LOS~CAUTI_HA, data=data)
wilcox.test(LOS~CLABSI_HA, data=data)
wilcox.test(LOS~CVC1, data=data)
wilcox.test(LOS~URCATH1, data=data)

#Estimate 95th and 99th percentiles for LOS by outcome

data %>%
  group_by(CLABSI_HA) %>%
  dplyr::summarize(q95 = quantile(LOS, .95, na.rm = TRUE),
                   q99 = quantile(LOS, .99, na.rm = TRUE))

#####
#Plot LOS by outcome (compare with and without)

#CDI
los.CDI <- ggplot(data, aes(x=LOS, fill=CDI_HA, colour=CDI_HA)) +
  geom_density(position="identity", alpha=0.25) +
  scale_y_sqrt() +
  theme_bw() +
  theme(text=element_text(family="serif", size=13)) +
  #xlim(0,77) + #truncation at 99th percentile
  labs(x="Length of hospital stay (days)", y="Density", fill="CDI", colour="CDI") +
  scale_color_paletteer_d("rcartocolor::Bold", labels=c("No","Yes")) +
  scale_fill_paletteer_d("rcartocolor::Bold", labels=c("No","Yes"))
ggsave("LOS_CDI_HA.png", los.CDI, dpi=600, height=6, width=6, units="in")

#CAUTI
los.CAUTI <- ggplot(data, aes(x=LOS, fill=CAUTI_HA, colour=CAUTI_HA)) +
  geom_density(position="identity", alpha=0.25) +
  scale_y_sqrt() +
  theme_bw() +
  theme(text=element_text(family="serif", size=13)) +
  xlim(0,86) + #truncation at 99th percentile
  labs(x="Length of hospital stay (days)", y="Density", fill="CAUTI", colour="CAUTI") +
  scale_color_paletteer_d("rcartocolor::Bold", labels=c("No","Yes")) +
  scale_fill_paletteer_d("rcartocolor::Bold", labels=c("No","Yes"))
ggsave("LOS_CAUTI_HA.png", los.CAUTI, dpi=600, height=6, width=6, units="in")

#CLABSI
los.CLABSI <- ggplot(data, aes(x=LOS, fill=CLABSI_HA, colour=CLABSI_HA)) +
  geom_density(position="identity", alpha=0.25) +
  scale_y_sqrt() +
  theme_bw() +
  theme(text=element_text(family="serif", size=13)) +
  xlim(0,205) +  #truncation at 99th percentile
  labs(x="Length of hospital stay (days)", y="Density", fill="CLABSI", colour="CLABSI") +
  scale_color_paletteer_d("rcartocolor::Bold", labels=c("No","Yes")) +
  scale_fill_paletteer_d("rcartocolor::Bold", labels=c("No","Yes"))
ggsave("LOS_CLABSI_HA.png", los.CLABSI, dpi=600, height=6, width=6, units="in")

#Candidiasis
los.CAN <- ggplot(data, aes(x=LOS, fill=CAN_HA, colour=CAN_HA)) +
  geom_density(position="identity", alpha=0.25) +
  scale_y_sqrt() +
  theme_bw() +
  theme(text=element_text(family="serif", size=13)) +
  xlim(0,127) +  #truncation at 99th percentile
  labs(x="Length of hospital stay (days)", y="Density", fill="CAN", colour="CAN") +
  scale_color_paletteer_d("rcartocolor::Bold", labels=c("No","Yes")) +
  scale_fill_paletteer_d("rcartocolor::Bold", labels=c("No","Yes"))
ggsave("LOS_CAN_HA.png", los.CAN, dpi=600, height=6, width=6, units="in")

#####
#Calculate person-time
data %>% summarise(LOS_tot = sum(LOS.num, na.rm=T))
data %>%
  group_by(CLABSI_HA) %>%
  summarise(LOS_CLABSI = sum(LOS.num, na.rm=T))
data %>%
  group_by(CDI_HA) %>%
  summarise(LOS_CDI = sum(LOS.num, na.rm=T))
data %>%
  group_by(CAUTI_HA) %>%
  summarise(LOS_CAUTI = sum(LOS.num, na.rm=T))
data %>%
  group_by(CAN_HA) %>%
  summarise(LOS_CAN = sum(LOS.num, na.rm=T))


#END
