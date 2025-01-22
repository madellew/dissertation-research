#Outcomes of interest by admission month and hospital subset analysis

#Created: 12 Nov 2024
#Last modified: 20 Jan 2025

#####
#Data subset - admission month

data_hosp <- data %>%
  convert(num(CDI_HA,CAN,CAN_HA,CLABSI_HA,CAUTI_HA,CVC,CVC1,URCATH,URCATH1))

amonth <- data_hosp %>%
  group_by(AMONTH) %>%
  summarise(count = n(), CDI = sum(CDI_HA, na.rm=T), CAUTI = sum(CAUTI_HA, na.rm=T),
            CLABSI = sum(CLABSI_HA, na.rm=T), CAN = sum(CAN_HA, na.rm=T))
amonth.plot <- amonth[,-2]
amonth.plot <- amonth.plot %>% slice(-13)
amonth.plot <- melt(amonth.plot)
colnames(amonth.plot) <- c("AMONTH", "OUTCOME", "COUNT")

#Discharge record counts ("cases") for outcomes of interest by admission month
plt4 <- ggplot(amonth.plot, aes(x=AMONTH, y=COUNT, group=OUTCOME, fill=OUTCOME)) +
  geom_bar(stat="identity", position=position_dodge(), alpha=.75) +
  theme_bw() +
  labs(x="Admission month", y="Cases", fill="Outcome") +
  scale_fill_paletteer_d("calecopal::superbloom3", labels=c("CDI","CAUTI", "CLABSI", "Candidiasis"))
ragg::agg_png("Cases_AMONTH.png", res=500, height=3.5, width=6, units="in", scaling=0.75)
plt4
dev.off()

#Discharge record counts over total record counts ("incidence") by admission month
amonth.plot <- amonth %>%
  slice(-13) %>%
  group_by(AMONTH) %>%
  mutate(CDI=CDI/count, CAUTI=CAUTI/count, CLABSI=CLABSI/count, CAN=CAN/count)
amonth.plot <- amonth.plot[,-2]
amonth.plot <- melt(amonth.plot)
colnames(amonth.plot) <- c("AMONTH", "OUTCOME", "INC")

plt5 <- ggplot(amonth.plot, aes(x=AMONTH, y=INC, group=OUTCOME, color=OUTCOME)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(text=element_text(family="serif", size=13)) +
  labs(x="Admission month", y="Cumulative Incidence", color="Outcome") +
  scale_color_paletteer_d("calecopal::superbloom3", labels=c("CDI","CAUTI", "CLABSI", "Candidiasis"))
ragg::agg_png("Cases_AMONTH-INC.png", res=500, height=3.5, width=6, units="in", scaling=0.6)
plt5
dev.off()

#####
#Data subset - hospitals

hosp_inf <- data_hosp %>%
  group_by(HOSP_NIS) %>%
  summarise(count = n(), CDI_HOSP = sum(CDI_HA, na.rm=T), 
            CAUTI_HOSP = sum(CAUTI_HA, na.rm=T), CLABSI_HOSP = sum(CLABSI_HA, na.rm=T), 
            CAN_HOSP = sum(CAN_HA, na.rm=T))

#Hospitals with at least one record of outcome of interest   
subset(hosp_inf, CDI_HOSP >= 1) %>%
Hmisc::describe()

##Discharge record counts ("cases") for outcomes of interest by NIS hospital
hosp.inf.plot <- hosp_inf[,-2]
colnames(hosp.inf.plot) <- c("HOSP_NIS","CDI","CAUTI", "CLABSI", "Candidiasis")
hosp.inf.plot <- melt(hosp.inf.plot, id.vars = "HOSP_NIS")
colnames(hosp.inf.plot) <- c("HOSP_NIS", "OUTCOME", "COUNT")

plt3 <- ggplot(hosp.inf.plot, aes(OUTCOME, COUNT)) +
  geom_jitter(aes(color=OUTCOME),alpha=.5) +
  theme_bw() +
  labs(x="Outcome", y="Cases", color="Outcome") +
  scale_color_paletteer_d("calecopal::superbloom3", labels=c("CDI","CAUTI", "CLABSI", "Candidiasis"))
ragg::agg_png("Cases_HOSP.png", res=500, height=3.5, width=6, units="in", scaling=0.75)
plt3
dev.off()

#Discharge record counts over total record counts ("incidence") by NIS hospital
hosp.inf.plot2 <- hosp_inf %>%
  group_by(HOSP_NIS) %>%
  mutate(CDI=CDI_HOSP/count, CAUTI=CAUTI_HOSP/count, CLABSI=CLABSI_HOSP/count, CAN=CAN_HOSP/count,count=count)

hosp.inf.plot2 <- hosp.inf.plot2[,-c(2:6)]
colnames(hosp.inf.plot2) <- c("HOSP_NIS","CDI","CAUTI", "CLABSI", "Candidiasis")
hosp.inf.plot2 <- melt(hosp.inf.plot2, id.vars = "HOSP_NIS")
colnames(hosp.inf.plot2) <- c("HOSP_NIS", "OUTCOME", "INC")

plt2 <- ggplot(hosp.inf.plot2, aes(OUTCOME, INC)) +
  geom_jitter(aes(color=OUTCOME),alpha=.5) +
  scale_y_log10() +
  theme_bw() +
  theme(text=element_text(family="serif", size=13)) +
  labs(x="Outcome", y="Cumulative Incidence", color="Outcome") +
  scale_color_paletteer_d("calecopal::superbloom3", labels=c("CDI","CAUTI", "CLABSI", "Candidiasis"))
ragg::agg_png("Cases_HOSP-INC.png", res=500, height=3.5, width=6, units="in", scaling=0.6)
plt2
dev.off()

#plot without zeros
hosp.inf.plot2[hosp.inf.plot2 == 0] <- NA

plt <- ggplot(hosp.inf.plot2, aes(OUTCOME, INC)) +
  geom_jitter(aes(color=OUTCOME),alpha=.5) +
  scale_y_log10() +
  theme_bw() +
  theme(text=element_text(family="serif", size=13)) +
  labs(x="Outcome", y="Cumulative Incidence", color="Outcome") +
  scale_color_paletteer_d("calecopal::superbloom3", labels=c("CDI","CAUTI", "CLABSI", "Candidiasis"))
ragg::agg_png("Cases_HOSP-INC-log10-no-zeroes.png", res=500, height=3.5, width=6, units="in", scaling=0.6)
plt
dev.off()

#####
#Subset - specific hospitals (hospitals with highest counts or incidence for each outcome)

#####
#CDI
CDI_HOSP <- subset(data_hosp, HOSP_NIS == 30496 | HOSP_NIS == 30709 | HOSP_NIS == 80175, 
                   select=c(LOS,CDI_HA,AMONTH,TOTAL_DISC,HOSP_NIS))

amonth_HOSP <- CDI_HOSP %>%
  group_by(AMONTH, HOSP_NIS) %>%
  summarise(CDI = sum(CDI_HA, na.rm=T)) %>%
  convert(fct(HOSP_NIS))

plt.CDI <- ggplot(amonth_HOSP, aes(x=AMONTH, y=CDI, group=HOSP_NIS, fill=HOSP_NIS)) +
  geom_bar(stat="identity", position=position_dodge(), alpha=.75) +
  theme_bw() +
  labs(x="Admission month", y="Cases", fill="Hospital") +
  scale_fill_paletteer_d("tvthemes::Alexandrite", labels=c("Hospital 1","Hospital 2", "Hospital 3"))
ragg::agg_png("CDI_HOSP.png", res=500, height=3.5, width=6, units="in", scaling=0.5)
plt.CDI
dev.off()

amonth_HOSP_INC <- CDI_HOSP %>%
  group_by(AMONTH, HOSP_NIS) %>%
  summarise(count = n(), CDI = sum(CDI_HA, na.rm=T)) %>%
  mutate(CDI_INC = CDI/count) %>%
  convert(fct(HOSP_NIS))
amonth_HOSP_INC <- amonth_HOSP_INC[-c(3,4)]

plt.CDI2 <- ggplot(amonth_HOSP_INC, aes(x=AMONTH, y=CDI_INC, group=HOSP_NIS, color=HOSP_NIS)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x="Admission month", y="Cumulative Incidence", color="Hospital") +
  scale_color_paletteer_d("tvthemes::Alexandrite", labels=c("Hospital 1","Hospital 2", "Hospital 3"))
ragg::agg_png("CDI_HOSP_INC.png", res=500, height=3.5, width=6, units="in", scaling=0.5)
plt.CDI2
dev.off()

#####
#CAUTI
CAUTI_HOSP <- subset(data_hosp, HOSP_NIS == 70082 | HOSP_NIS == 50665 | HOSP_NIS == 30443, 
                   select=c(LOS,CAUTI_HA,AMONTH,TOTAL_DISC,HOSP_NIS))

amonth_HOSP <- CAUTI_HOSP %>%
  group_by(AMONTH, HOSP_NIS) %>%
  summarise(CAUTI = sum(CAUTI_HA, na.rm=T)) %>%
  convert(fct(HOSP_NIS))

plt.CAUTI <- ggplot(amonth_HOSP, aes(x=AMONTH, y=CAUTI, group=HOSP_NIS, fill=HOSP_NIS)) +
  geom_bar(stat="identity", position=position_dodge(), alpha=.75) +
  theme_bw() +
  labs(x="Admission month", y="Cases", fill="Hospital") +
  scale_fill_paletteer_d("tvthemes::Alexandrite", labels=c("Hospital 1","Hospital 2", "Hospital 3"))
ragg::agg_png("CAUTI_HOSP.png", res=500, height=3.5, width=6, units="in", scaling=0.5)
plt.CAUTI
dev.off()

amonth_HOSP_INC <- CAUTI_HOSP %>%
  group_by(AMONTH, HOSP_NIS) %>%
  summarise(count = n(), CAUTI = sum(CAUTI_HA, na.rm=T)) %>%
  mutate(CAUTI_INC = CAUTI/count) %>%
  convert(fct(HOSP_NIS))
amonth_HOSP_INC <- amonth_HOSP_INC[-c(3,4)]

plt.CAUTI2 <- ggplot(amonth_HOSP_INC, aes(x=AMONTH, y=CAUTI_INC, group=HOSP_NIS, color=HOSP_NIS)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x="Admission month", y="Cumulative Incidence", color="Hospital") +
  scale_color_paletteer_d("tvthemes::Alexandrite", labels=c("Hospital 1","Hospital 2", "Hospital 3"))
ragg::agg_png("CAUTI_HOSP_INC.png", res=500, height=3.5, width=6, units="in", scaling=0.5)
plt.CAUTI2
dev.off()

#####
#CLABSI
CLABSI_HOSP <- subset(data_hosp, HOSP_NIS == 70743 | HOSP_NIS == 40472 | HOSP_NIS == 90006, 
                     select=c(LOS,CLABSI_HA,AMONTH,TOTAL_DISC,HOSP_NIS))

amonth_HOSP <- CLABSI_HOSP %>%
  group_by(AMONTH, HOSP_NIS) %>%
  summarise(CLABSI = sum(CLABSI_HA, na.rm=T)) %>%
  convert(fct(HOSP_NIS))

plt.CLABSI <- ggplot(amonth_HOSP, aes(x=AMONTH, y=CLABSI, group=HOSP_NIS, fill=HOSP_NIS)) +
  geom_bar(stat="identity", position=position_dodge(), alpha=.75) +
  theme_bw() +
  labs(x="Admission month", y="Cases", fill="Hospital") +
  scale_fill_paletteer_d("tvthemes::Alexandrite", labels=c("Hospital 1","Hospital 2", "Hospital 3"))
ragg::agg_png("CLABSI_HOSP.png", res=500, height=3.5, width=6, units="in", scaling=0.5)
plt.CLABSI
dev.off()

amonth_HOSP_INC <- CLABSI_HOSP %>%
  group_by(AMONTH, HOSP_NIS) %>%
  summarise(count = n(), CLABSI = sum(CLABSI_HA, na.rm=T)) %>%
  mutate(CLABSI_INC = CLABSI/count) %>%
  convert(fct(HOSP_NIS))
amonth_HOSP_INC <- amonth_HOSP_INC[-c(3,4)]

plt.CLABSI2 <- ggplot(amonth_HOSP_INC, aes(x=AMONTH, y=CLABSI_INC, group=HOSP_NIS, color=HOSP_NIS)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x="Admission month", y="Cumulative Incidence", color="Hospital") +
  scale_color_paletteer_d("tvthemes::Alexandrite", labels=c("Hospital 1","Hospital 2", "Hospital 3"))
ragg::agg_png("CLABSI_HOSP_INC.png", res=500, height=3.5, width=6, units="in", scaling=0.5)
plt.CLABSI2
dev.off()

#####
#CAN
CAN_HOSP <- subset(data_hosp, HOSP_NIS == 70079 | HOSP_NIS == 90048 | HOSP_NIS == 50088 | HOSP_NIS == 70540,
                      select=c(LOS,CAN_HA,AMONTH,TOTAL_DISC,HOSP_NIS))

amonth_HOSP <- CAN_HOSP %>%
  group_by(AMONTH, HOSP_NIS) %>%
  summarise(CAN = sum(CAN_HA, na.rm=T)) %>%
  convert(fct(HOSP_NIS))

plt.CAN <- ggplot(amonth_HOSP, aes(x=AMONTH, y=CAN, group=HOSP_NIS, fill=HOSP_NIS)) +
  geom_bar(stat="identity", position=position_dodge(), alpha=.75) +
  theme_bw() +
  labs(x="Admission month", y="Cases", fill="Hospital") +
  scale_fill_paletteer_d("tvthemes::Alexandrite", labels=c("Hospital 1","Hospital 2", "Hospital 3", "Hospital 4"))
ragg::agg_png("CAN_HOSP.png", res=500, height=3.5, width=6, units="in", scaling=0.5)
plt.CAN
dev.off()

amonth_HOSP_INC <- CAN_HOSP %>%
  group_by(AMONTH, HOSP_NIS) %>%
  summarise(count = n(), CAN = sum(CAN_HA, na.rm=T)) %>%
  mutate(CAN_INC = CAN/count) %>%
  convert(fct(HOSP_NIS))
amonth_HOSP_INC <- amonth_HOSP_INC[-c(3,4)]

plt.CAN2 <- ggplot(amonth_HOSP_INC, aes(x=AMONTH, y=CAN_INC, group=HOSP_NIS, color=HOSP_NIS)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x="Admission month", y="Cumulative Incidence", color="Hospital") +
  scale_color_paletteer_d("tvthemes::Alexandrite", labels=c("Hospital 1","Hospital 2", "Hospital 3", "Hospital 4"))
ragg::agg_png("CAN_HOSP_INC.png", res=500, height=3.5, width=6, units="in", scaling=0.5)
plt.CAN2
dev.off()

#END

