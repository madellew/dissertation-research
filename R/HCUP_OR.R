#Calculate odds ratios by length of hospital stay and device utilization

#Created: 12 Nov 2024
#Last modified: 13 Nov 2024

#####
#Length of stay

#CDI
logit_CDI <- glm(CDI_HA ~ LOS.num, data=data, family = "binomial")
summary(logit_CDI)
exp(coef(logit_CDI)) 
exp(confint(logit_CDI))

ggpredict(logit_CDI, terms=c("LOS.num[all]")) |>
  plot(colors="#5D69B1FF") +
  theme_bw() +
  labs(title=NULL, x="Length of stay", y="Predicted probability")
ggsave("LOS_CDI_predict.png", dpi=300, height=6, width=7, units="in")

#CLABSI
logit_CLABSI <- glm(CLABSI_HA ~ LOS.num, data=data, family = "binomial")
summary(logit_CLABSI)
exp(coef(logit_CLABSI)) 
exp(confint(logit_CLABSI))

ggpredict(logit_CLABSI, terms=c("LOS.num[all]")) |>
  plot(colors="#5D69B1FF") +
  theme_bw() +
  labs(title=NULL, x="Length of stay", y="Predicted probability")
ggsave("LOS_CLABSI_predict.png", dpi=300, height=6, width=7, units="in")

#CAUTI
logit_CAUTI <- glm(CAUTI_HA ~ LOS.num, data=data, family = "binomial")
summary(logit_CAUTI)
exp(coef(logit_CAUTI)) 
exp(confint(logit_CAUTI))

ggpredict(logit_CAUTI, terms=c("LOS.num[all]")) |>
  plot(colors="#5D69B1FF") +
  theme_bw() +
  labs(title=NULL, x="Length of stay", y="Predicted probability")
ggsave("LOS_CAUTI_predict.png", dpi=300, height=6, width=7, units="in")

#Candidiasis
logit_CAN <- glm(CAN_HA ~ LOS.num, data=data, family = "binomial")
summary(logit_CAN)
exp(coef(logit_CAN)) 
exp(confint(logit_CAN))

ggpredict(logit_CAN, terms=c("LOS.num[all]")) |>
  plot(colors="#5D69B1FF") +
  theme_bw() +
  labs(title=NULL, x="Length of stay", y="Predicted probability")
ggsave("LOS_CAN_predict.png", dpi=300, height=6, width=7, units="in")

#####
#Devices

#central venous catheter
logit_CVC <- glm(CLABSI_HA ~ CVC1, data=data, family = "binomial")
summary(logit_CVC)
exp(coef(logit_CVC)) 
exp(confint(logit_CVC))

logit_CVC_los <- glm(CLABSI_HA ~ CVC1 + LOS.num, data=data, family = "binomial")
summary(logit_CVC_los)
exp(coef(logit_CVC_los)) 
exp(confint(logit_CVC_los))

dat <- predict_response(logit_CVC_los, terms=c("LOS.num","CVC1")) #plotting option #1

ggpredict(logit_CVC_los, terms=c("LOS.num[all]","CVC1")) |>  #plotting option #2
  plot() +
  theme_bw() +
  labs(title=NULL, x="Length of stay", y="Predicted probability", colour="CVC", fill="CVC") +
  scale_color_carto_d(palette="Vivid", labels=c("No","Yes")) +
  scale_fill_carto_d(palette="Vivid", labels=c("No","Yes")) +
  guides(color=guide_legend(override.aes = list(fill=NA)))
ggsave("LOS_CVC_predict.png", dpi=300, height=6, width=7, units="in")

#urinary catheter
logit_URCATH <- glm(CAUTI ~ URCATH1, data=data, family = "binomial")
summary(logit_URCATH)
exp(coef(logit_URCATH)) 
exp(confint(logit_URCATH))

logit_URCATH_los <- glm(CAUTI ~ URCATH1 + LOS, data=data, family = "binomial")
summary(logit_URCATH_los)
exp(coef(logit_URCATH_los)) 
exp(confint(logit_URCATH_los))

#END
