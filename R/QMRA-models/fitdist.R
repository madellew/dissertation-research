#Fitting distributions to Candida auris concentration data from Zhu et al 2020

#Created: 9 Aug 2024
#Last modified: 14 Sep 2024

#####
#Require necessary packages
library(tidyverse)
library(dplyr)
library(fitdistrplus)
library(VGAM)
library(mixtools)
library(mclust)

#Import data

##Patient skin colonization concentration data - C. auris
zhu_axilla <- read.csv("zhu_axilla_conc.csv")
zhu_nares <- read.csv("zhu_nares_conc.csv")

##Environmental colonization concentration data - C. auris
zhu_non_porous <- read.csv("zhu_non-porous_conc.csv")
zhu_porous <- read.csv("zhu_porous_conc.csv")

#####
#Standardize concentration data 

#Transform (imported data is log10)
zhu_axilla <- zhu_axilla %>%
  mutate(conc_swab = 10^y_swab)
zhu_nares <- zhu_nares %>%
  mutate(conc_swab = 10^y_swab)
zhu_non_porous <- zhu_non_porous %>%
  mutate(conc_surface = 10^y_surface)
zhu_porous <- zhu_porous %>%
  mutate(conc_surface <- 10^y)

#Transform units from CFU/swab or surface to CFU/cm^2
zhu_axilla <- zhu_axilla %>%
  mutate(conc = (conc_swab/100)) 
zhu_nares <- zhu_nares %>%
  mutate(conc = (conc_swab/150))
zhu_non_porous <- zhu_non_porous %>%
  mutate(conc = (conc_surface/100))

#Now back to log-transformed#
zhu_axilla <- zhu_axilla %>%
  mutate(y = log10(conc))
zhu_nares <- zhu_nares %>%
  mutate(y = log10(conc))
zhu_non_porous <- zhu_non_porous %>%
  mutate(y = log10(conc))

#####
#First look at data

plotdist(zhu_axilla$conc, histo = TRUE, demp = TRUE)
plotdist(zhu_nares$conc, histo = TRUE, demp = TRUE)
plotdist(zhu_non_porous$conc, histo = TRUE, demp = TRUE)
plotdist(zhu_porous$conc, histo = TRUE, demp = TRUE)
#And log transformed data
plotdist(zhu_axilla$y, histo = TRUE, demp = TRUE)
plotdist(zhu_nares$y, histo = TRUE, demp = TRUE)
plotdist(zhu_non_porous$y, histo = TRUE, demp = TRUE)
plotdist(zhu_porous$y, histo = TRUE, demp = TRUE)

#Combine some datasets

#wilcox.test(zhu_axilla$conc, zhu_nares$conc)
zhu_skin <- rbind(zhu_axilla, zhu_nares)
zhu_env <- rbind(zhu_non_porous, zhu_porous)

#Look at combined data

plotdist(zhu_skin$conc, histo = TRUE, demp = TRUE)
plotdist(zhu_env$conc, histo = TRUE, demp = TRUE)
##log10
plotdist(zhu_skin$y, histo = TRUE, demp = TRUE)
plotdist(zhu_env$y, histo = TRUE, demp = TRUE)

#####
#Calculate skewness and kurtosis to investigate which distributions to use

descdist(zhu_axilla$conc, boot = 1000)
descdist(zhu_nares$conc, boot = 1000)
descdist(zhu_skin$conc, boot = 1000)
#log10
descdist(zhu_axilla$y, boot = 1000, discrete=FALSE)
descdist(zhu_nares$y, boot = 1000)
descdist(zhu_skin$y, boot = 1000)

descdist(zhu_non_porous$conc, boot = 1000)
descdist(zhu_porous$conc, boot = 1000)
descdist(zhu_env$conc, boot = 1000)
#log10
descdist(zhu_non_porous$y, boot = 1000)
descdist(zhu_porous$y, boot = 1000)
descdist(zhu_env$y, boot = 1000)

#####
#Fit data to distributions

fit_func <- function(x, fits2try){
  #x is a list of fitdist objects
  #fits2try is a list of names
  cdfcomp(x, legendtext=fits2try)
  denscomp(x, legendtext=fits2try)
  qqcomp(x, legendtext=fits2try)
  ppcomp(x, legendtext=fits2try)
  gofstat(x, fitnames=fits2try)
}


#Zhu Axilla Data ####

fit_lnorm_aix <- fitdist(zhu_axilla$conc, "lnorm")
#fit_gam_aix <- fitdist(zhu_axilla$conc, "gamma") 
  #mle failed to estimate the parameters
#fit_wei_aix <- fitdist(zhu_axilla$conc, "weibull") 
  #NaNs produced

summary(fit_lnorm_aix)

#graph
par(mfrow = c(2, 2)) 
denscomp(fit_lnorm_aix) 
qqcomp(fit_lnorm_aix) 
cdfcomp(fit_lnorm_aix) 
ppcomp(fit_lnorm_aix)

#Using log-transformed data
fit_log10norm_aix <- fitdist(zhu_axilla$y, "norm")
fit_log10unif_aix <- fitdist(zhu_axilla$y, "unif")
fit_log10wei_aix <- fitdist(zhu_axilla$y, "weibull") 

summary(fit_log10norm_aix)
summary(fit_log10unif_aix)
summary(fit_log10wei_aix)

fits_compare_aix <- list(fit_log10norm_aix, fit_log10gam_aix, fit_log10wei_aix)
fit_names_aix <- c("normal", "gamma", "weibull")
fit_func(fits_compare_aix, fit_names_aix)

#Zhu Nares Data ####

fit_lnorm_nare <- fitdist(zhu_nares$conc, "lnorm")
summary(fit_lnorm_nare)

#graph
par(mfrow = c(2, 2)) 
denscomp(fit_log10unif_nare) 
qqcomp(fit_log10unif_nare) 
cdfcomp(fit_log10unif_nare) 
ppcomp(fit_log10unif_nare)

#Using log-transformed data
fit_log10norm_nare <- fitdist(zhu_nares$y, "norm")
fit_log10unif_nare <- fitdist(zhu_nares$y, "unif")
fit_log10gam_nare <- fitdist(zhu_nares$y, "gamma")
fit_log10wei_nare <- fitdist(zhu_nares$y, "weibull") 

summary(fit_log10norm_nare)
summary(fit_log10gam_nare)
summary(fit_log10wei_nare)

fits_compare_nare <- list(fit_log10norm_nare, fit_log10gam_nare, fit_log10wei_nare)
fit_names_nare <- c("normal", "gamma", "weibull")
fit_func(fits_compare_nare, fit_names_nare)

#Zhu skin combined data ####

fit_lnorm_zhu_skin <- fitdist(zhu_skin$conc, "lnorm")

summary(fit_lnorm_zhu_skin)

#graph
par(mfrow = c(2, 2)) 
denscomp(fit_lnorm_zhu_skin) 
qqcomp(fit_lnorm_zhu_skin) 
cdfcomp(fit_lnorm_zhu_skin) 
ppcomp(fit_lnorm_zhu_skin)

##Using log-transformed data
fit_log10norm_zhu_skin <- fitdist(zhu_skin$y, "norm")
fit_log10gam_zhu_skin <- fitdist(zhu_skin$y, "gamma")
fit_log10wei_zhu_skin <- fitdist(zhu_skin$y, "weibull") 

summary(fit_log10norm_zhu_skin)
summary(fit_log10gam_zhu_skin)
summary(fit_log10wei_zhu_skin)

fits_compare_zhu_skin <- list(fit_log10norm_zhu_skin, fit_log10gam_zhu_skin, fit_log10wei_zhu_skin)
fit_names_zhu_skin <- c("normal", "gamma", "weibull")
fit_func(fits_compare_zhu_skin, fit_names_zhu_skin)

#Zhu non-porous Data ####

fit_lnorm_non_porous <- fitdist(zhu_non_porous$conc, "lnorm")

summary(fit_lnorm_non_porous)

#graph
par(mfrow = c(2, 2)) 
denscomp(fit_log10norm_non_porous) 
qqcomp(fit_log10norm_non_porous) 
cdfcomp(fit_log10norm_non_porous) 
ppcomp(fit_log10norm_non_porous)

##Using log-transformed data
fit_log10norm_non_porous <- fitdist(zhu_non_porous$y, "norm")
fit_log10unif_non_porous <- fitdist(zhu_non_porous$y, "unif")
fit_log10gam_non_porous <- fitdist(zhu_non_porous$y, "gamma")
fit_log10wei_non_porous <- fitdist(zhu_non_porous$y, "weibull") 

summary(fit_log10norm_non_porous)
summary(fit_log10gam_non_porous)
summary(fit_log10wei_non_porous)

fits_compare_non_porous <- list(fit_log10norm_non_porous, fit_log10gam_non_porous, fit_log10wei_non_porous)
fit_names_non_porous <- c("normal", "gamma", "weibull")
fit_func(fits_compare_non_porous, fit_names_non_porous)

#Zhu env combined data ####

fit_lnorm_zhu_env <- fitdist(zhu_env$conc, "lnorm")

summary(fit_lnorm_zhu_env)

#graph
par(mfrow = c(2, 2)) 
denscomp(fit_lnorm_zhu_env) 
qqcomp(fit_lnorm_zhu_env) 
cdfcomp(fit_lnorm_zhu_env) 
ppcomp(fit_lnorm_zhu_env)

#Using log-transformed data
fit_log10norm_zhu_env <- fitdist(zhu_env$y, "norm")
fit_log10gam_zhu_env <- fitdist(zhu_env$y, "gamma")
fit_log10wei_zhu_env <- fitdist(zhu_env$y, "weibull") 

summary(fit_log10norm_zhu_env)
summary(fit_log10gam_zhu_env)
summary(fit_log10wei_zhu_env)

fits_compare_zhu_env <- list(fit_log10norm_zhu_env, fit_log10gam_zhu_env, fit_log10wei_zhu_env)
fit_names_zhu_env <- c("normal", "gamma", "weibull")
fit_func(fits_compare_zhu_env, fit_names_zhu_env)

#####
#Fitting distributions one at a time

#Lognormal
fln = fitdist(zhu_axilla$conc, "lnorm", method="mle") 
summary(fln)
#graph
par(mfrow = c(2, 2)) 
denscomp(fln) 
qqcomp(fln) 
cdfcomp(fln) 
ppcomp(fln)

#Normal using log-transformed data
fn_log10 = fitdist(zhu_axilla$y, "norm", method="mle")
summary(fn_log10)
#graph
par(mfrow = c(2, 2)) 
denscomp(fn_log10) 
qqcomp(fn_log10) 
cdfcomp(fn_log10) 
ppcomp(fn_log10)

#Weibull and uniform for comparison
fw = fitdist(zhu_axilla$y, "weibull", method="mle") 
summary(fw)
funi = fitdist(zhu_axilla$y, "unif", method="mle") 
summary(funi)

#Graph them all
par(mfrow = c(2, 2)) 
plot.legend = c("Weibull", "lognormal", "uniform")
denscomp(list(fw, fln, funi), legendtext = plot.legend) 
qqcomp(list(fw, fln, funi), legendtext = plot.legend) 
cdfcomp(list(fw, fln, funi), legendtext = plot.legend) 
ppcomp(list(fw, fln, funi), legendtext = plot.legend)

par(mfrow = c(2, 2)) 
denscomp(funi) 
qqcomp(funi) 
cdfcomp(funi) 
ppcomp(funi)

#END

