#Using shinySIR for Interactive plotting of EITS models
#Ref: Sinead E. Morris (2020). shinySIR: Interactive Plotting for Mathematical Models of Infectious Disease Spread. R package version 0.1.2.

#Created: 1 Sep 2024
#Last modified: 27 Nov 2024


#Load SIR Shiny App
library(shinySIR)
#Modified functions
source('shiny-app-mod.R')

#####

#run_shiny(model = #model title 
#          neweqns = #specify function with model equations 
#          tstart= , timestep= , tmax=  #specify time (default tstart=0, timestep=1, tmax=365)
#          ics = #specify initial conditions
#          parm0 = #specify starting parameter values
#          parm_names = #parameter variable names
#          parm_min = #specify minimum parameter values
#          parm_max = #specify maximum parameter values
#          sigfigs = #specify sig figs (optional; default=4)
#          remove_cols = #specify compartments to not be plotted (optional)

#####
#EITS model with single infectious state - pick-up rate simplified

#Define model equations
EITS <- function(t, y, parms) {
  
  with(as.list(c(y, parms)),{
    
    dSdt <- (m*(S+R2) + g*I)*(1-th.r) - 0.001*r*p*S*E - m*S
    dIdt <- 0.001*r*p*S*E - g*I
    dRdt <- g*I
    dR2dt <- (m*(S+R2) + (g*I))*th.r - m*R2
    dEdt <- a*I - E*((100*0.001*r)+l)
  
    return(list(c(dSdt, dIdt, dRdt, dR2dt, dEdt)))
  })
}

#Interactive plotting 
#Parameter specification for C. difficile infections
run_shiny(model = "EITS Model", 
          neweqns = EITS, 
          tstart=0, timestep=1, tmax=365,
          ics = c(S=69, I=1, R=0, R2=30, E=0),
          parm0 = c(th.r=0.3, m=0.08, p=5.312759E-5, r=(14.2*34.8*0.4099*0.57), g=1/14.7, a=1900, l=1/28),
          parm_names = c("Proportion Resistant", "Discharge rate", "Infectivity","Pick-up rate",
                          "Recovery rate", "Deposition rate", "Decay rate"),
          parm_min = c(th.r=0, m=0, p=5.312759E-7, r=0.115457, g=1/20, a=1, l=1/360),
          parm_max = c(th.r=1, m=0.5, p=0.003, r=11545.7, g=1/2, a=19000, l=1),
          sigfigs = 10) #,
          #remove_cols = c('E'))


#####
#EITS model with single infectious state - pick-up rate as variables (can specify contact rates and transfer efficiencies)

EITS2 <- function(t, y, parms) {
  
  with(as.list(c(y, parms)),{
    
    dSdt <- (m*(S+R2) + g*I)*(1-th.r) - 0.001*(n.m*k.m*n.e*k.e)*p*S*E - m*S
    dIdt <- 0.001*(n.m*k.m*n.e*k.e)*p*S*E - g*I
    dRdt <- g*I
    dR2dt <- (m*(S+R2) + (g*I))*th.r - m*R2
    dEdt <- a*I - E*((100*0.001*n.m*k.m*n.e*k.e)+l)
    
    return(list(c(dSdt, dIdt, dRdt, dR2dt, dEdt)))
  })
}

run_shiny(model = "EITS Model", 
          neweqns = EITS2,
          tstart=0, timestep=1, tmax=365,
          ics = c(S=69, I=1, R=0, R2=30, E=0),
          parm0 = c(th.r=0.3, m=0.08, p=5.312759E-5, k.e=14.2, k.m=34.8, n.m=0.4099, n.e=0.57, g=1/14.7, a=1900, l=1/28),
          parm_names = c("Proportion Resistant", "Discharge rate", "Infectivity","Fomite contact rate", "Self-innoculation rate",
                         "Fomite transfer efficiency", "Self-innoculation transfer efficiency", "Recovery rate", 
                         "Deposition rate", "Decay rate"),
          parm_min = c(th.r=0, m=0, p=5.312759E-7, k.e=1, k.m=1, n.m=0, n.e=0, g=1/20, a=1, l=1/360),
          parm_max = c(th.r=1, m=0.5, p=0.003, k.e=100, k.m=100, n.m=1, n.e=1, g=1/2, a=19000, l=1),
          sigfigs = 10) #,
          #remove_cols = c('E')))


#####
#EITS model with multiple infectious states

EITS.Inf <- function(t, y, parms) {
  
  with(as.list(c(y, parms)),{
    
    dSdt <- (m*(S+C+R2) + g*Tr)*(1-th.r-th.c) - 0.001*r*p*S*E - m*S
    dCdt <- (m*(S+C+R2) + g*Tr)*th.c + 0.001*r*p*S*E - (m+o)*C
    dIdt <- o*C - n*I
    dTrdt <- n*I - g*Tr
    dRdt <- g*Tr
    dR2dt <- (m*(S+C+R2) + g*Tr)*th.r - m*R2
    dEdt <- (a1*C + a2*I + a3*Tr) - E*((100*0.001*r)+l)
    
    return(list(c(dSdt,dCdt,dIdt,dTrdt,dRdt,dR2dt,dEdt)))
  })
}

run_shiny(model = "EITS Infectious Compartments Model", 
          neweqns = EITS.Inf, 
          tstart=0, timestep=1, tmax=365,
          ics = c(S=69, C=0, I=1, Tr=0, R=0, R2=30, E=0),
          parm0 = c(th.r=0.3, th.c=0.1, m=0.08, p=5.312759E-5, r=(14.2*34.8*0.4099*0.57), o=1/4, n=1/1.75, 
                    g=1/8.25, a1=2500, a2=1900, a3=950, l=1/28),
          parm_names = c("Proportion Resistant", "Proportion Colonized", "Discharge rate", "Infectivity","Pick-up rate",
                         "Clinical Disease Rate", "Treatment Rate", "Recovery rate", "Deposition rate - colonized", 
                         "Deposition rate - infected", "Deposition rate - treated", "Decay rate"),
          parm_min = c(th.r=0, th.c=0, m=0, p=5.312759E-7, r=0.115457, o=1/7, n=1/10, g=1/10, a1=1, a2=1, a3=1, l=1/360),
          parm_max = c(th.r=1, th.c=0.25, m=0.5, p=0.003, r=11545.7, o=1, n=1, g=1/2, a1=180000, a2=19000, a3=19000, l=1),
          sigfigs = 10,
          remove_cols = c('E','R'))


#END

