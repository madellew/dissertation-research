#EITS with single infectious compartment to characterize transmission  dynamics for C. difficile 
#infection (CDI) and C. auris device-related infections (parameterized for catheter-associated 
#urinary tract infections (CAUTI) and central line-associated bloodstream infections CLABSI).

#Created: 1 Feb 2024
#Last modified: 22 Jan 2025

#State variables: Susceptible (S), Infectious (I), Recovered/removed (R), and Resistant (Rs)
#For tracking purposes N sums demography (admissions and discharge) and In sums incident cases


library(deSolve)
library(ragg)

#####
#Constants

#N  - population size
#k  - kappa; contact rate
#th - theta; proportion admitted in X state
#m  - mu; discharge rate
#p  - pi; infectivity
#r  - rho; amount picked up from environment per day
#g  - gamma; recovery rate
#a  - alpha; deposition rate (CFU/day)
#l  - lambda; pathogen decay rate

#General constants
N = 100      
scenarios <- c('C-Diff', 'CAUTI', 'CLABSI')
k.env <- 14.2 #rate at which patients contact the environment
k.cath <- 12 #self-inoculation rate, urinary catheter (CAUTI)
k.cath2 <- 2 #self-inoculation rate, CVC (CLABSI)
k.mouth <- 34.8 #self-inoculation rate, mouth

#Constants for C. difficile
th.r.diff <- 0.3 #non-susceptible
m.diff <- 0.08 
p.diff <- 5.312759E-5
r.diff <- 0.57*k.env*0.4099*k.mouth
g.diff <- 1/14.7 
a.diff <- 1900
l.diff <- 1/28

#Constants for C. auris; CAUTI
th.r.cauti <- 0.292
m.cauti <- 0.08
p.cauti <- 9.086E-10 
#r.cauti <- 0.001757*k.env*0.008131*k.cath
r.cauti <- 0.3*k.env*0.1026*k.cath
g.cauti <- 1/17.6 
a.cauti <- 76*10000 
l.cauti <- 1/14
  
#Constants for C. auris; CLABSI --
th.r.clabsi <- 0.44
m.clabsi <- 0.08
p.clabsi <- 4.529E-5 
#r.clabsi <- 0.001757*k.env*0.008131*k.cath
r.clabsi <- 0.0075*k.env*0.008131*k.cath2
g.clabsi <- 1/22.5 
a.clabsi <- 76*10000
l.clabsi <- 1/14
  
#####
#Model equations
  
  EITSModel <- function(t, y, param){
    S  <- y[1]
    I  <- y[2]
    R  <- y[3] 
    R2 <- y[4] 
    E  <- y[5]
    N  <- y[6] 
    In <- y[7]
    
    th.r <- param['th.r']
    m <- param['m']
    p <- param['p']
    r <- param['r']
    g <- param['g']
    a <- param['a']
    l <- param['l']
    
    with(as.list(param),{
      dSdt <- (m*(S+R2) + g*I)*(1-th.r) - 0.001*r*p*S*E - m*S
      dIdt <- 0.001*r*p*S*E - g*I
      dRdt <- g*I
      dR2dt <- (m*(S+R2) + (g*I))*th.r - m*R2
      dEdt <- a*I - E*((100*0.001*r)+l)
      dNdt <- m*(S+R2) + g*I
      dIndt <- 0.001*r*p*S*E
    
    dydt <- c(dSdt,dIdt,dRdt,dR2dt,dEdt,dNdt,dIndt)
    
    list(dydt)
    })
  }
  
#####
#ODE solution
  
#Set initial conditions
  
ini.cond.diff <- c(S=69, I=1, R=0, R2=30, E=0, N=0, In=0) 
ini.cond.cauti <- c(S=70, I=1, R=0, R2=29, E=0, N=0, In=0) 
ini.cond.clabsi <- c(S=55, I=1, R=0, R2=44, E=0, N=0, In=0)
  
times <- seq(0,365,1)

#Define parameters

param.diff <- c(th.r=th.r.diff, m=m.diff, p=p.diff, r=r.diff, g=g.diff, a=a.diff, l=l.diff)
param.cauti <- c(th.r=th.r.cauti, m=m.cauti, p=p.cauti, r=r.cauti, g=g.cauti, a=a.cauti, l=l.cauti)
param.clabsi <- c(th.r=th.r.clabsi, m=m.clabsi, p=p.clabsi, r=r.clabsi, g=g.clabsi, a=a.clabsi, l=l.clabsi)
  
#deSolve package for ode solution
  
EITS.diff <- ode(ini.cond.diff, times, EITSModel, param.diff)
EITS.cauti <- ode(ini.cond.cauti, times, EITSModel, param.cauti)
EITS.clabsi <- ode(ini.cond.clabsi, times, EITSModel, param.clabsi)

#R0

  R_0 <- function(N, param){
  
    th.r <- param['th.r']
    m <- param['m']
    p <- param['p']
    r <- param['r']
    g <- param['g']
    a <- param['a']
    l <- param['l']
  
    with(as.list(param),{
    
      #Completely susceptible population
      R0 <- (a/g)*((0.001*r*N)/((0.001*N*r)+l))*p
      
      #Population with proportion resistant
      R0_rs <- (a/g)*((0.001*r*N*(1-th.r))/((0.001*N*r)+l))*p
      
      R <- c(R0, R0_rs)
      R
      })
  }

R_0(100,param.diff)
R_0(100,param.cauti)
R_0(100,param.clabsi)

#####
#Outputs
  
#Summary statistics
  
diff.summary <- summary(EITS.diff)
colnames(diff.summary) <- c("Susceptible", "Infected", "Recovered", "Resistant", "Environment", "Discharges", "Incidence")
write.csv(diff.summary, file = "C-diff-Summary-Stats.csv")

cauti.summary <- summary(EITS.cauti)
colnames(cauti.summary) <- c("Susceptible", "Infected", "Recovered", "Resistant", "Environment", "Discharges", "Incidence")
write.csv(cauti.summary, file = "CAUTI-Summary-Stats.csv")

clabsi.summary <- summary(EITS.clabsi)
colnames(clabsi.summary) <- c("Susceptible", "Infected", "Recovered", "Resistant", "Environment", "Discharges", "Incidence")
write.csv(clabsi.summary, file = "CLABSI-Summary-Stats.csv")
  
#Plots

#Total infectious patients for CDI and C. auris (CLABSI parameterization)
ragg::agg_png("Inc_both.png", res=500, height=3.5, width=5, units="in", scaling=0.6)
par(mar = c(5, 4, 4, 4) + 0.2, family='serif')
leg.text <- c(expression(italic('C. difficile')),expression(italic('C. auris')))
plot(EITS.diff[,1], (EITS.diff[,3]+EITS.diff[,4]), type ="l", xlab = "Time (days)",
     ylab = "Total Infectious Patients", col = '#4B4B8FFF', lwd=3, cex.lab=1.2)
lines(EITS.clabsi[,1], (EITS.clabsi[,3]+EITS.clabsi[,4]), col = '#4B4B8FFF', lwd=3, lty=4)
grid()
legend('topleft',leg.text, col = c('#4B4B8FFF','#4B4B8FFF'), lty = c(1,4), lwd=c(2,2))
dev.off()

#Size of E
ragg::agg_png("ENV_both.png", res=500, height=3.5, width=5, units="in", scaling=0.6)
par(mar = c(5, 4, 4, 4) + 0.2, family='serif')
leg.text <- c(expression(italic('C. difficile')),expression(italic('C. auris')))
plot(EITS.diff[,1],EITS.diff[,6],type="l",xlab="Time (days)", ylab="Environment (CFU)", col="#80BA5AFF", lwd=3, cex.lab=1.2)
# set parameter new=True for a new axis
par(new = TRUE)
# Draw second plot using axis y2
plot(EITS.clabsi[,1],EITS.clabsi[,6], type="l", axes = FALSE, xlab = "", ylab = "", col="#80BA5AFF", lwd=3, lty=4, cex.lab=1.2)
axis(side = 4, at = pretty(range(EITS.clabsi[,6])))
grid()
legend('topleft',leg.text, col = c('#80BA5AFF','#80BA5AFF'), lty=c(1,4), lwd=c(2,2))
dev.off()

#Automated plotting for all scenarios
  
  for(i in 1:length(scenarios)){
    I <- scenarios[i]
    if(i==1){Z <- EITS.diff}; if(i==2){Z <- EITS.cauti}; if(i==3){Z <- EITS.clabsi}
    
    #plot I
    png(sprintf('%s-Prevalence.png',I),500,500)
    plot(Z[,1], Z[,3], type ="l", xlab = "Time (days)", ylab = "Number of Patients", col = '#3969ACFF', lwd=3)
    grid()
    dev.off()
    
    #plot E
    png(sprintf('%s-Environment.png',I),500,500)
    plot(Z[,1], Z[,6], type ="l", xlab = "Time (days)", ylab = "Environment (CFU)", col = '#80BA5AFF', lwd=3)
    grid()
    dev.off()
    
    #plot I+R
    png(sprintf("%s-Cumulative-Incidence.png",I),500,500) 
    plot(Z[,1], (Z[,3]+Z[,4]), type ="l", xlab = "Time (days)", ylab = "Total Infectious Patients", col = '#4B4B8FFF', lwd=3)
    grid()
    dev.off()
    
    #plot S, I+R, I
    leg.text <- c('Susceptible','Cumulative Incidence','Prevalence')
    png(sprintf("%s-Combined-Plot.png",I),500,500)
    plot(Z[,1], Z[,2], ylim=c(min(Z[,3]), max(100)), type ="l", xlab = "Time (days)", ylab = "Number of Patients", col = '#E73F74FF', lwd=3, lty=2)
    lines(Z[,1], (Z[,3]+Z[,4]), col = '#4B4B8FFF', lwd=3)
    lines(Z[,1], Z[,3], col = '#3969ACFF', lwd =3, lty=4)
    grid()
    legend('topleft',leg.text, col = c('#E73F74FF','#4B4B8FFF','#3969ACFF'), lty = c(2,1,4), lwd=(2,2,2))
    dev.off()
  }

#Interactive plotting
source('shiny-app-mod.R')   


#END
