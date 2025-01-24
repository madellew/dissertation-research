#EITS with multiple infectious compartments (C,I,T) to characterize transmission  
#dynamics for C. difficile infection (CDI).

#Created: 1 Feb 2024
#Last modified: 22 Jan 2025

#State variables: Susceptible (S), Colonized (C), Symptomatic (I), Treated (T), 
#Recovered/removed (R), and Resistant (Rs)
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
#o  - omega; disease progression rate
#n  - eta; treatment rate
#g  - gamma; recovery rate
#a  - alpha; deposition rate (CFU/day)
#l  - lambda; pathogen decay rate

#General constants
N = 100         # Fixed population
scenarios <- 'C-Diff'
#scenarios <- c('C-Diff', 'C-Auris')
k.env <- 14.2 #rate at which patients contact the environment
#k.cath <- 2 #self-inoculation rate, catheter
k.mouth <- 34.8 #self-inoculation rate, mouth

#Constants for C. difficile
th.r.diff <- 0.3 #non-susceptible
th.c.diff <- 0 #colonized
m.diff <- 0.08 
p.diff <- 5.312759E-5
r.diff <- 0.57*k.env*0.4099*k.mouth
o.diff <- 0.25
n.diff <- 1/1.75
g.diff <- 1/8.25
#a1.diff <- 0.00196*10000 #asymptomatic
a2.diff <- 0.19*10000 #symptomatic
#a3.diff <- 1.479E-4*10000 #treated
l.diff <- 1/28

#test increasing shedding rate
a1.diff <- 0.7*10000
a3.diff <- 0.095*10000

#Constants for C. auris
#th.r.auris <- 0.3
#th.c.auris <- 0
#m.auris <- 0.08
#p.auris <- 4.529E-5 
#r.auris <- 0.001757*k.env*0.008131*k.cath
#r.auris <- 0.3*k.env*0.008131*k.cath
#o.auris <- 1/7
#n.auris <- 1/2
#g.auris <- 1/10 
#a1.auris <- 76*10000 
#a2.auris <- 220.8*10000
#a3.auris <- 1.919*10000
#l.auris <- 1/14

#####
#Model equations
    
  EITSModel.Inf <- function(t, y, param){
    S  <- y[1]
    C  <- y[2]
    I  <- y[3]
    Tr <- y[4]
    R  <- y[5]
    R2 <- y[6]
    E  <- y[7]
    N  <- y[8]
    In <- y[9]
    
    
    th.r <- param['th.r']
    th.c <- param['th.c']
    m <- param['m']
    p <- param['p']
    r <- param['r']
    o <- param['o']
    n <- param['n']
    g <- param['g']
    a1 <- param['a1']
    a2 <- param['a2']
    a3 <- param['a3']
    l <- param['l']
    
    with(as.list(param),{  
      dSdt <- (m*(S+C+R2) + g*Tr)*(1-th.r-th.c) - 0.001*r*p*S*E - m*S
      dCdt <- (m*(S+C+R2) + g*Tr)*th.c + 0.001*r*p*S*E - (m+o)*C
      dIdt <- o*C - n*I
      dTrdt <- n*I - g*Tr
      dRdt <- g*Tr
      dR2dt <- (m*(S+C+R2) + g*Tr)*th.r - m*R2
      dEdt <- (a1*C + a2*I + a3*Tr) - E*((100*(0.001*r))+l)
      dNdt <- m*(S+C+R2) + g*Tr
      dIndt <- 0.001*r*p*S*E
      
      dydt <- c(dSdt,dCdt,dIdt,dTrdt,dRdt,dR2dt,dEdt,dNdt,dIndt)

      
      list(dydt)
    })
  }
    
#####
#ODE solution

#Set initial conditions
    
ini.cond.diff <- c(S=69, C=0, I=1, Tr=0, R=0, R2=30, E=0, N=0, In=0)
#ini.cond.auris <- c(S=69, C=0, I=1, Tr=0, R=0, R2=30, E=0, N=0, In=0)

times <- seq(0,365,1)

#Define parameters
param.diff <- c(th.r=th.r.diff, th.c=th.c.diff, m=m.diff, p=p.diff, r=r.diff, o=o.diff, 
                n=n.diff, g=g.diff, a1=a1.diff, a2=a2.diff, a3=a3.diff, l=l.diff)
#param.auris <- c(th.r=th.r.auris, th.c=th.c.auris, m=m.auris, p=p.auris, r=r.auris, o=o.auris, 
#                n=n.auris, g=g.auris, a1=a1.auris, a2=a2.auris, a3=a3.auris, l=l.auris)   
    
#deSolve package for ode solution
    
EITS.Inf.diff <- ode(ini.cond.diff,times,EITSModel.Inf,param.diff)
#EITS.Inf.auris <- ode(ini.cond.auris,times,EITSModel.Inf,param.auris)


#####
#Outputs

#Summary statistics
    
diff.summary <- summary(EITS.Inf.diff)
colnames(diff.summary) <- c("Susceptible", "Col", "Symp", "Treated", "Recovered", "Resistant", "Environment", "Discharges", "Inc")
write.csv(diff.summary, file = "C-diff-Summary-Stats-inf.csv")

#auris.summary <- summary(EITS.Inf.auris)
#colnames(auris.summary) <- c("Susceptible", "Col", "Symp", "Treated", "Recovered", "Resistant", Environment", "Discharges", "Inc")
#write.csv(auris.summary, file = "C-auris-Summary-Stats-inf.csv")
    
#Plots
    
#Total infectious patients
ragg::agg_png("I_CDI-INF.png", res=500, height=3.5, width=5, units="in", scaling=0.6)
par(family='serif')
leg.text <- c("Infectious", "Colonized", "Symptomatic", "Treated")
plot(EITS.Inf.diff[,1], (EITS.Inf.diff[,3]+EITS.Inf.diff[,4]+EITS.Inf.diff[,5]), ylim=c(0,6),
     type ="l", xlab = "Time (days)", ylab = "Total Patients", col = '#4B4B8FFF', lwd=3, cex.lab=1.2)
lines(EITS.Inf.diff[,1], EITS.Inf.diff[,3], type ="l",col = '#E68310FF', lwd=3, lty=2)
lines(EITS.Inf.diff[,1], EITS.Inf.diff[,4], type ="l",col = '#7F3C8DFF', lwd=3, lty=3)
lines(EITS.Inf.diff[,1], EITS.Inf.diff[,5], type ="l",col = '#11A579FF', lwd=3, lty=4)
grid()
legend('topleft',leg.text, col = c('#4B4B8FFF','#E68310FF', '#7F3C8DFF', '#11A579FF'), lty = c(1,2,3,4), lwd=c(2,2,2,2))
dev.off()

#Size of E
ragg::agg_png("E_CDI-INF.png", res=500, height=3.5, width=5, units="in", scaling=0.6)
par(family='serif')
plot(EITS.Inf.diff[,1],EITS.Inf.diff[,8], type="l", xlab="Time (days)", ylab="Environment (CFU)", col="#80BA5AFF", lwd=3, cex.lab=1.2)
grid()
dev.off()
    

#Automated plotting for all scenarios    
    
  for(i in 1:length(scenarios)){
    I <- scenarios[i]
    if(i==1){Z <- EITS.Inf.diff}; if(i==2){Z <- EITS.Inf.auris}
    
    #plot infectious (C, I, T)
    leg.text <- c('Colonized','Infected','Treated')
    png(sprintf("%s-Inf-Plot-Inf.png",I),500,500)
    plot(Z[,1],Z[,3],ylim=c(0, max(Z[,5])),type="l",xlab="Time (days)", ylab="Number of patients", col='#E68310FF', lwd=3, lty=2) 
    lines(Z[,1], Z[,4], col = '#7F3C8DFF', lwd = 3, lty=3)
    lines(Z[,1], Z[,5], col = '#11A579FF', lwd =3, lty=3)
    grid()
    legend('topleft',leg.text, col = c('#E68310FF', '#7F3C8DFF', '#11A579FF'), lty = c(2,3,4), lwd=c(2,2,2))
    dev.off()
    
    #plot E
    png(sprintf('%s-Environment-Inf.png',I),500,500)
    plot(Z[,1], Z[,8], type ="l", xlab = "Time (days)", ylab = "Environment (CFU)", col = '#80BA5AFF', lwd=3)
    grid()
    dev.off()
    
    #plot C+I+T+R
    png(sprintf("%s-Cumulative-Incidence-Inf.png",I),500,500) 
    plot(Z[,1], (Z[,3]+Z[,4]+Z[,5]+Z[,6]), type ="l", xlab = "Time (days)", ylab = "Total Infectious Patients", col = '#4B4B8FFF', lwd=3)
    grid()
    dev.off()
    
    #plot S, C++T+R, C+I+T
    leg.text <- c('Susceptible','Cumulative Incidence','Prevalence')
    png(sprintf("%s-Combined-Plot-Inf.png",I),500,500)
    plot(Z[,1], Z[,2], ylim=c(0, max((Z[,3]+Z[,4]+Z[,5]+Z[,6]))), type ="l", xlab = "Time (days)", ylab = "Number of patients", col = '#E73F74FF', lwd=3, lty=2)
    lines(Z[,1], (Z[,3]+Z[,4]+Z[,5]+Z[,6]), col = '#4B4B8FFF', lwd=3)
    lines(Z[,1], (Z[,3]+Z[,4]+Z[,5]), col = '#3969ACFF', lwd =3, lty=4)
    grid()
    legend('topleft',leg.text, col = c('#E73F74FF','#4B4B8FFF','#3969ACFF'), lty = c(2,1,4), lwd=c(2,2,2))
    dev.off()
    }


#Interactive plotting
source('shiny-app-mod.R')

    
#END
    