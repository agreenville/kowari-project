
library("MARSS")

###################################################
### code chunk number 6: Cs1_Exercise1
# for graphing each sim run.
###################################################
# par(mfrow=c(3,3))
# sim.u = mean.u #-0.05 
# sim.Q = mean.Q #0.02 
# sim.R = mean.R #0.05 
# nYr= 50 
# fracmiss = 0.1 
# init = kow.x.ZC.1.TN.r[1] #7 
# years = seq(1:nYr)
# for(i in 1:9){
#   x = rep(NA,nYr) # vector for ts w/o measurement error 
#   y = rep(NA,nYr) # vector for ts w/ measurement error
#   x[1]=init			
#   for(t in 2:nYr){	
#     x[t] = x[t-1]+ sim.u + rnorm(1, mean=0, sd=sqrt(sim.Q)) }
#   for(t in 1:nYr){ 
#     y[t]= x[t] + rnorm(1,mean=0,sd=sqrt(sim.R)) }
#   missYears = 
#     sample(years[2:(nYr-1)],floor(fracmiss*nYr),replace = FALSE)
#   y[missYears]=NA
#   plot(years, y,
#        xlab="",ylab="log abundance",lwd=2,bty="l")
#   lines(years,x,type="l",lwd=2,lty=2)
#   title(paste("simulation ",i) )
# }
# legend("topright", c("Observed","True"),
#        lty = c(-1, 2), pch = c(1, -1))


###################################################
### code chunk number 16: Cs1_Exercise2
# simulate projected population parameters with growth rate,
# process and obs error from MARSS 1-state model in kowari_script.
# Dennis is process error only.
# MARSS has both process and obs error.
# up-dated for Bayesian methods. See HPC pop sim script to run
# 1000s of simulations.
###################################################
sim.u = mean.u #-0.05   # growth rate
sim.Q = mean.Q #0.02    # process error variance
sim.R = mean.R #0.05    # non-process error variance
nYr= 100         # number of years of data to generate
fracmiss = 0.1  # fraction of years that are missing
init = kow.x.ZC.1.TN[1] #7        # log of initial pop abundance 
nsim = 1000
years = seq(1:nYr)  # col of years
params = matrix(NA, nrow=(nsim+2), ncol=5, 
                dimnames=list(c(paste("sim",1:nsim),"mean sim","true"),
                              c("kem.U","den91.U","kem.Q","kem.R", "den91.Q")))
x.ts = matrix(NA,nrow=nsim,ncol=nYr)  # ts w/o measurement error 
y.ts = matrix(NA,nrow=nsim,ncol=nYr)  # ts w/ measurement error
for(i in 1:nsim){ 
  x.ts[i,1]=init			
  for(t in 2:nYr){	
    x.ts[i,t] = x.ts[i,t-1]+sim.u+rnorm(1,mean=0,sd=sqrt(sim.Q))}
  for(t in 1:nYr){ 
    y.ts[i,t] = x.ts[i,t]+rnorm(1,mean=0,sd=sqrt(sim.R))}
  missYears = sample(years[2:(nYr-1)], floor(fracmiss*nYr),
                     replace = FALSE) 
  y.ts[i,missYears]=NA
  
  #MARSS estimates 
  kem = bayes.marss.1.fn(y.ts[i,])   #MARSS(y.ts[i,], silent=TRUE)
  #type=vector outputs the estimates as a vector instead of a list
  params[i,c(1,3,4)] = unlist(kem$BUGSoutput$mean[c(3,5,6)]) #params[i,c(1,3,4)] = coef(kem,type="vector")[c(2,3,1)]
  
  #Dennis et al 1991 estimates
  den.years = years[!is.na(y.ts[i,])]  # the non missing years
  den.yts = y.ts[i,!is.na(y.ts[i,])]   # the non missing counts
  den.n.yts = length(den.years) 	
  delta.pop = rep(NA, den.n.yts-1 ) # transitions
  tau = rep(NA, den.n.yts-1 )       # time step lengths
  for (t in 2:den.n.yts ){
    delta.pop[t-1] = den.yts[t] - den.yts[t-1] # transitions
    tau[t-1] =  den.years[t]-den.years[t-1]    # time step length
  } # end i loop
  den91 <- lm(delta.pop ~ -1 + tau) # -1 specifies no intercept
  params[i,c(2,5)] = c(den91$coefficients, var(resid(den91))) 
}
params[nsim+1,]=apply(params[1:nsim,],2,mean)
params[nsim+2,]=c(sim.u,sim.u,sim.Q,sim.R,sim.Q)


###################################################
### code chunk number 20: Cs1_Exercise3
# Dennis diffusion process esitmates of extinction risk
###################################################
#Needs Example 2 to be run first

########################################################
## function for diffusion process esitmates of extinction risk.
extinc.risk.fn <- function(pd, te,sim.u, sim.Q, params){
  xd = -log(pd) 
  tyrs = 1:te 
  marss.sim <-matrix(nrow = te, ncol = nsim+1)
  denn.sim <-matrix(nrow = te, ncol = nsim+1)
  for(j in c(nsim+1,1:nsim)){ # c(10, 1:8)
  real.ex = denn.ex = kal.ex = matrix(nrow=te) 
  
  #MARSS parameter estimates
  u=params[j,1];   Q=params[j,3]
  if(Q==0) Q=1e-4  #just so the extinction calc doesn't choke
  p.ever = ifelse(u<=0,1,exp(-2*u*xd/Q))
  for (i in 1:te){      
    if(is.finite(exp(2*xd*abs(u)/Q))){
      sec.part = exp(2*xd*abs(u)/Q)*pnorm((-xd-abs(u)* tyrs[i])/sqrt(Q*tyrs[i]))
    }else sec.part=0      
    kal.ex[i]=p.ever*pnorm((-xd+abs(u)*tyrs[i])/sqrt(Q*tyrs[i]))+sec.part
    
    marss.sim[,j] <- kal.ex #matrix to store sim runs. 10th (i.e. nsim +1) is average
    
  } # end i loop
  
  #Dennis et al 1991 parameter estimates
  u=params[j,2];   Q=params[j,5]
  p.ever = ifelse(u<=0,1,exp(-2*u*xd/Q)) 
  for (i in 1:te){      
    denn.ex[i]=p.ever*pnorm((-xd+abs(u)*tyrs[i])/sqrt(Q*tyrs[i]))+
      exp(2*xd*abs(u)/Q)*pnorm((-xd-abs(u)*tyrs[i])/sqrt(Q*tyrs[i]))
    
    denn.sim[,j] <- denn.ex #matrix to store sim runs as colns. nsim+1 is the average of runs 
    
  } # end i loop
  
  #True parameter values
  u=sim.u;   Q=sim.Q
  p.ever = ifelse(u<=0,1,exp(-2*u*xd/Q)) 
  for (i in 1:te){      
    real.ex[i]=p.ever*pnorm((-xd+abs(u)*tyrs[i])/sqrt(Q*tyrs[i]))+
      exp(2*xd*abs(u)/Q)*pnorm((-xd-abs(u)*tyrs[i])/sqrt(Q*tyrs[i]))
  } # end i loop
}
return(list(tyrs, marss.sim, denn.sim, real.ex ))
} # end function

########################################
# calc PVA's

# 90% prob of pop decline over 100 years
extinct.100.90 <- extinc.risk.fn(pd=0.1, te = 100, sim.u=sim.u, sim.Q = sim.Q, params = params)

# IUCN criterion A: % reductions over 10 years
extinct.90 <- extinc.risk.fn(pd=0.1, te = 10, sim.u=sim.u, sim.Q = sim.Q, params = params)
# critically endangered
extinct.80 <- extinc.risk.fn(pd=0.2, te = 10, sim.u=sim.u, sim.Q = sim.Q, params = params)
# endangered
extinct.50 <- extinc.risk.fn(pd=0.5, te = 10, sim.u=sim.u, sim.Q = sim.Q, params = params)
# vulnerable
extinct.30 <- extinc.risk.fn(pd=0.7, te = 10, sim.u=sim.u, sim.Q = sim.Q, params = params)

# IUCN criterion E: 99% reduction (Extinction) prob over years
# critically endangered 50% chance of extinction in 10 years
extinct.E.CR <- extinc.risk.fn(pd=0.01, te = 10, sim.u=sim.u, sim.Q = sim.Q, params = params)
# endangered 20% chance of extinction in 20 years
extinct.E.EN <- extinc.risk.fn(pd=0.01, te = 20, sim.u=sim.u, sim.Q = sim.Q, params = params)
# vulnerable 10% chance of extinction in 100 years
extinct.E.VU <- extinc.risk.fn(pd=0.01, te = 100, sim.u=sim.u, sim.Q = sim.Q, params = params)

#############################################
# graphing
#########################################

# function for graphing

extinct.graph.fn <-function(out, lab){
  plot(out[[1]], out[[4]], xlab="", 
       ylab="", ylim=c(0,1), bty="l")
  title(paste(lab))
  lines(out[[1]],out[[2]][,nsim+1],type="l",col="red",lwd=2,lty=1)
  lines(out[[1]], apply(out[[2]][,1:nsim],1,quantile,0.975), lty=2, col="red")
  lines(out[[1]], apply(out[[2]][,1:nsim],1,quantile,0.025),lty=2, col="red")
  
  lines(out[[1]],out[[3]][,nsim+1],type="l",col="black",lwd=2,lty=1)
  lines(out[[1]], apply(out[[3]][,1:nsim],1,quantile,0.975), lty=2, col="black")
  lines(out[[1]], apply(out[[3]][,1:nsim],1,quantile,0.025),lty=2, col="black")
}


# IUCN criterion A

#png(filename = "output/fig_PVA_Bayes_CrA.png", width = 120, height = 170, units = 'mm', res = 300) 

par(mfrow=c(3,1), mar=c(5.1, 4.1, 4.1, 9.5))
extinct.graph.fn(extinct.80, "CR: 80% population decline")
  abline(h=1, lty=2, col="blue", xpd=F)
  legend("topright",c("1-state model","KalmanEM", "Dennis", "KalmanEM 95% CI", "Dennis 95% CI"),
       pch=c(1,-1,-1, -1, -1), col=c(1,"red", "black", "red", "black"),
       lty=c(-1,1,1,2,2),lwd=c(-1,2,2,1,1),bty="n",xpd=T, inset=c(-0.44,-0.1), box.col="white")
  
extinct.graph.fn(extinct.50, "EN: 50% population decline")
  abline(h=1, lty=2, col="blue", xpd=F)
extinct.graph.fn(extinct.30, "VU: 30% population decline")
  abline(h=1, lty=2, col="blue", xpd=F)
par(mfrow = c(1,1), mar=c(5, 4, 4, 2) + 0.1)
mtext("Probability of population decline", 2, line=3)
mtext("Years", 1, line=3)

#dev.off()  


# criterion E

#png(filename = "output/fig_PVA_bayes_CrE.png", width = 120, height = 170, units = 'mm', res = 300) 

par(mfrow=c(3,1), mar=c(5.1, 4.1, 4.1, 9.5))
extinct.graph.fn(extinct.E.CR, "CR: 50% chance of extinction in 10 years")
  legend("topright",c("1-state model","KalmanEM", "Dennis", "KalmanEM 95% CI", "Dennis 95% CI"),
       pch=c(1,-1,-1, -1, -1), col=c(1,"red", "black", "red", "black"),
       lty=c(-1,1,1,2,2),lwd=c(-1,2,2,1,1),bty="n", xpd=T, inset=c(-0.44,-0.1))
    abline(h=0.5, lty=2, col="blue",  xpd=F)
extinct.graph.fn(extinct.E.EN, "EN: 20% chance of extinction in 20 years")
  abline(h=0.2, lty=2, col="blue", xpd=F)
extinct.graph.fn(extinct.E.VU, "VU: 10% chance of extinction in 100 years")
  abline(h=0.1, lty=2, col="blue", xpd=F)
par(mfrow = c(1,1), mar=c(5, 4, 4, 2) + 0.1)
mtext("Probability of 99% population decline", 2, line=3)
mtext("Years", 1, line=3)

#dev.off()  






