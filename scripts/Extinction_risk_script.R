
library("MARSS")

###################################################
### code chunk number 6: Cs1_Exercise1
###################################################
par(mfrow=c(3,3))
sim.u = mean.u #-0.05 
sim.Q = mean.Q #0.02 
sim.R = mean.R #0.05 
nYr= 50 
fracmiss = 0.1 
init = kow.x.ZC.1.TN.r[1] #7 
years = seq(1:nYr)
for(i in 1:9){
  x = rep(NA,nYr) # vector for ts w/o measurement error 
  y = rep(NA,nYr) # vector for ts w/ measurement error
  x[1]=init			
  for(t in 2:nYr){	
    x[t] = x[t-1]+ sim.u + rnorm(1, mean=0, sd=sqrt(sim.Q)) }
  for(t in 1:nYr){ 
    y[t]= x[t] + rnorm(1,mean=0,sd=sqrt(sim.R)) }
  missYears = 
    sample(years[2:(nYr-1)],floor(fracmiss*nYr),replace = FALSE)
  y[missYears]=NA
  plot(years, y,
       xlab="",ylab="log abundance",lwd=2,bty="l")
  lines(years,x,type="l",lwd=2,lty=2)
  title(paste("simulation ",i) )
}
legend("topright", c("Observed","True"),
       lty = c(-1, 2), pch = c(1, -1))


###################################################
### code chunk number 16: Cs1_Exercise2
###################################################
sim.u = mean.u #-0.05   # growth rate
sim.Q = mean.Q #0.02    # process error variance
sim.R = mean.R #0.05    # non-process error variance
nYr= 100         # number of years of data to generate
fracmiss = 0.1  # fraction of years that are missing
init = kow.x.ZC.1.TN[1] #7        # log of initial pop abundance (~1100 individuals)
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
  kem = MARSS(y.ts[i,], silent=TRUE)
  #type=vector outputs the estimates as a vector instead of a list
  params[i,c(1,3,4)] = coef(kem,type="vector")[c(2,3,1)]
  
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
###################################################
#Needs Example 2 to be run first
# par(mfrow=c(3,3))
pd = 0.1; xd = -log(pd)   # decline threshold 90%
te = 100; tyrs = 1:te   # extinction time horizon 
marss.sim <-matrix(nrow = te, ncol = nsim+1)
denn.sim <-matrix(nrow = te, ncol = nsim+1)
for(j in c(nsim+1,1:nsim)){ # c(10, 1:8)
  real.ex = denn.ex = kal.ex = matrix(nrow=te) 
  
  #MARSS parameter estimates
  u=params[j,1];   Q=params[j,3]
  if(Q==0) Q=1e-4  #just so the extinction calc doesn't choke
  p.ever = ifelse(u<=0,1,exp(-2*u*xd/Q))
  for (i in 1:100){      
    if(is.finite(exp(2*xd*abs(u)/Q))){
      sec.part = exp(2*xd*abs(u)/Q)*pnorm((-xd-abs(u)* tyrs[i])/sqrt(Q*tyrs[i]))
    }else sec.part=0      
    kal.ex[i]=p.ever*pnorm((-xd+abs(u)*tyrs[i])/sqrt(Q*tyrs[i]))+sec.part
    
    marss.sim[,j] <- kal.ex #matrix to store sim runs. 10th (i.e. nsim +1) is average
    
  } # end i loop
  
  #Dennis et al 1991 parameter estimates
  u=params[j,2];   Q=params[j,5]
  p.ever = ifelse(u<=0,1,exp(-2*u*xd/Q)) 
  for (i in 1:100){      
    denn.ex[i]=p.ever*pnorm((-xd+abs(u)*tyrs[i])/sqrt(Q*tyrs[i]))+
      exp(2*xd*abs(u)/Q)*pnorm((-xd-abs(u)*tyrs[i])/sqrt(Q*tyrs[i]))
    
    denn.sim[,j] <- denn.ex #matrix to store sim runs as colns. nsim+1 is the average of runs 
    
  } # end i loop
  
  #True parameter values
  u=sim.u;   Q=sim.Q
  p.ever = ifelse(u<=0,1,exp(-2*u*xd/Q)) 
  for (i in 1:100){      
    real.ex[i]=p.ever*pnorm((-xd+abs(u)*tyrs[i])/sqrt(Q*tyrs[i]))+
      exp(2*xd*abs(u)/Q)*pnorm((-xd-abs(u)*tyrs[i])/sqrt(Q*tyrs[i]))
  } # end i loop
  
#   #plot it
#   plot(tyrs, real.ex, xlab="time steps into future", 
#        ylab="probability of extinction", ylim=c(0,1), bty="l")
#   if(j<=8) title(paste("simulation ",j) )
#   if(j==10) title("average over sims")
#   lines(tyrs,denn.ex,type="l",col="red",lwd=2,lty=1) 
#   lines(tyrs,kal.ex,type="l",col="blue",lwd=2,lty=2)
   }
# legend("bottomright",c("True","Dennis","KalmanEM"),pch=c(1,-1,-1),
#        col=c(1,2,"blue"),lty=c(-1,1,2),lwd=c(-1,2,2),bty="n")
# par(mfrow=c(1,1))

#############################################
plot(tyrs, real.ex, xlab="Years", 
     ylab="probability of extinction", ylim=c(0,1), bty="l")
 title("Average over 1000 sims")
#lines(tyrs,denn.ex,type="l",col="red",lwd=2,lty=1) 
lines(tyrs,marss.sim[,nsim+1],type="l",col="blue",lwd=2,lty=1)
lines(tyrs, apply(marss.sim[,1:nsim],1,quantile,0.975), lty=2, col="red")
lines(tyrs, apply(marss.sim[,1:nsim],1,quantile,0.025),lty=2, col="red")

lines(tyrs,denn.sim[,nsim+1],type="l",col="black",lwd=2,lty=1)
lines(tyrs, apply(denn.sim[,1:nsim],1,quantile,0.975), lty=2, col="black")
lines(tyrs, apply(denn.sim[,1:nsim],1,quantile,0.025),lty=2, col="black")

legend("bottomright",c("1-state model","KalmanEM", "Dennis", "KalmanEM 95% CI", "Dennis 95% CI"),
       pch=c(1,-1,-1, -1, -1), col=c(1,"blue", "black", "red", "black"),
       lty=c(-1,1,1,2,2),lwd=c(-1,2,2,1,1),bty="n")


