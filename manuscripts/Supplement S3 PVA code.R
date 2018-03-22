################################################################################
################################################################################
# Example MARSS PVA simulation script for High Performance Computer.
#
# Greenville et al. (2018). Dynamics, habitat use and extinction risk of a 
#   carnivorous desert marsupial.
################################################################################
################################################################################

# Acknowledgements: Holmes, E.E., E.J. Ward & M.D. Scheuerell (2014). Analysis of multivariate
# time-series using the MARSS package. Seattle, USA.: NOAA Fisheries.

source("manuscripts/Supplement S2 JAGS code for MARSS models.R") #bayesian models

# load capture data (logged)

################################################################################
# Step 1: Run MARSS population to calculate growth rate, process and 
# observation error
# 1-state MARSS model
################################################################################
whichPop.1 <- c(1,1)
n.states.1 <- max(whichPop.1)

#We can use the vector 'whichPop' to construct our Z matrix (following MARSS notation):
# Use whichPop to construct your Z matrix
Z.1 <- matrix(0,2,n.states.1)
for(i in 1:length(whichPop.1)) {
  Z.1[i,whichPop.1[i]] = 1
}


# Y captures/ 100 TN
Y <-kowari.log

# setting up vector for number of sites and years
n.pop <- dim(Y)[1]
n.yrs <- dim(Y)[2]

model <- MARSS.1.B # model function from Supplement S2

#Then we can run the model using the code
# Specify the parameters / data list / model file name
jags.params <- c("x", "sigmaQ","sigmaR", "B", "U", "A") 
jags.data <- list("Y","n.pop","n.yrs","n.states.1","Z.1") #     


kowariMARSS <- jags.parallel(jags.data, inits = NULL, parameters.to.save= jags.params,
                    model.file=model, n.chains = 3, n.thin = 25,
                    n.burnin=60000, n.iter =100000, DIC = TRUE)

kowariMARSS

attach.jags(kowariMARSS)

#kowari population estimates
kow.x.ZC.1.TN <- apply(x[,1,],2,mean)     #means for each yr

detach.jags()

# growth rate and errors from Bayes MARSS model
mean.u <- kowariMARSS$BUGSoutput$summary[4] 
mean.Q <- kowariMARSS$BUGSoutput$summary[6]
mean.R <- kowariMARSS$BUGSoutput$summary[7]

################################################################################
# Step 2: Define JAGS Bayesian 1 state model with
#  one sub-population for simulation
################################################################################

MARSS.bayes.1 <- function(){
  
  # PRIORS
  tauQ[1] ~ dgamma(0.01,0.01);
  sigmaQ[1] <- 1/sqrt(tauQ[1]); # convert from precision to sd
  tauR[1] ~ dgamma(0.01,0.01); 
  sigmaR[1] <- 1/sqrt(tauR[1]); # convert from precision to sd
  
  for(i in 2:n.states.1) {
    # uncomment out these next two lines for independent Qs
    tauQ[i] ~ dgamma(0.01,0.01);
    sigmaQ[i] <- 1/sqrt(tauQ[i]);
    
  }
  
  for(i in 2:n.pop) {
    # diagonal and equal
    tauR[i] <- tauR[1];               
    sigmaR[i] <- 1/sqrt(tauR[i]);
  }
  
  # PROCESS MODEL
  for(i in 1:n.states.1) {
    x[i,1] ~ dnorm(0,0.1);                    # independent priors on initial states
    U[i] ~ dnorm(0,1);
    B[i] <- 1;                                # setting B=1 mean no density dependence
   
    for(t in 2:n.yrs) {
      predx[i,t] <- B[i]*x[i,t-1] + U[i];       #X = BX + U
      x[i,t] ~ dnorm(predx[i,t],tauQ[i]);
    }
  }
  
  # DATA / LIKELIHOOD MODEL (Obs model)
  A[1] <- 0;
  
  for(i in 1:n.pop) {
    for(t in 1:n.yrs) {
      # inprod does matrix multiplication
      predy[t] <- inprod(Z.1[i,],x[,t]) + A[i];               #Y=ZX+A        
      lambda[t] ~ dnorm(predy[t], tauR[i]);                #obs error 
      
      Y[t] ~ dnorm(lambda[t], tauR[i]);
    }
  }
}



################################################################################
# Step 3: Define Bayesian R function for simulations
################################################################################
bayes.marss.1.fn <- function(y){
  Y <- y
  n.pop <- 1 
  n.yrs <- nYr 
  whichPop <- 1 
  n.states.1 <- 1 
  
  Z.1 <- matrix(1) 
  
  model <- MARSS.bayes.1 
  jags.params <- c("x", "sigmaQ","sigmaR", "B", "U", "A") 
  jags.data <- list("Y","n.pop","n.yrs","n.states.1","Z.1")
  
  # Set MCMC parameters
  mcmcchains <- 3
  mcmcthin <- 25
  mcmcburn <- 60000 
  samples2Save <- 40000
  
  jags.parallel(jags.data, inits = NULL, parameters.to.save= jags.params,
                model.file=model, n.chains = mcmcchains, n.thin = mcmcthin,
                n.burnin=mcmcburn, n.iter =(mcmcburn+samples2Save), DIC = TRUE)
}

################################################################################
# Step 4: Simulate projected population parameters with growth rate,
# process and obs error from MARSS 1-state model.
# 
# MARSS has both process and obs error.
# Note: 10 000 simulations runs takes a lot of time and was run on
# a High Performance Computer. 
################################################################################
sim.u = mean.u                  # growth rate
sim.Q = mean.Q                  # process error variance
sim.R = mean.R                  # non-process error variance
nYr= 10000                      # number of years of data to generate
fracmiss = 0.1                  # fraction of years that are missing
init = kow.x.ZC.1.TN[1]         # log of initial pop abundance 
nsim = 10000                    # Set number of simulations
years = seq(1:nYr)              # col of years
params = matrix(NA, nrow=(nsim+2), ncol=3, 
                dimnames=list(c(paste("sim",1:nsim),"mean sim","true"),
                              c("kem.U","kem.Q","kem.R")))

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
  kem = bayes.marss.1.fn(y.ts[i,])   # Bayesian model function
  params[i,c(1,2,3)] = unlist(kem$BUGSoutput$mean[c(3,5,6)]) 
}
params[nsim+1,]=apply(params[1:nsim,],2,mean) # mean of simulation runs
params[nsim+2,]=c(sim.u,sim.Q,sim.R) # true values


################################################################################
# Step 5: function for diffusion process estimates of extinction risk.
################################################################################

# where:
# pd = population reduction. i.e. 0.01 = 99% reduction
# te = years
# sim.u = simulated growth rate
# sim.Q  = simulated process error
# params = simulated parameters from above

extinc.risk.fn <- function(pd, te,sim.u, sim.Q, params){
  xd = -log(pd) 
  tyrs = 1:te 
  marss.sim <-matrix(nrow = te, ncol = nsim+1)
  #denn.sim <-matrix(nrow = te, ncol = nsim+1)
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
      
      marss.sim[,j] <- kal.ex #matrix to store sim runs. i.e. nsim +1 is average
      
    } # end i loop
    
   
    #True parameter values
    u=sim.u;   Q=sim.Q
    p.ever = ifelse(u<=0,1,exp(-2*u*xd/Q)) 
    for (i in 1:te){      
      real.ex[i]=p.ever*pnorm((-xd+abs(u)*tyrs[i])/sqrt(Q*tyrs[i]))+
        exp(2*xd*abs(u)/Q)*pnorm((-xd-abs(u)*tyrs[i])/sqrt(Q*tyrs[i]))
    } # end i loop
  }
  return(list(tyrs, marss.sim, real.ex ))
} 



# IUCN criterion E: 99% reduction (Extinction) prob over years
# critically endangered 50% chance of extinction in 10 years
extinct.E.CR <- extinc.risk.fn(pd=0.01, te = 10, sim.u=sim.u, sim.Q = sim.Q, params = params)
# endangered 20% chance of extinction in 20 years
extinct.E.EN <- extinc.risk.fn(pd=0.01, te = 20, sim.u=sim.u, sim.Q = sim.Q, params = params)
# vulnerable 10% chance of extinction in 100 years
extinct.E.VU <- extinc.risk.fn(pd=0.01, te = 100, sim.u=sim.u, sim.Q = sim.Q, params = params)

################################################################################
# Step 5: function for graphing estimates of extinction risk.
################################################################################


extinct.graph.paper.fn <-function(out, lab){
  plot(out[[1]], out[[2]][,nsim+1], xlab="", 
       ylab="", ylim=c(0,1), typ="l", bty="l")
  title(paste(lab))
  lines(out[[1]], apply(out[[2]][,1:nsim],1,quantile,0.975), lty=2, col="black")
  lines(out[[1]], apply(out[[2]][,1:nsim],1,quantile,0.025),lty=2, col="black")
  }

# Figure 4
par(mfrow=c(3,1), mar=c(5.1, 4.1, 4.1, 2))
extinct.graph.paper.fn(extinct.E.CR, "")
  mtext("a)", 3, line=1, adj=-0.01)
  abline(h=0.5, lty=2, col="red",  xpd=F)
extinct.graph.paper.fn(extinct.E.EN, "")
  mtext("b)", 3, line=1, adj=-0.01)
  abline(h=0.2, lty=2, col="red", xpd=F)
extinct.graph.paper.fn(extinct.E.VU, "")
  mtext("c)", 3, line=1, adj=-0.01)
  abline(h=0.1, lty=2, col="red", xpd=F)
par(mfrow = c(1,1), mar=c(5, 4, 4, 2) + 0.1)
  mtext("Probability of 99% population decline", 2, line=3)
  mtext("Years", 1, line=3)
