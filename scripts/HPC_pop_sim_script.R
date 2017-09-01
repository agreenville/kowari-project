#######################################
# MARSS PVA pop sim script for HPC
#
# Aaron Greenville
######################################

source("R/MARSS_functions.R") #bayesian models

load("data/HPC_kowari_log.Rdata") #load capture data logged and captures/100TN

#############################################################################
# kowari 100TN data
# 1-state MARSS model
########################################################################
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

model.location.1ZC.TN <- MARSS.1.B #model function


#Then we can run the model using the code
# Specify the parameters / data list / model file name
jags.params <- c("x", "sigmaQ","sigmaR", "B", "U", "A") 
jags.data <- list("Y","n.pop","n.yrs","n.states.1","Z.1") #     


# Set MCMC parameters
mcmcchains <- 3
mcmcthin <- 25
mcmcburn <- 60000 
samples2Save <- 40000
kowariMARSS <- jags.parallel(jags.data, inits = NULL, parameters.to.save= jags.params,
                    model.file=model.location.1ZC.TN, n.chains = mcmcchains, n.thin = mcmcthin,
                    n.burnin=mcmcburn, n.iter =(mcmcburn+samples2Save), DIC = TRUE)


kowariMARSS

attach.jags(kowariMARSS)

#kowari 1 state
kow.x.ZC.1.TN <- apply(x[,1,],2,mean) #means for each yr for 1 state mulgara
kow.x.ZC.1.TN.LC <- apply(x[,1,],2,quantile,0.025)
kow.x.ZC.1.TN.UC <- apply(x[,1,],2,quantile,0.975)

#Back-transform z-scoring of captures (note values are in log space still)
# for converting x out of log space (logged captures at start)
kow.x.ZC.1.TN.r <-  exp(kow.x.ZC.1.TN)-1 
kow.x.ZC.1.TN.LC.r <- exp(kow.x.ZC.1.TN.LC)-1
kow.x.ZC.1.TN.UC.r <- exp(kow.x.ZC.1.TN.UC)-1
detach.jags()

# growth rate and errors from Bayes MARSS model
mean.u <- kowariMARSS$BUGSoutput$summary[4] 
mean.Q <- kowariMARSS$BUGSoutput$summary[6]
mean.R <- kowariMARSS$BUGSoutput$summary[7]

##############################################################################
## simulate projected population parameters with growth rate,
# process and obs error from MARSS 1-state model in kowari_script.
# Dennis is process error only.
# MARSS has both process and obs error.
###################################################
sim.u = mean.u    # growth rate
sim.Q = mean.Q     # process error variance
sim.R = mean.R     # non-process error variance
nYr= 100         # number of years of data to generate
fracmiss = 0.1  # fraction of years that are missing
init = kow.x.ZC.1.TN[1] #7        # log of initial pop abundance 
nsim = 10
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

###################################################3
save.image("workspaces/HPC_pop_sims.RData")


