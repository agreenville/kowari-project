
########################################################################################
#Functions for MARSS Bayesian models for kowari project
# 
# Y ~ Normal
# Aaron Dec 2016
######################################################################################


#Load packages
library(MASS)
library(R2jags)
library(rjags)
library(runjags)
#library(mcmcplots)



#################################################################################
# 1 state model
# two sub-pops
# 
# 
##############################################################################

MARSS.1.B <- function(){
  
  # PRIORS
  tauQ[1] ~ dgamma(0.01,0.01);
  sigmaQ[1] <- 1/sqrt(tauQ[1]); # convert from precision to sd
  tauR[1] ~ dgamma(0.01,0.01); 
  sigmaR[1] <- 1/sqrt(tauR[1]); # convert from precision to sd
  
  for(i in 2:n.states.1) {
    # uncomment out these next two lines for independent Qs
    tauQ[i] ~ dgamma(0.01,0.01);
    sigmaQ[i] <- 1/sqrt(tauQ[i]);
    #tauQ[i] <- tauQ[1];             #diagonal and equal
    #sigmaQ[i] <- 1/sqrt(tauQ[i]);
  }
  
  for(i in 2:n.pop) {
    # uncomment out these next two lines for independent Rs
    #tauR[i] ~ dgamma(0.01,0.01);
    #sigmaR[i] <- 1/sqrt(tauR[i]);
    tauR[i] <- tauR[1];               #diagonal and equal
    sigmaR[i] <- 1/sqrt(tauR[i]);
  }
  
  # PROCESS MODEL
  for(i in 1:n.states.1) {
    x[i,1] ~ dnorm(0,0.1);                    # independent priors on initial states
    U[i] ~ dnorm(0,1);
    B[i] <- 1;                                # setting B=1
    #comment out for density denpendence
    #B[i] ~ dunif(0,1)               #B is b/t 0 and 1. Density depend
    
    
    for(t in 2:n.yrs) {
      predx[i,t] <- B[i]*x[i,t-1] + U[i];       #X = BX + U
      x[i,t] ~ dnorm(predx[i,t],tauQ[i]);
    }
  }
  
  # DATA / LIKELIHOOD MODEL (Obs model)
  A[1] <- 0;
  A[2] ~ dnorm(0,1);
  
  for(i in 1:n.pop) {
    for(t in 1:n.yrs) {
      # inprod does matrix multiplication
      predy[i,t] <- inprod(Z.1[i,],x[,t]) + A[i];               #Y=ZX+A        
      #pred_lambda[i,t] <- predy[i,t] + logTrapEffort.s[i,t];    #adding offset 
      lambda[i,t] ~ dnorm(predy[i,t], tauR[i]);                #obs error only
      #log.lambda[i,t] ~ dnorm(pred_lambda[i,t], tauR[i]);       #obs error and offset in log space as x is in log-space
      
      #comment out next 4 rows for Y~neg bin
      #lambda.nb[i,t] <- exp(log.lambda[i,t])                      #lambda out of log-space
      #p[i,t] <- r[i,t]/(r[i,t]+lambda.nb[i,t])                # calc p, r=size
      #Y[i,t] ~ dnegbin(p[i,t],r[i,t])                         # Y from bin neg dist
      #r[i,t] ~ dunif(0,100)                                  # setting r parameter - work on this
      
      #comment out row below for Y~Poisson
      #Y[i,t] ~ dpois(exp(log.lambda[i,t]));                      #obs from a poisson dist. exp is putting lambda on Y scale as x in log space
      
      #comment out row below for Y~normal and comment in log.lambda line above
      Y[i,t] ~ dnorm(lambda[i,t], tauR[i]);
      
      #offset trap effort distr. ~normal
      #TrapEffort.mu[i,t] <- 6.733511                          # mean of mean(logTrapEffort.s, na.rm=T)
      #TrapEffort.sigma[i,t] <-1/(1.028217)^2;               # sd from sd(logTrapEffort.s, na.rm=T) 
      #logTrapEffort.s[i,t] ~  dgamma(64.359377, 9.019594); #dnorm(TrapEffort.mu[i,t], TrapEffort.sigma[i,t]); #logTrap effort follows a normal dist (or trap effort follows a log-normal dist - from hist)
      
      #For Posterior predictive loss
      #Generate replicated data
      #Y.rep[i,t] ~ dpois(exp(log.lambda[i,t])); 
      
      #Compute squared difference b/t observed and replicated data
      #sqdiff[i,t] <- pow(Y[i,t]-Y.rep[i,t], 2);
    }
    
  }
  #For posterior predictive loss, sum squared differences across all years
  #Dsum <- sum(sqdiff);
  
 }

####################\
# 2 state model
#
# two sub-pops
##############################################################################

MARSS.2.B <- function(){
  
  # PRIORS
  tauQ[1] ~ dgamma(0.01,0.01);
  sigmaQ[1] <- 1/sqrt(tauQ[1]); # convert from precision to sd
  tauR[1] ~ dgamma(0.01,0.01); 
  sigmaR[1] <- 1/sqrt(tauR[1]); # convert from precision to sd
  
  for(i in 2:n.states.1) {
    # uncomment out these next two lines for independent Qs
    tauQ[i] ~ dgamma(0.01,0.01);
    sigmaQ[i] <- 1/sqrt(tauQ[i]);
    #tauQ[i] <- tauQ[1];             #diagonal and equal
    #sigmaQ[i] <- 1/sqrt(tauQ[i]);
  }
  
  for(i in 2:n.pop) {
    # uncomment out these next two lines for independent Rs
    #tauR[i] ~ dgamma(0.01,0.01);
    #sigmaR[i] <- 1/sqrt(tauR[i]);
    tauR[i] <- tauR[1];               #diagonal and equal
    sigmaR[i] <- 1/sqrt(tauR[i]);
  }
  
  # PROCESS MODEL
  for(i in 1:n.states.1) {
    x[i,1] ~ dnorm(0,0.1);          # independent priors on initial states
    U[i] ~ dnorm(0,1);
    B[i] <- 1;                      # B =1
    #comment out for density denpendence
    #B[i] ~ dunif(0,1)               #B is b/t 0 and 1. Density depend
    
    for(t in 2:n.yrs) {
      predx[i,t] <- B[i]*x[i,t-1] + U[i]; #X = BX + U
      x[i,t] ~ dnorm(predx[i,t],tauQ[i]);
    }
  }
  
  # DATA / LIKELIHOOD MODEL (Obs model)
  A[1] <- 0;
  A[2] <- 0;
  
  
  
  for(i in 1:n.pop) {
    for(t in 1:n.yrs) {
      # inprod does matrix multiplication
      predy[i,t] <- inprod(Z.1[i,],x[,t]) + A[i];               #Y=ZX+A        
      #pred_lambda[i,t] <- predy[i,t] + logTrapEffort.s[i,t];    #adding offset 
      lambda[i,t] ~ dnorm(predy[i,t], tauR[i]);                #obs error only
      #log.lambda[i,t] ~ dnorm(pred_lambda[i,t], tauR[i]);       #obs error and offset in log space as x is in log-space
      
      #comment out next 4 rows for Y~neg bin
      #lambda.nb[i,t] <- exp(log.lambda[i,t])                      #lambda out of log-space
      #p[i,t] <- r[i,t]/(r[i,t]+lambda.nb[i,t])                # calc p, r=size
      #Y[i,t] ~ dnegbin(p[i,t],r[i,t])                         # Y from bin neg dist
      #r[i,t] ~ dunif(0,100)                                  # setting r parameter - work on this
      
      #comment out row below for Y~Poisson
      #Y[i,t] ~ dpois(exp(log.lambda[i,t]));                      #obs from a poisson dist. exp is putting lambda on Y scale as x in log space
      
      #comment out row below for Y~normal and comment in log.lambda line above
      Y[i,t] ~ dnorm(lambda[i,t], tauR[i]);
      
      #offset trap effort distr. ~normal
      #TrapEffort.mu[i,t] <- 6.733511                          # mean of mean(logTrapEffort.s, na.rm=T)
      #TrapEffort.sigma[i,t] <-1/(1.028217)^2;               # sd from sd(logTrapEffort.s, na.rm=T) 
      #logTrapEffort.s[i,t] ~  dgamma(64.359377, 9.019594); #dnorm(TrapEffort.mu[i,t], TrapEffort.sigma[i,t]); #logTrap effort follows a normal dist (or trap effort follows a log-normal dist - from hist)
      
      #For Posterior predictive loss
      #Generate replicated data
      #Y.rep[i,t] ~ dpois(exp(log.lambda[i,t])); 
      
      #Compute squared difference b/t observed and replicated data
      #sqdiff[i,t] <- pow(Y[i,t]-Y.rep[i,t], 2);
    }
    
  }
  #For posterior predictive loss, sum squared differences across all years
  #Dsum <- sum(sqdiff);
  
}


#################################################################################
# 1 state model
# one sub-pops
# 
# 
##############################################################################

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
    #tauQ[i] <- tauQ[1];             #diagonal and equal
    #sigmaQ[i] <- 1/sqrt(tauQ[i]);
  }
  
  for(i in 2:n.pop) {
    # uncomment out these next two lines for independent Rs
    #tauR[i] ~ dgamma(0.01,0.01);
    #sigmaR[i] <- 1/sqrt(tauR[i]);
    tauR[i] <- tauR[1];               #diagonal and equal
    sigmaR[i] <- 1/sqrt(tauR[i]);
  }
  
  # PROCESS MODEL
  for(i in 1:n.states.1) {
    x[i,1] ~ dnorm(0,0.1);                    # independent priors on initial states
    U[i] ~ dnorm(0,1);
    B[i] <- 1;                                # setting B=1
    #comment out for density denpendence
    #B[i] ~ dunif(0,1)               #B is b/t 0 and 1. Density depend
    
    
    for(t in 2:n.yrs) {
      predx[i,t] <- B[i]*x[i,t-1] + U[i];       #X = BX + U
      x[i,t] ~ dnorm(predx[i,t],tauQ[i]);
    }
  }
  
  # DATA / LIKELIHOOD MODEL (Obs model)
  A[1] <- 0;
  #A[2] ~ dnorm(0,1);
  
  for(i in 1:n.pop) {
    for(t in 1:n.yrs) {
      # inprod does matrix multiplication
      predy[t] <- inprod(Z.1[i,],x[,t]) + A[i];               #Y=ZX+A        
      lambda[t] ~ dnorm(predy[t], tauR[i]);                #obs error only
      
      Y[t] ~ dnorm(lambda[t], tauR[i]);
      
      #For Posterior predictive loss
      #Generate replicated data
      #Y.rep[i,t] ~ dpois(exp(log.lambda[i,t])); 
      
      #Compute squared difference b/t observed and replicated data
      #sqdiff[i,t] <- pow(Y[i,t]-Y.rep[i,t], 2);
    }
    
  }
  #For posterior predictive loss, sum squared differences across all years
  #Dsum <- sum(sqdiff);
  
}

##################################
# function for bays marss sim
# 1 state model
# no density dependency

bayes.marss.1.fn <- function(y){
  Y <- y
  n.pop <- 1 #dim(Y)[1]
  n.yrs <- nYr #dim(Y)[2]
  whichPop <- 1 #rep(1,n.pop)
  n.states.1 <- 1 #max(whichPop)
  
  Z.1 <- matrix(1) #matrix(0,2,n.states.1)
  #for(i in 1:length(whichPop)) {
  #  Z.1[i,whichPop[i]] = 1
  #}
  
  model.location.1ZC.TN <- MARSS.bayes.1 #MARSS.1.B #model function
  jags.params <- c("x", "sigmaQ","sigmaR", "B", "U", "A") 
  jags.data <- list("Y","n.pop","n.yrs","n.states.1","Z.1")
  
  # Set MCMC parameters
  mcmcchains <- 3
  mcmcthin <- 25
  mcmcburn <- 60000 
  samples2Save <- 40000
  
  jags.parallel(jags.data, inits = NULL, parameters.to.save= jags.params,
       model.file=model.location.1ZC.TN, n.chains = mcmcchains, n.thin = mcmcthin,
       n.burnin=mcmcburn, n.iter =(mcmcburn+samples2Save), DIC = TRUE)
}

#################################################################################
# 1 state model. Split growth U into birth and survival
# 1 sub-pop
# No density dependance
# 
##############################################################################

MARSS.1.bs <- function(){
  
  # PRIORS
  tauQ[1] ~ dgamma(0.01,0.01);
  sigmaQ[1] <- 1/sqrt(tauQ[1]); # convert from precision to sd
  #tauR[1] ~ dgamma(0.01,0.01); #poisson won't run. error in sigmaR
  tauR[1] ~ dgamma(0.03,0.03);
  sigmaR[1] <- 1/sqrt(tauR[1]); # convert from precision to sd
  
  for(i in 2:n.states.1) {
    # uncomment out these next two lines for independent Qs
    tauQ[i] ~ dgamma(0.01,0.01);
    sigmaQ[i] <- 1/sqrt(tauQ[i]);
    #tauQ[i] <- tauQ[1];             #diagonal and equal
    #sigmaQ[i] <- 1/sqrt(tauQ[i]);
  }
  
  for(i in 2:n.pop) {
    # uncomment out these next two lines for independent Rs
    #tauR[i] ~ dgamma(0.01,0.01);
    #sigmaR[i] <- 1/sqrt(tauR[i]);
    tauR[i] <- tauR[1];               #diagonal and equal
    sigmaR[i] <- 1/sqrt(tauR[i]);
  }
  
  # PROCESS MODEL
  for(i in 1:n.states.1) {
    x[i,1] ~ dnorm(0,0.1);                    # independent priors on initial states
    #U[i] ~ dnorm(0,1);
    #s[i] ~ dunif(0,1); # survival prob
    #b[i] ~ dpois(3); # mean births
    b[i] ~ dunif(0,1); # mean birth rate
    #d[i] ~ dpois(1); #deaths
    d[i] ~ dunif(0,1); #death rate
    # density denpendence
    B[i] <- 1;                                # setting B=1
    #comment out for density denpendence
    #B[i] ~ dunif(0,1)               #B is b/t 0 and 1. Density depend
    
    
    for(t in 2:n.yrs) {
      predx[i,t] <- B[i]*x[i,t-1] + b[i]*x[i,t-1]-d[i]*x[i,t-1]; #X = BX + U: U split into b[i] (birth rate) and d (death rate)
      x[i,t] ~ dnorm(predx[i,t],tauQ[i]);
    }
  }
  
  # DATA / LIKELIHOOD MODEL (Obs model)
  A[1] <- 0;
  A[2] ~ dnorm(0,1);
  
  for(i in 1:n.pop) {
    for(t in 1:n.yrs) {
      # inprod does matrix multiplication
      predy[i,t] <- inprod(Z.1[i,],x[,t]) + A[i];               #Y=ZX+A        
      lambda[i,t] ~ dnorm(predy[i,t], tauR[i]);                #obs error only
      
      #comment out next 4 rows for Y~neg bin
      # lambda.nb[i,t] <- exp(predy[i,t])                      #lambda out of log-space
      # p[i,t] <- r[i,t]/(r[i,t]+lambda.nb[i,t])                # calc p, r=size
      # Y[i,t] ~ dnegbin(p[i,t],r[i,t])                         # Y from bin neg dist
      # r[i,t] ~ dunif(0,100)
      
      #comment out row below for Y~Poisson
      # Y[i,t] ~ dpois(exp(predy[i,t]));                      #obs from a poisson dist. exp is putting lambda on Y scale as x in log space
      
      #comment out row below for Y~normal and comment in log.lambda line above
      Y[i,t] ~ dnorm(lambda[i,t], tauR[i]);
      
      
      #For Posterior predictive loss
      #Generate replicated data
      #Y.rep[i,t] ~ dpois(exp(log.lambda[i,t])); 
      
      #Compute squared difference b/t observed and replicated data
      #sqdiff[i,t] <- pow(Y[i,t]-Y.rep[i,t], 2);
    }
    
  }
  #For posterior predictive loss, sum squared differences across all years
  #Dsum <- sum(sqdiff);
  
}


