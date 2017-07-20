
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
library(mcmcplots)



#################################################################################
# 1 state model
# 
# 
# density dependence, B b/w 0 and 1
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
    #B[i] <- 1;                                # setting B=1
    #comment out for density denpendence
    B[i] ~ dunif(0,1)               #B is b/t 0 and 1. Density depend
    
    
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
# density dependence, B b/w 0 and 1
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
    #B[i] <- 1;                      # B =1
    #comment out for density denpendence
    B[i] ~ dunif(0,1)               #B is b/t 0 and 1. Density depend
    
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
