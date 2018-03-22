################################################################################
################################################################################
# JAGS Functions for MARSS Bayesian models for kowari dynamics
# 
# 
# Greenville et al. (2018). Dynamics, habitat use and extinction risk of
#     a carnivorous desert marsupial.
################################################################################
################################################################################

# Note need to define n.pop, n.yrs, Z and B matrix first in R.
# Capture data (Y) is in a matrix (2 x no. of years).

# setting up vector for number of sub-populations and years
# n.pop <- dim(Y)[1]
# n.yrs <- dim(Y)[2]

# z-matrix
# whichPop <- c(1,2) # assign sub-pop structure. E.g. 2-state model.
# n.states <- max(whichPop) # number of states/sub-population structures

# Construct Z matrix 
# Z <- matrix(0,n.pop,n.states)
# for(i in 1:length(whichPop)) {
#  Z[i,whichPop[i]] = 1
# }

# Acknowledgements: Ward, E., Scheuerell, M. & Holmes, E. (2015). Fish 507:
# Applied Time Series Analysis.https://catalyst.uw.edu/workspace/fish203/35553/243766 

#Load packages
library(MASS)
library(R2jags)
library(rjags)
library(runjags)

################################################################################
# 1 state model (synchronous model)
# two sub-pops
################################################################################

MARSS.1.B <- function(){
  
  # PRIORS
  # Q (process errors) are not equal for all states
  tauQ[1] ~ dgamma(0.01,0.01);
  sigmaQ[1] <- 1/sqrt(tauQ[1]); # convert from precision to sd
  
  # R (observation error)
  tauR[1] ~ dgamma(0.01,0.01); 
  sigmaR[1] <- 1/sqrt(tauR[1]); # convert from precision to sd
  
  for(i in 2:n.states.1) {
    # independent Qs
    tauQ[i] ~ dgamma(0.01,0.01);
    sigmaQ[i] <- 1/sqrt(tauQ[i]);
    }
  
  for(i in 2:n.pop) {
  # diagonal R terms equal 
    tauR[i] <- tauR[1];               #diagonal and equal
    sigmaR[i] <- 1/sqrt(tauR[i]);
  }
  
  # PROCESS MODEL
  for(i in 1:n.states.1) {
    x[i,1] ~ dnorm(0,0.1);                    # independent priors on initial states
    U[i] ~ dnorm(0,1);
    B[i] <- 1;                                # setting B=1
    
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
      lambda[i,t] ~ dnorm(predy[i,t], tauR[i]);                #obs error only
      
      # Y~normal
      Y[i,t] ~ dnorm(lambda[i,t], tauR[i]);
    }
   }
 }

################################################################################
# 2 state model (Asynchronous model)
#
# two sub-pops
################################################################################

MARSS.2.B <- function(){
  
  # PRIORS
  tauQ[1] ~ dgamma(0.01,0.01);
  sigmaQ[1] <- 1/sqrt(tauQ[1]); # convert from precision to sd
  tauR[1] ~ dgamma(0.01,0.01); 
  sigmaR[1] <- 1/sqrt(tauR[1]); # convert from precision to sd
  
  for(i in 2:n.states.1) {
    # independent Qs
    tauQ[i] ~ dgamma(0.01,0.01);
    sigmaQ[i] <- 1/sqrt(tauQ[i]);
   }
  
  for(i in 2:n.pop) {
  # Diagonal and equal R
    tauR[i] <- tauR[1];               
    sigmaR[i] <- 1/sqrt(tauR[i]);
  }
  
  # PROCESS MODEL
  for(i in 1:n.states.1) {
    x[i,1] ~ dnorm(0,0.1);          # independent priors on initial states
    U[i] ~ dnorm(0,1);
    B[i] <- 1;                      # B =1
    
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
      lambda[i,t] ~ dnorm(predy[i,t], tauR[i]);                #obs error only
     
      # Y~normal 
      Y[i,t] ~ dnorm(lambda[i,t], tauR[i]);
    }
  }
 }

