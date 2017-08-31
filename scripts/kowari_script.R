####################################
# Kowari pop dyn and habitat use
# MARSS analysis
# Aaron
#####################################

# load packages and graphing functions

library(ggplot2)
library(grid)
library(gridExtra)
library(psych)
require(reshape)
library(plyr)
library(mcmcplots)

source("R/MARSS_functions.R") #bayesian models

########################
# load data
#######################

# mammal data
kowari <- read.csv("data/DEWNRSA_Kowari_captures.csv", header=T) 
head(kowari)

# add year and month
kowari.yr <- cbind(kowari, Year=format(as.Date(kowari$Obs.Capt.date, format="%d/%m/%Y"),"%Y"),
                   Month=format(as.Date(kowari$Obs.Capt.date, format="%d/%m/%Y"),"%m.%Y"))
head(kowari.yr)

# order by date and site
kowari.yr <- kowari.yr[order(kowari.yr$Year, as.Date(kowari.yr$Obs.Capt.date), kowari.yr$Site),]
# subset grids
kowari.grid.yr <- subset(kowari.yr, Site=="PAN GRID"|Site=="WAL GRID")

# id recapts in the same trip 
kowari.grid.recap.yr <- cbind(kowari.grid.yr, 
                              recapt = ifelse(duplicated(kowari.grid.yr$Month) & 
                                                duplicated(kowari.grid.yr$Microchip)==TRUE, 1,"no" ))
# remove recapts
kowari.grid.no.recap.yr <- subset(kowari.grid.recap.yr, recapt == "no")

# trapping effort
Trap.effort <- read.csv("data/trappingEffort.csv" , header=TRUE)
head(Trap.effort)


# trapping effort and captures per year
trap.effort.yr <- aggregate(Trap.Nights~Year+Site, data=Trap.effort, sum)
kowari.count.yr <- aggregate(Species~Site+Year, data=kowari.grid.no.recap.yr, length)

head(kowari.count.yr)

#kowari.count.grid.yr <- subset(kowari.count.yr, Site=="PAN GRID"|Site=="WAL GRID")

kowari.effort <- merge(kowari.count.yr, trap.effort.yr, by=c("Year", "Site"))
kowari.effort <- cbind(kowari.effort, kowari.TN = kowari.effort$Species/kowari.effort$Trap.Nights*100)
kowari.effort.site <- data.frame(Year=kowari.effort$Year, Site = kowari.effort$Site, kowari.TN=kowari.effort$kowari.TN)

#adding missing years
Years <- data.frame(Year=rep(seq(min(as.numeric(as.character(kowari.effort.site$Year))), 
             max(as.numeric(as.character(kowari.effort.site$Year))),1),2), 
             Site = c(rep("WAL GRID", 12), rep("PAN GRID", 12)))

kowari.effort.site.t <- join(kowari.effort.site, Years, type="full")
# sort by year and site
kowari.effort.site.t <-kowari.effort.site.t[order(kowari.effort.site.t$Year, kowari.effort.site.t$Site),]

#first transpose so site name is col with captures as values.
# Captures per 100 trap night data for graphing

kowari.l <- reshape(kowari.effort.site.t, direction="wide",idvar="Site", timevar = "Year" )
kowari.l <- kowari.l[,2:ncol(kowari.l)] #removes site col
kowari.l <- as.matrix(kowari.l) #converts data frame to matrix
is.numeric(kowari.l) # needs to equal true


years.m <-unique(Years$Year)
sites <- droplevels(unique(kowari.effort.site$Site))

#############################################################################
# kowari 100TN data
# 1-state
########################################################################
whichPop.1 <- c(1,1)
n.states.1 <- max(whichPop.1)

#We can use the vector 'whichPop' to construct our Z matrix (following MARSS notation):
# Use whichPop to construct your Z matrix
Z.1 <- matrix(0,2,n.states.1)
for(i in 1:length(whichPop.1)) {
  Z.1[i,whichPop.1[i]] = 1
}
print(Z.1)


# Y captures/ 100 TN
kowari.log <- log(kowari.l+1)


#z score captures/100 N
#mulgaraTN.mean = apply(mulgaralm.log,1,mean, na.rm=T)
#the.sigma.mulgaraTN = sqrt(apply(mulgaralm.log,1,var,na.rm=TRUE))
#mulgaraTN.z=(mulgaralm.log-mulgaraTN.mean)*(1/the.sigma.mulgaraTN)

#plot(table(rodentTN.z), type="h")


#For alexis and mulgara covariate
#Y.TN.z = rbind(mulgaraTN.z, rodentTN.z)

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
kowariMARSS <- jags(jags.data, inits = NULL, parameters.to.save= jags.params,
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


###############################################
# adding symbols for sites. Used in paper
###############################################

#png(filename = "output/fig_kowari_1State.png", width = 200, height = 160, units = 'mm', res = 300) 
# 1 state
par(mfrow = c(1,1), xpd=F, mar=c(5.1, 4.1, 4.1, 11), family = "serif")

matplot(years.m, (kow.x.ZC.1.TN.r), xlim=c(1999, 2012), ylim = c(0,8), type="l",lwd=1, xlab="Year", ylab="Captures (100 trap nights)",
        bty="l", font.lab=2)
polygon(c(years.m,rev(years.m)),c(t(kow.x.ZC.1.TN.LC.r), rev( t(kow.x.ZC.1.TN.UC.r))) ,
        col = adjustcolor("grey90", 0.9), border = NA)
matlines(years.m, (kow.x.ZC.1.TN.r), lwd=1)
matpoints(years.m,t(kowari.l),pch=15:23, cex =1, col=c("black", "blue"))
#mtext("(c)",3, adj=0, line=2)

#mtext("Year",1, adj=0.5, line=3.5, font=2)
#mtext("Captures per 100 trap nights", 2, adj=0.4,line=3.2, font=2)

legend("topright",xpd=T, legend=sites, pch=15:23, cex=1, pt.cex = 1,
       col=c("black","blue"),box.col=NA,inset=c(-.33,0))

par(mfrow = c(1,1), xpd=NA, mar=c(5, 4, 4, 2) + 0.1, family = "")
#dev.off()


detach.jags()


mcmcplot(kowariMARSS)

caterplot(kowariMARSS, "x")
caterplot(kowariMARSS, "A")
caterplot(kowariMARSS, parms = c("U", "B","sigmaQ", "sigmaR"),
          style= "gray")

#############################################################################
# kowari 100TN data
# 2-state
########################################################################
whichPop.1 <- c(1,2)
n.states.1 <- max(whichPop.1)

#We can use the vector 'whichPop' to construct our Z matrix (following MARSS notation):
# Use whichPop to construct your Z matrix
Z.1 <- matrix(0,2,n.states.1)
for(i in 1:length(whichPop.1)) {
  Z.1[i,whichPop.1[i]] = 1
}
print(Z.1)

model.location.1ZC.TN <- MARSS.2.B #model function

kowariMARSS.2 <- jags(jags.data, inits = NULL, parameters.to.save= jags.params,
                    model.file=model.location.1ZC.TN, n.chains = mcmcchains, n.thin = mcmcthin,
                    n.burnin=mcmcburn, n.iter =(mcmcburn+samples2Save), DIC = TRUE)

kowariMARSS.2

attach.jags(kowariMARSS.2)

kow.x.2.TN <- apply(x[,1,],2,mean) #means for each yr for 1 state mulgara
kow.x.2.TN.2 <- apply(x[,2,],2,mean) #means for each yr for 1 state mulgara

kow.x.2.TN.LC <- apply(x[,1,],2,quantile,0.025)
kow.x.2.TN.UC <- apply(x[,1,],2,quantile,0.975)
kow.x.2.TN.LC.2 <- apply(x[,2,],2,quantile,0.025)
kow.x.2.TN.UC.2 <- apply(x[,2,],2,quantile,0.975)

detach.jags()

#Back-transform z-scoring of captures (note values are in log space still)
# for converting x out of log space (logged captures at start)
kow.x.2.TN.r <-  exp(kow.x.2.TN)-1 
kow.x.2.TN.2.r <-  exp(kow.x.2.TN.2)-1 

kow.x.2.TN.LC.r <- exp(kow.x.2.TN.LC)-1
kow.x.2.TN.UC.r <- exp(kow.x.2.TN.UC)-1
kow.x.2.TN.2.LC.r <- exp(kow.x.2.TN.LC.2)-1
kow.x.2.TN.2.UC.r <- exp(kow.x.2.TN.UC.2)-1

# 2 states

#png(filename = "output/fig_kowari_2State.png", width = 200, height = 160, units = 'mm', res = 300) 
par(mfrow = c(1,1), xpd=F, mar=c(5.1, 4.1, 4.1, 11), family = "serif")

matplot(years.m, (kow.x.2.TN.r), ylim = c(0,15), type="l",lwd=1, xlab="Year", ylab="Captures (100 trap nights)",
        bty="l", font.lab=2)
  polygon(c(years.m,rev(years.m)),c(t(kow.x.2.TN.LC.r), rev( t(kow.x.2.TN.UC.r))) ,
        col = adjustcolor("grey90", 0.9), border = NA)
  polygon(c(years.m,rev(years.m)),c(t(kow.x.2.TN.2.LC.r), rev( t(kow.x.2.TN.2.UC.r))) ,
        col = adjustcolor("blue", 0.2), border = NA)

  matlines(years.m, (kow.x.2.TN.r), lwd=2)
  matlines(years.m, (kow.x.2.TN.2.r), lwd=2, col="blue")
  matpoints(years.m,t(kowari.l),pch=15:23, cex =1, col=c("black", "blue"))

  legend("topright",xpd=T, legend=sites, pch=15:23, cex=1, pt.cex = 1,
       col=c("black","blue"),box.col=NA,inset=c(-.33,0))
par(mfrow = c(1,1), xpd=NA, mar=c(5, 4, 4, 2) + 0.1, family = "")
#dev.off()  
  
  
# diagnositics
  mcmcplot(kowariMARSS.2)
  
  caterplot(kowariMARSS.2, "x")
  caterplot(kowariMARSS.2, "A")
  caterplot(kowariMARSS.2, parms = c("U", "B","sigmaQ", "sigmaR"),
            style= "gray")  

############################################################
# PVA
# The probability of hitting a threshold, denoted P(xd; te), is typically presented
# as a curve showing the probabilities of hitting the threshold (y-axis)
# over diferent time horizons (te) on the x-axis. Extinction probabilities can be
# computed through Monte Carlo simulations or analytically using Equation 16
# in Dennis et al. (1991). From MARSS user guide PVA chapter.
  
# xe is the threshold and is dened as xe = log(N0=Ne). N0 is the current
# population estimate and Ne is the threshold. If we are using fractional declines
# then xe =log(N0=(pdN0))=􀀀log(pd). p(u) is the probability that the threshold
# is eventually hit  

#library("PVAClone")  
library("MARSS")
    
kowari.pred <-data.frame(years = years, pop = kow.x.ZC.1.TN.r)  
kowari.data <- data.frame(years = years,
                          kowari.TN = t(kowari.l)[,2])
  
mean.u <- kowariMARSS$BUGSoutput$summary[4] 
mean.Q <- kowariMARSS$BUGSoutput$summary[6]
mean.R <- kowariMARSS$BUGSoutput$summary[7]
Pi <-NA
    
pd = 0.1 #means a 90 percent decline
tyrs = 1:100
xd = -log(pd)
p.ever = ifelse(mean.u<=0,1,exp(-2*mean.u*xd/mean.Q)) #Q=sigma2
for (i in 1:100){
  Pi[i] = p.ever * pnorm((-xd+abs(mean.u)*tyrs[i])/sqrt(mean.Q*tyrs[i]))+
    exp(2*xd*abs(mean.u)/mean.Q)*pnorm((-xd-abs(mean.u)*tyrs[i])/sqrt(mean.Q*tyrs[i]))
  }  
  
plot(tyrs, Pi, ylim = c(0,1),bty="l", type="l", col="red",
     xlab = "Years", ylab="Probabilty of 90% decline")


CSEGtmufigure(N = 50, u = mean.u, s2p = mean.Q)
CSEGriskfigure(kowari.pred, te = 100,  threshold = 0.1, 
               datalogged = FALSE, CI.method = "hessian", CI.sim = 1000)


# m1 <- pva(kowari.data$kowari.TN, gompertz("normal"), c(5,10))
# summary(m1)
# plot(m1)
# coef(m1)
