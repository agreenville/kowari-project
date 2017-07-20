####################################
# Kowari pop dyn and habitat use
#
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

# add year
kowari.yr <- cbind(kowari, Year=format(as.Date(kowari$Obs.Capt.date, format="%d/%m/%Y"),"%Y"))
head(kowari.yr)

# trapping effort
Trap.effort <- read.csv("data/trappingEffort.csv" , header=TRUE)
head(Trap.effort)



# trapping effort and captures per year
trap.effort.yr <- aggregate(Trap.Nights~Year+Site, data=Trap.effort, sum)
kowari.count.yr <- aggregate(Species~Site+Year, data=kowari.yr, length)

head(kowari.count.yr)

kowari.count.grid.yr <- subset(kowari.count.yr, Site=="PAN GRID"|Site=="WAL GRID")

kowari.effort <- merge(kowari.count.grid.yr, trap.effort.yr, by=c("Year", "Site"))
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


years <-unique(Years$Year)
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

#empty B matrix
#B <- matrix( nrow=2,ncol=1)
#B <-NA

model.location.1ZC.TN <- MARSS.1.B #model function


#Then we can run the model using the code
# Specify the parameters / data list / model file name
jags.params <- c("x", "sigmaQ","sigmaR", "B", "U", "A") 
jags.data <- list("Y","n.pop","n.yrs","n.states.1","Z.1") #     


# Set MCMC parameters
mcmcchains <- 3
mcmcthin <- 10
mcmcburn <- 6000 
samples2Save <- 4000
kowariMARSS <- jags(jags.data, inits = NULL, parameters.to.save= jags.params,
                            model.file=model.location.1ZC.TN, n.chains = mcmcchains, n.thin = mcmcthin,
                            n.burnin=mcmcburn, n.iter =(mcmcburn+samples2Save), DIC = TRUE)

# #parallel processing
# system.time(mulgaraMARSS.1.ZC.TN <- jags.parallel(jags.data, inits = NULL, parameters.to.save= jags.params,
#                                                   model.file=model.location.1ZC.TN, n.chains = 3, n.thin = 10,
#                                                   n.burnin=280000, n.iter =300000, DIC = TRUE))

attach.jags(kowariMARSS)

kowariMARSS


#kowari
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
#tiff(filename = "output/fig_mulgara.tiff", width = 200, height = 160, units = 'mm', res = 300) 

# 1 state
par(mfrow = c(1,1), xpd=F, mar=c(5.1, 4.1, 4.1, 11), family = "serif")

matplot(years, (kow.x.ZC.1.TN.r), ylim = c(0,15), type="l",lwd=1, xlab="Year", ylab="Captures (100 trap nights)",
        bty="l", font.lab=2)
polygon(c(years,rev(years)),c(t(kow.x.ZC.1.TN.LC.r), rev( t(kow.x.ZC.1.TN.UC.r))) ,
        col = adjustcolor("grey90", 0.9), border = NA)
matlines(years, (kow.x.ZC.1.TN.r), lwd=1)
matpoints(years,t(kowari.l),pch=15:23, cex =1, col=1)
#mtext("(c)",3, adj=0, line=2)

#mtext("Year",1, adj=0.5, line=3.5, font=2)
#mtext("Captures per 100 trap nights", 2, adj=0.4,line=3.2, font=2)

legend("topright",xpd=T, legend=sites, pch=15:23, cex=1, pt.cex = 1,
       col=c(1,1),box.col=NA,inset=c(-.33,0))

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

#empty B matrix
#B <- matrix( nrow=2,ncol=2)

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
par(mfrow = c(1,1), xpd=F, mar=c(5.1, 4.1, 4.1, 11), family = "serif")

matplot(years, (kow.x.2.TN.r), ylim = c(0,15), type="l",lwd=1, xlab="Year", ylab="Captures (100 trap nights)",
        bty="l", font.lab=2)
  polygon(c(years,rev(years)),c(t(kow.x.2.TN.LC.r), rev( t(kow.x.2.TN.UC.r))) ,
        col = adjustcolor("grey90", 0.9), border = NA)
  polygon(c(years,rev(years)),c(t(kow.x.2.TN.2.LC.r), rev( t(kow.x.2.TN.2.UC.r))) ,
        col = adjustcolor("blue", 0.4), border = NA)

  matlines(years, (kow.x.2.TN.r), lwd=2)
  matlines(years, (kow.x.2.TN.2.r), lwd=2, col="blue")
  matpoints(years,t(kowari.l),pch=15:23, cex =1, col=c("black", "blue"))

  legend("topright",xpd=T, legend=sites, pch=15:23, cex=1, pt.cex = 1,
       col=c("black","blue"),box.col=NA,inset=c(-.33,0))

# dia
  mcmcplot(kowariMARSS.2)
  
  caterplot(kowariMARSS.2, "x")
  caterplot(kowariMARSS.2, "A")
  caterplot(kowariMARSS.2, parms = c("U", "B","sigmaQ", "sigmaR"),
            style= "gray")  


