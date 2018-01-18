####################################
# Kowari pop dyn 
# MARSS analysis
# Aaron
#####################################

# load packages and graphing functions
require(reshape)
library(plyr)
library(mcmcplots)
library(plotrix)
library(ggplot2)
library(gridExtra)
library(grid)

source("R/MARSS_functions.R") #bayesian models
source("R/summarySE.R")

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
             Site = c(rep("WAL GRID", 16),
                      rep("PAN GRID", 16)))

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
######################
# load rainfall data from Birdsville
# rain mostly from airport, but with 2 values from police station to fill missing values

rain <- read.csv("data/BirdvilleAirportRain.csv", header=T) 
head(rain)

# add 1 yr lag coln
rain$lag.1 <- c(NA, rain$Annual[-length(rain$Annual)])

ggplot(rain, aes(Year, Annual))+
  geom_bar(stat="identity")+
  theme_classic()

# load fractual cover from grids
Pan.frac <- read.csv("data/Pan_Year_frac.csv", header=T) 
head(Pan.frac)

# add 1 yr lag coln
Pan.frac$Green_Mean.lag.1 <- c(NA, Pan.frac$Green_Mean[-length(Pan.frac$Green_Mean)])
Pan.frac$NonGreen_Mean.lag.1 <- c(NA, Pan.frac$NonGreen_Mean[-length(Pan.frac$NonGreen_Mean)])
Pan.frac$Bare_Mean.lag.1 <- c(NA, Pan.frac$Bare_Mean[-length(Pan.frac$Bare_Mean)])

Wal.frac <- read.csv("data/Wal_Year_frac.csv", header=T) 
head(Wal.frac)

# add 1 yr lag coln
Wal.frac$Green_Mean.lag.1 <- c(NA, Wal.frac$Green_Mean[-length(Wal.frac$Green_Mean)])
Wal.frac$NonGreen_Mean.lag.1 <- c(NA, Wal.frac$NonGreen_Mean[-length(Wal.frac$NonGreen_Mean)])
Wal.frac$Bare_Mean.lag.1 <- c(NA, Wal.frac$Bare_Mean[-length(Wal.frac$Bare_Mean)])


# combine both sites together
frac.cover <- rbind(Pan.frac, Wal.frac)

  
######################################################
# data for repro and condition 
head(kowari.grid.no.recap.yr) # capture data with recapts removed
kowari.grid.no.recap.yr["Num_young"][is.na(kowari.grid.no.recap.yr["Num_young"])] <- 0

# mean testes size/year
teste.size <- aggregate(Testes_width~Year,data=kowari.grid.no.recap.yr, mean)
teste.size <- data.frame(Year = as.numeric(as.character(teste.size$Year)),
                                           Testes.width = teste.size$Testes_width)


# mean pouch young/year
py <- aggregate(Num_young~Year,data=subset(kowari.grid.no.recap.yr, SEX=="F"), mean)

py <- summarySE(subset(kowari.grid.no.recap.yr, SEX=="F"),
                measurevar="Num_young", groupvars="Year", na.rm=TRUE)
# convert year to numeric
py["Year"]<-as.numeric(as.character(py$Year))

#py <- data.frame(Year = as.numeric(as.character(py$Year)),
#                         py = py$Num_young)

# % of pop breeding/year

no.breeding.f <- aggregate(Species~Year,data=subset(kowari.grid.no.recap.yr, Num_young>0) , length)

tot.f <- aggregate(Species~Year,data=subset(kowari.grid.no.recap.yr, SEX=="F") , length)

prop.breeding.f <- merge(tot.f,no.breeding.f, by="Year", all=T )
prop.breeding.f[is.na(prop.breeding.f)] <- 0
prop.breeding.f <- data.frame(Year = as.numeric(as.character(prop.breeding.f$Year)),
                              total=prop.breeding.f$Species.x, breeding=prop.breeding.f$Species.y,
                              pro.breed = prop.breeding.f$Species.y/prop.breeding.f$Species.x)

# weights with females with py removed
mass <- aggregate(Weight~Year+SEX,
                  data=subset(kowari.grid.no.recap.yr, Num_young<1), mean)

# removing females with py (no recapts)

kowari.no.py <- subset(kowari.grid.no.recap.yr, Num_young<1)
kowari.cond.m <- subset(kowari.no.py, SEX=="M")
kowari.cond.f <- subset(kowari.no.py, SEX=="F")

################################################################################
# ANALYSIS
#
################################################################################
# body condition
#
################################################################################

plot(kowari.cond.m$Head_length, kowari.cond.m$Weight)
 abline(lm(kowari.cond.m$Weight~kowari.cond.m$Head_length))

plot(kowari.cond.f$Head_length, kowari.cond.f$Weight)
 abline(lm(kowari.cond.f$Weight~kowari.cond.f$Head_length))

# Residual deviations from the linear regression line
body.lm.m <- lm(log(Weight)~log(Head_length), data=kowari.cond.m)
 summary(body.lm.m)
 
body.lm.f <- lm(log(Weight)~log(Head_length), data=kowari.cond.f)
  summary(body.lm.f)
  
# residuals
body.res.m <-  body.lm.m$residuals
body.res.f <-  body.lm.f$residuals
# add residuals back to dataset. Note there is an na in head_length 
kowari.cond.m$res[!is.na(kowari.cond.m$Head_length)] <- body.res.m
kowari.cond.f$res <- body.res.f

#plot(kowari.cond.m$Year, kowari.cond.m$res)
 
#plot(kowari.cond.f$Year, kowari.cond.f$res) 

# summary stats 

cond.m <- summarySE(kowari.cond.m, measurevar="res", #males
          groupvars="Year", na.rm=TRUE)

cond.f <- summarySE(kowari.cond.f, measurevar="res", #females
                    groupvars="Year", na.rm=TRUE)
# convert year to numeric
cond.m["Year"]<-as.numeric(as.character(cond.m$Year))
cond.f["Year"]<-as.numeric(as.character(cond.f$Year))

# all years surveyed data.frame
years.df <- data.frame(Year = years.m)

# add missing years
cond.m <- join(cond.m, years.df, type ="full")
cond.f <- join(cond.f, years.df, type ="full")
py.y <- join(py, years.df, type ="full") # no. of py for graphing later on


male.cond.g <- ggplot(cond.m, aes(Year, res))+
                geom_point()+
                geom_errorbar(aes(ymin=res-ci, ymax=res+ci), width=.1) +
                geom_hline(yintercept = 0, linetype=2) + 
                theme_classic()+
                #scale_x_continuous(breaks = seq(2000, 2015, 1))+
                scale_x_continuous(limits=c(1999, 2016))+
                ggtitle("a)")+
                theme(plot.title = element_text(face="plain"),
                      axis.text.x=element_blank())+
                ylab("Residual value")+ xlab("")


female.cond.g <- ggplot(cond.f, aes(Year, res))+
                  geom_point()+
                  geom_errorbar(aes(ymin=res-ci, ymax=res+ci), width=.1) +
                  geom_hline(yintercept = 0, linetype=2) + 
                  theme_classic()+
                  scale_x_continuous(limits=c(1999, 2016))+
                  ggtitle("b)")+
                  theme(plot.title = element_text(face="plain"),
                        axis.text.x=element_blank())+
                  ylab("Residual value")+ xlab("")

rain.g <- ggplot(subset(rain, Year >= 2000 & Year <=2015), aes(Year, Annual))+
  geom_bar(stat="identity")+
  ggtitle("d)")+
  scale_x_continuous(limits=c(1999, 2016))+
  ylab("Rainfall (mm)")+ xlab("")+
  theme_classic()


#png(filename = "output/bodyCondition.png", width = 120, height = 100, units = 'mm', res = 600) 
# 
# grid.arrange(male.cond.g, female.cond.g, rain.g,
#              left=textGrob("",gp=gpar(fontsize=12),  rot=90),
#              bottom=textGrob("Year", gp=gpar(fontsize=12)))


#dev.off()

####################################
# lm of individual animals boby condition (residuals) and rain and fractual cover
# join rain to kowari data
kowari.cond.f.rain <- join(kowari.cond.f, rain, "Year") # female
kowari.cond.m.rain <- join(kowari.cond.m, rain, "Year") # male

# join fractual cover data

kowari.cond.f.rain.frac <- join(kowari.cond.f.rain, frac.cover, c("Year", "Site")) # female
kowari.cond.m.rain.frac <- join(kowari.cond.m.rain, frac.cover, c("Year", "Site")) # male


female.cond.rain.lm <- lm(res~Annual, data=kowari.cond.f.rain)
  summary(female.cond.rain.lm)

female.cond.rain.lag.lm <- lm(res~lag.1, data=kowari.cond.f.rain)
summary(female.cond.rain.lag.lm)

par(mfrow = c(2,2))
plot(female.cond.rain.lm)
plot(female.cond.rain.lag.lm)
mfrow = c(1,1)

male.cond.rain.lm <- lm(res~Annual, data=kowari.cond.m.rain)
summary(male.cond.rain.lm)

male.cond.rain.lag.lm <- lm(res~lag.1, data=kowari.cond.m.rain)
summary(male.cond.rain.lag.lm)

par(mfrow = c(2,2))
plot(male.cond.rain.lm)
plot(male.cond.rain.lag.lm)
par(mfrow = c(1,1))

# frac cover
female.cond.frac.lm <- lm(res~Green_Mean*NonGreen_Mean, data=kowari.cond.f.rain.frac)
summary(female.cond.frac.lm)

female.cond.frac.lag.lm <- lm(res~Green_Mean.lag.1*NonGreen_Mean.lag.1, data=kowari.cond.f.rain.frac)
summary(female.cond.frac.lag.lm)

par(mfrow = c(2,2))
plot(female.cond.frac.lm)
plot(female.cond.frac.lag.lm)
mfrow = c(1,1)

male.cond.frac.lm <- lm(res~Green_Mean*NonGreen_Mean, data=kowari.cond.m.rain.frac)
summary(male.cond.frac.lm)

male.cond.frac.lag.lm <- lm(res~Green_Mean.lag.1*NonGreen_Mean.lag.1, data=kowari.cond.m.rain.frac)
summary(male.cond.frac.lag.lm)

par(mfrow = c(2,2))
plot(male.cond.frac.lm)
plot(male.cond.frac.lag.lm)
mfrow = c(1,1)



########################################
# Reproductive condition
#
# testes size condition

# plot(kowari.cond.m$Head_length, kowari.cond.m$Testes_width)
# abline(lm(kowari.cond.m$Testes_width~kowari.cond.m$Head_length))

# Residual deviations from the linear regression line
testes.lm <- lm(log(Testes_width)~log(Head_length), data=kowari.cond.m)
summary(testes.lm)

# residuals
testes.res <-  testes.lm$residuals

# add residuals back to dataset. Note there is an na in head_length 
kowari.cond.m$teste.res[!is.na(kowari.cond.m$Head_length)] <- testes.res

#plot(kowari.cond.m$Year, kowari.cond.m$teste.res)

# summary stats 

teste.m <- summarySE(kowari.cond.m, measurevar="teste.res", #males
                    groupvars="Year", na.rm=TRUE)

# convert year to numeric
teste.m["Year"]<-as.numeric(as.character(teste.m$Year))

# add missing years
teste.m <- join(teste.m, years.df, type ="full")

male.teste.g <- ggplot(teste.m, aes(Year, teste.res))+
  geom_point()+
  geom_errorbar(aes(ymin=teste.res-ci, ymax=teste.res+ci), width=.1) +
  geom_hline(yintercept = 0, linetype=2) + 
  theme_classic()+
  scale_x_continuous(limits=c(1999, 2016))+
  ggtitle("c)")+ theme(plot.title = element_text(face="plain"),
  axis.text.x=element_blank())+
  ylab("Residual value")+xlab("")



# body condition and male repro condition fig for paper
#png(filename = "output/fig_bodyCondition.png", width = 100, height = 170, units = 'mm', res = 600) 

grid.arrange(male.cond.g, female.cond.g, male.teste.g, rain.g, nrow=4, 
             left=textGrob("",gp=gpar(fontsize=12),  rot=90),
             bottom=textGrob("Year", gp=gpar(fontsize=12)))

# dev.off()

# lm of individual teste condition (residuals) and rain
# join rain to kowari data
kowari.teste.rain <- join(kowari.cond.m, rain, "Year") # male

# join fractual cover data
kowari.teste.rain.frac <- join(kowari.teste.rain, frac.cover, c("Year", "Site")) # male


teste.rain.lag.lm <- lm(teste.res~lag.1, data=kowari.teste.rain)
summary(teste.rain.lag.lm)

teste.rain.lm <- lm(teste.res~Annual, data=kowari.teste.rain)
summary(teste.rain.lm)


par(mfrow = c(2,2))
plot(teste.rain.lag.lm)
plot(teste.rain.lm)
par(mfrow = c(1,1))

# frac cover
teste.frac.lag.lm <- lm(teste.res~Green_Mean.lag.1*NonGreen_Mean.lag.1,
                        data=kowari.teste.rain.frac)
summary(teste.frac.lag.lm)

teste.frac.lm <- lm(teste.res~Green_Mean*NonGreen_Mean,
                    data=kowari.teste.rain.frac)
summary(teste.frac.lm)


par(mfrow = c(2,2))
plot(teste.frac.lag.lm)
plot(teste.frac.lm)
par(mfrow = c(1,1))




# proportion of females breeding/yr
# join with rainfall
prop.breeding.f.rain <- join(prop.breeding.f, rain, "Year") # 

# binomal regression - proportional odds. Proportion of females with py each year
breed.f.p <- round(prop.breeding.f.rain$pro.breed*100)
breed.f.mirror <- 100-breed.f.p
breed.f <- cbind(breed.f.p, breed.f.mirror)

breed.f.glm = glm(breed.f ~ Annual, family = "quasibinomial", 
                      data = prop.breeding.f.rain)
summary(breed.f.glm)

breed.f.lag.glm = glm(breed.f ~ lag.1, family = "quasibinomial", 
                  data = prop.breeding.f.rain)
summary(breed.f.lag.glm)

par(mfrow = c(2,2))
plot(breed.f.glm)
plot(breed.f.lag.glm)
par(mfrow = c(1,1))

# number of py per female and rainfall
# female data
K.females <- subset(kowari.grid.no.recap.yr, SEX=="F")
# join with rainfall
K.females.rain <- join(K.females, rain, "Year") # 

# join with fractual cover
K.females.rain.frac <- join(K.females.rain, frac.cover, c("Year", "Site"))

plot(K.females.rain$lag.1, K.females.rain$Num_young)
abline(lm(K.females.rain$Num_young~K.females.rain$lag.1))

breed.py.lag.glm = glm(Num_young ~ lag.1,family = "quasipoisson", data = K.females.rain)
summary(breed.py.lag.glm)

breed.py.glm = glm(Num_young ~ Annual, family = "quasipoisson", data = K.females.rain)
summary(breed.py.glm)

par(mfrow = c(2,2))
plot(breed.py.lag.glm)
plot(breed.py.glm)
par(mfrow = c(1,1))

# plot model
new.data <- data.frame(lag.1=seq(min(K.females.rain$lag.1, na.rm=TRUE),
                max(K.females.rain$lag.1,na.rm=TRUE), by=1))
breed.py.fit <-predict(breed.py.lag.glm, newdata=new.data, se=TRUE, type = "response")

plot(K.females.rain$lag.1, K.females.rain$Num_young, bty="l",
     pch=16, xlab="Antecedent annual rainfall (mm)", ylab="Number of pouch young")
  lines(new.data$lag.1, breed.py.fit$fit )
  lines(new.data$lag.1, breed.py.fit$fit+1.96*breed.py.fit$se.fit, lty=2 )
  lines(new.data$lag.1, breed.py.fit$fit-1.96*breed.py.fit$se.fit, lty=2 )

# frac cover

breed.py.lag.frac.glm = glm(Num_young ~ Green_Mean.lag.1*NonGreen_Mean.lag.1,
                            family = "quasipoisson", data = K.females.rain.frac)
summary(breed.py.lag.frac.glm)
  
breed.py.frac.glm = glm(Num_young ~ Green_Mean*NonGreen_Mean,
                        family = "quasipoisson", data = K.females.rain.frac)
summary(breed.py.frac.glm)
  
par(mfrow = c(2,2))
plot(breed.py.lag.frac.glm)
plot(breed.py.frac.glm)
par(mfrow = c(1,1))  
  
  
  
  
  
  
  
###############################################################################
# MARSS temporal and spatial dynamics
#
################################################################################
# kowari 100TN data
# 1-state
################################################################################
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

# save data for HPC
#save(kowari.log, file="data/HPC_kowari_log.Rdata")

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

detach.jags()

###############################################
# adding symbols for sites. Used in paper
###############################################

#png(filename = "output/fig_kowari_1State.png", width = 200, height = 160, units = 'mm', res = 300) 
# 1 state
# par(mfrow = c(1,1), xpd=F, mar=c(5.1, 4.1, 4.1, 11), family = "serif")
# 
# matplot(years.m, (kow.x.ZC.1.TN.r), xlim=c(1999, 2016), ylim = c(0,8), type="l",lwd=1, xlab="Year", ylab="Captures (100 trap nights)",
#         bty="l", font.lab=2)
# polygon(c(years.m,rev(years.m)),c(t(kow.x.ZC.1.TN.LC.r), rev( t(kow.x.ZC.1.TN.UC.r))) ,
#         col = adjustcolor("grey90", 0.9), border = NA)
# matlines(years.m, (kow.x.ZC.1.TN.r), lwd=1)
# matpoints(years.m,t(kowari.l),pch=15:16, cex =1, col=c("black", "blue"))
# 
# par(xpd=T)
# for(i in 1:length(prop.breeding.f$pro.breed)){
#   ifelse(prop.breeding.f$breeding[i]>0, 
#   
#   floating.pie(prop.breeding.f$Year[i],8.5,
#              c(prop.breeding.f$breeding[i],
#                (prop.breeding.f$total[i]- prop.breeding.f$breeding[i])),
#              radius=0.4,col=c("darkgrey","white")),
#   
#   floating.pie(prop.breeding.f$Year[i],8.5,
#                c(prop.breeding.f$breeding[i]+0.01,
#                  (prop.breeding.f$total[i]- prop.breeding.f$breeding[i])),
#                radius=0.4,col=c("darkgrey","white")))
#   
#   
# }
# 
# legend("topright",xpd=T, legend=c("Breeding females", paste(sites)), pch=c(21,15,16),
#        pt.bg="darkgrey",
#        cex=1, pt.cex = c(3.8,1,1),
#        col=c("black","black","blue"),box.col=NA,inset=c(-.33,0))
# 
#par(mfrow = c(1,1), xpd=NA, mar=c(5, 4, 4, 2) + 0.1, family = "")
#######

# fig for female reproduction and long-term dynamics for paper

#png(filename = "output/fig_kowari_1State.png", width = 200, height = 160, units = 'mm', res = 300) 

par(mfrow = c(2,1), xpd=F, mar=c(2, 4.1, 3.1, 11), family = "serif")
plot(py.y$Year, py.y$Num_young, xaxt="n", xlab=NA, ylab="Pouch young",bty="l",
     pch=18, xlim=c(1999, 2016),ylim = c(0,6), font.lab=2, cex=1.5)
  axis(1, at=c(2000,2005,2010,2015) ,labels=NA)
  arrows(py.y$Year, py.y$Num_young-py.y$ci, py.y$Year, py.y$Num_young+py.y$ci,
         length=0.05, angle=90, code=3)

  par(xpd=T)
  for(i in 1:length(prop.breeding.f$pro.breed)){
    ifelse(prop.breeding.f$breeding[i]>0, 
         
         floating.pie(prop.breeding.f$Year[i],7,
                      c(prop.breeding.f$breeding[i],
                        (prop.breeding.f$total[i]- prop.breeding.f$breeding[i])),
                      radius=0.4,col=c("darkgrey","white")),
         
         floating.pie(prop.breeding.f$Year[i],7,    #prop.breeding.f$Year
                      c(prop.breeding.f$breeding[i]+0.01,
                        (prop.breeding.f$total[i]- prop.breeding.f$breeding[i])),
                      radius=0.4,col=c("darkgrey","white")))
  } 
  mtext(side = 3, line = 1, 'a)', adj = -.09,font=2)
  
  legend("topright",xpd=T, legend=c("Breeding females", "Pouch young" ,paste(sites)),
         pch=c(21,18,15,16),
       pt.bg="darkgrey",
       cex=1, pt.cex = c(3.8,1.5,1,1),
       col=c("black","black","black","blue"),box.col=NA,inset=c(-.3,0))  

par(mar=c(4.1, 4.1, 0, 11))
matplot(years.m, (kow.x.ZC.1.TN.r), xlim=c(1999, 2016), ylim = c(0,8), type="l",
        lwd=1, xlab="Year", ylab="Captures (100 trap nights)",
        bty="l", font.lab=2)
  polygon(c(years.m,rev(years.m)),c(t(kow.x.ZC.1.TN.LC.r), rev( t(kow.x.ZC.1.TN.UC.r))) ,
        col = adjustcolor("grey90", 0.9), border = NA)
  matlines(years.m, (kow.x.ZC.1.TN.r), lwd=1)
  matpoints(years.m,t(kowari.l),pch=15:16, cex =1, col=c("black", "blue"))

  mtext(side = 3, line = 1, 'b)', adj = -.09,font=2)
par(mfrow = c(1,1), xpd=NA, mar=c(5, 4, 4, 2) + 0.1, family = "")


# dev.off()




# mcmcplot(kowariMARSS)
# 
# caterplot(kowariMARSS, "x")
# caterplot(kowariMARSS, "A")
# caterplot(kowariMARSS, parms = c("U", "B","sigmaQ", "sigmaR"),
#           style= "gray")

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
  # mcmcplot(kowariMARSS.2)
  # 
  # caterplot(kowariMARSS.2, "x")
  # caterplot(kowariMARSS.2, "A")
  # caterplot(kowariMARSS.2, parms = c("U", "B","sigmaQ", "sigmaR"),
  #           style= "gray")  

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

library("MARSS")
    
kowari.pred <-data.frame(years = years.m, pop = kow.x.ZC.1.TN.r)  
kowari.data <- data.frame(years = years.m,
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
  
#plot(tyrs, Pi, ylim = c(0,1),bty="l", type="l", col="red",
#     xlab = "Years", ylab="Probabilty of 90% decline")



# CSEGtmufigure(N = 50, u = mean.u, s2p = mean.Q)
# CSEGriskfigure(kowari.pred, te = 100,  threshold = 0.1,
#               datalogged = FALSE, CI.method = "hessian", CI.sim = 1000)


# see extinction_risk_script for simulations, bootstrapping and PVA's

############
# MARSS package
# 
# Z=factor(c(1,1))
# kow.marss.1 <- MARSS(kowari.log, model=list(Z=Z))
# MARSSparamCIs(kow.marss.1)
