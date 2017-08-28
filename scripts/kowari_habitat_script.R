####################################
# Kowari pop dyn and habitat use
# habitat analysis
# Aaron
#####################################

# load packages and graphing functions
library(ggplot2)

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
head(kowari.grid.recap.yr)

kowari.count.id <- aggregate(Species~WPT, data=kowari.grid.recap.yr, length)
head(kowari.count.id) #species is total captures not adjusted for effort
#add coln for species present
kowari.count.id <- cbind(kowari.count.id, present=rep(1, length(kowari.count.id$WPT)))

# load habitat data around traps
habitat <- read.csv("data/DEWNRSA_Kowari_trap_habitat.csv", header =TRUE)
head(habitat)

habitat.grid <- subset(habitat, Site=="PAN GRID"|Site=="WAL GRID")


# data clean up

# Habitat_LU_cover	
# cover code	cover description
# 1	< 5%
# 2	5 - 25%
# 3	25 - 50%
# 4	50 - 75%
# 5	> 75%
# 
# Gibber Size	
# 1	pebble <50mm
# 2	cobble 50-250mm
# 3	boulder >250mm

unique(habitat.grid$Gibber.size)
unique(habitat.grid$Gibber.pavement.cover)
unique(habitat.grid$Sand.Mound.number)
unique(habitat.grid$Sand.Mound.cover)
unique(habitat.grid$Sand.Spread.cover)
unique(habitat.grid$Hard.Darainage.Depression.cover)
unique(habitat.grid$Cracking.Clay.Drainage.Depression.cover)
unique(habitat.grid$Dune.cover)
unique(habitat.grid$Site)

# # replace bad values with NA
# habitat[habitat==-1] <- NA
# habitat[habitat==0] <- NA

# merage kowari captures with habitat
kowari.habitat <- merge(habitat.grid, kowari.count.id, by.x="GPS_wpt", by.y="WPT" ,all.x = TRUE)
head(kowari.habitat)
kowari.habitat <- droplevels(kowari.habitat)
kowari.habitat[is.na(kowari.habitat)] <- 0
unique(kowari.habitat$present)
# replace Na with 0
kowari.habitat["present"][is.na(kowari.habitat["present"])] <- 0

kowari.habitat.m <- kowari.habitat[,c("GPS_wpt", "Gibber.size", "Gibber.pavement.cover", "Sand.Mound.number",
                                      "Sand.Mound.cover", "Sand.Spread.cover", "Hard.Darainage.Depression.cover",
                                      "Sand.Mound.cover", "Sand.Spread.cover", "Hard.Darainage.Depression.cover",
                                      "Cracking.Clay.Drainage.Depression.cover", "Dune.cover", "present")]
head(kowari.habitat.m)

# replace NA with 0 for MDS
kowari.habitat.m[is.na(kowari.habitat.m)] <- 0

# remove present coln and set it up as a factor for graphing
kowari.habitat.m <- kowari.habitat.m[,-13]
kowari.habitat.m <- as.matrix(kowari.habitat.m[,-1])
rownames(kowari.habitat.m) <- kowari.habitat$GPS_wpt

kowari.factor <-kowari.habitat$present
kowari.factor.df <-data.frame(kowari.factor, present = ifelse(kowari.factor==0, "Absent", "Present"))

####################################################
# Analysis
# nMDS 

library(vegan)
library(parallel)
set.seed(2)

 
NMDS.hab <- metaMDS(kowari.habitat.m, distance="euclidean",  # Our community-by-species matrix
                 trymax = 20, k=3,  parallel = 6) # The number of reduced dimensions
NMDS.hab
 
plot(NMDS.hab)
stressplot(NMDS.hab)

# ordiplot(NMDS)
# 
 ordiplot(NMDS.hab,type="n")
#ordihull(NMDS.hab,groups=kowari.factor,draw="polygon",col="grey90",label=T)
 orditorp(NMDS.hab,display="species",col="red",air=0.01)

 ordiellipse(NMDS.hab, groups=kowari.factor,label=T, col=c(1:6))
# 
ordiplot(NMDS.hab,type="n")
orditorp(NMDS.hab,display="sites", col="red",air=0.01)
ordiellipse(NMDS.hab, groups=kowari.factor,label=T, col=c(1:2))
# 
# ordiplot(NMDS,type="n")
# orditorp(NMDS,display="species",col="red",air=0.01)
# ordiellipse(NMDS, groups=no.threat.sp.status$no.threats,label=T, col=c(1:6))

# ggplot version
ggplot.MDS.fn <- function(NMDS, df.factor, factor){
  NMDS.df <- data.frame(x=NMDS$point[,1],y=NMDS$point[,2],df.factor)
  ggplot(NMDS.df, aes(x,y, colour= factor))+
    geom_point()+
    theme_classic()+
    theme(legend.title=element_blank())+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))+
    xlab("NMDS1")+ylab("NMDS2")
}


ggplot.MDS.fn(NMDS.hab, kowari.factor.df, kowari.factor.df$present)

#######
# binomal glm
library(MuMIn)

kowari.glm <- glm(present ~ Site + scale(Gibber.size) + scale(Gibber.pavement.cover) +
                    scale(Sand.Mound.number) + scale(Sand.Mound.cover) +
                    scale(Sand.Spread.cover) + scale(Hard.Darainage.Depression.cover)   
                     , data = kowari.habitat, family="binomial", na.action = na.fail)

summary(kowari.glm)

par(mfrow = c(2,2))
plot(kowari.glm)
par(mfrow = c(1,1))


dd<- dredge(kowari.glm, rank = "AICc")
plot(subset(dd, delta < 2))
subset(dd, delta < 2)

top.model <- get.models(dd, subset = delta < 2)

summary(model.avg(top.model))

confint(model.avg(top.model))



