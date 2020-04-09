#This code uses output from SimulatedWorld_ROMS_TrophicInteraction_function to build a GAM
#Species A: distribution and abundance drivn by SST
### SAMPLING BIAS###
#Note only using GFDL for now

#----Directories----
#Set your working directory
setwd("~/DisMAP project/Location, Location, Location/Location Workshop")

#----Load Library & Function----
library(dismo)
library(mgcv)
library(ggplot2)
library(viridis)
library(gbm)
library(BBmisc)
source("~/DisMAP project/Location, Location, Location/Location Workshop/SimulatedWorld_ROMS_SamplingBias.R") #load ROMS simulation function

#Set parameters for functions
dir <- "~/DisMAP project/Location, Location, Location/Location Workshop/ROMS" #directory where ROMS data is stored (on dropbox, email steph for access)

#Run this function
dat <- SimulateWorld_ROMS_SB(dir=dir) #takes a few mins
colnames(dat)[1:2] <- c("Lon","Lat")
names(dat)[names(dat) == 'sst'] <- 'temp' #matching roms names. Quick temporary fix. 
  # head(dat)
  # tail(dat)
#Create dataframe with historical/forecast data
dat_hist <- dat[dat$year<=2020,]
dat_fcast <- dat[dat$year>2020,]

#Make some quick plots to explore the data
#All Years
par(mfrow=c(2,2))
plot(aggregate(suitability~year,dat,FUN="mean"),type="l", lwd=2, ylab="Suitability",col="dark grey")
lines(aggregate(suitability~year,dat[dat$year<=2020,],FUN="mean"),col="blue")
plot(aggregate(abundance~year,dat,FUN="sum"),type="l",  lwd=2,ylab="Abundance", col="dark grey")
lines(aggregate(abundance~year,dat[dat$year<=2020,],FUN="sum"),col="blue")
plot(aggregate(temp~year,dat,FUN="min"),type="l",ylab="Temperature",ylim=c(8,30), col="dark grey")

#----Build GAM Models----
#Run if lognormal response was simulated
dat_hist$log_abundance <- log(dat_hist$abundance)

gam1.p <- gam(pres ~ s(temp,bs='gp'), data=dat_hist, family=binomial)
gam1.a <- gam(log_abundance ~ s(temp,bs='gp'), data=dat_hist[dat_hist$abundance>0,], family=gaussian)
summary(gam1.p)
summary(gam1.a)
plot(gam1.p)
plot(gam1.a)

#Historical predictions
dat_hist$gam1.p <- predict(gam1.p,dat_hist,type='response')
dat_hist$gam1.a <- predict(gam1.a,dat_hist,type="response")
dat_hist$gam1 <- dat_hist$gam1.p*exp(dat_hist$gam1.a)

#Future predictions
dat_fcast$gam1.p <- predict(gam1.p,dat_fcast,type='response')
dat_fcast$gam1.a <- predict(gam1.a,dat_fcast,type="response")
dat_fcast$gam1 <- dat_fcast$gam1.p*exp(dat_fcast$gam1.a)

##----Build BRT Model -----#
#Make sure >1000 trees fitted

#function to extract explained deviance from BRT
dev_eval=function(model_object){
  null <- model_object$self.statistics$mean.null
  res <- model_object$self.statistics$mean.resid
  dev=((null - res)/null)*100 
  return(dev)
}

#Run if lognormal response was simulated
brt1.a <- gbm.fixed(data=dat_hist[dat_hist$abundance>0,], gbm.x = 6 ,gbm.y = 8,family = "gaussian",tree.complexity = 3, learning.rate = 0.01, n.trees=1000, bag.fraction = 0.6)
brt1.p <- gbm.fixed(data=dat_hist, gbm.x = 6,gbm.y = 4,family = "bernoulli",tree.complexity = 3, learning.rate = 0.01, n.trees=1000, bag.fraction = 0.6)
# saveRDS(brt1.a,paste0(Sim1,'BRT_Sim1_lognorm.rds'))
# saveRDS(brt1.p,paste0(Sim1,'BRT_Sim1_binom.rds'))
# brt1.a <- readRDS(paste0(Sim1,'BRT_Sim1_lognorm.rds'))#read in model object if required
# brt1.p <- readRDS(paste0(Sim1,'BRT_Sim1_binom.rds'))
dev_eval(brt1.p)
dev_eval(brt1.a)
plot(brt1.p)
plot(brt1.a)

  #Historiacl predictions 
dat_hist$brt1.p <- predict(brt1.p,dat_hist,n.trees=brt1.p$gbm.call$best.trees,type='response')
dat_hist$brt1.a <- predict(brt1.a,dat_hist,n.trees=brt1.a$gbm.call$best.trees,type='response')
dat_hist$brt1 <- dat_hist$brt1.p*exp(dat_hist$brt1.a)

  #Future predictions
dat_fcast$brt1.p <- predict(brt1.p,dat_fcast,n.trees=brt1.p$gbm.call$best.trees,type='response')
dat_fcast$brt1.a <- predict(brt1.a,dat_fcast,n.trees=brt1.a$gbm.call$best.trees,type='response')
dat_fcast$brt1 <- dat_fcast$brt1.p*exp(dat_fcast$brt1.a)

##----Plot results for present and future and compare performance of models-----##
#Quick plot: historical
plot(aggregate(abundance~year,dat_hist,FUN="sum"),type="l", lwd=2, ylab="Abundance")
lines(aggregate(gam1~year,dat_hist,FUN="sum"),type="l", lwd=2, ylab="Abundance",col="blue")
lines(aggregate(brt1~year, dat_hist, FUN="sum"), type="l", lwd=2, ylab="Abundance", col="green")

#Quick plot: future
plot(aggregate(abundance~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance")
lines(aggregate(gam1~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='blue')
lines(aggregate(brt1~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='green')


#-----Calculate and plot centre of gravity-----

#Run if lognormal response was simulated
#Historical COG
cog_hist_lat <- as.data.frame(matrix(NA,nrow=20,ncol=6))
colnames(cog_hist_lat) <- c("year","truth","gam1","gam1.p","brt1","brt1.p")
counter=1
for (y in 1980:2020){
  cog_hist_lat[counter,1] <- y
  cog_hist_lat[counter,2] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$abundance[dat_hist$year==y])
  cog_hist_lat[counter,3] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$gam1[dat_hist$year==y])
  cog_hist_lat[counter,4] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$gam1.p[dat_hist$year==y])
  cog_hist_lat[counter,5] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$brt1[dat_hist$year==y])
  cog_hist_lat[counter,6] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$brt1.p[dat_hist$year==y])
  counter = counter + 1
}
head(cog_hist_lat)
plot(cog_hist_lat$year,cog_hist_lat$truth, type='b')
lines(cog_hist_lat$year,cog_hist_lat$gam1, type='b', col="blue")
lines(cog_hist_lat$year,cog_hist_lat$brt1, type='b', col="green")

#Future COG
cog_fcast_lat <- as.data.frame(matrix(NA,nrow=80,ncol=6))
colnames(cog_fcast_lat) <- c("year","truth","gam1","gam1.p","brt1","brt1.p")
counter=1
for (y in 2021:2100){
  cog_fcast_lat[counter,1] <- y
  cog_fcast_lat[counter,2] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$abundance[dat_fcast$year==y])
  cog_fcast_lat[counter,3] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$gam1[dat_fcast$year==y])
  cog_fcast_lat[counter,4] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$gam1.p[dat_fcast$year==y])
  cog_fcast_lat[counter,5] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$brt1[dat_fcast$year==y])
  cog_fcast_lat[counter,6] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$brt1.p[dat_fcast$year==y])
  counter = counter + 1
}

head(cog_fcast_lat)
plot(cog_fcast_lat$year,cog_fcast_lat$truth, type='b')
lines(cog_fcast_lat$year,cog_fcast_lat$gam1, type='b', col="blue")
lines(cog_fcast_lat$year,cog_fcast_lat$brt1, type='b', col="green")


#-----Plot Surface Predictions-----
#Future
y=2100
#Truth
par(mfrow=c(1,1))
ggplot(dat_fcast[dat_fcast$year==y,],aes(Lon,Lat))+
  geom_tile(aes(fill=abundance)) +
  theme_classic() +
  ggtitle("Truth")+
  labs(y="Latitude") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust=0.5),
        # axis.text.x=element_blank(),
        # axis.ticks=element_blank(),
        # axis.title.x=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_viridis()

#Gam
plot(gam1)
ggplot(dat_fcast[dat_fcast$year==y,],aes(Lon,Lat))+
  geom_tile(aes(fill=gam1)) +
  theme_classic() +
  ggtitle("GAM")+
  labs(y="Latitude") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust=0.5),
        # axis.text.x=element_blank(),
        # axis.ticks=element_blank(),
        # axis.title.x=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_viridis()



  