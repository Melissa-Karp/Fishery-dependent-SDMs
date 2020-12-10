
setwd("~/DisMAP project/Location, Location, Location/Location Workshop")

#----Load Library & Function----
library(mgcv)
library(ggplot2)
library(viridis)
library(raster)
library(gbm)
library(dismo)
source("~/DisMAP project/Location, Location, Location/Location Workshop/SimulatedWorld_ROMS_FishDep_UnequalCoverage.R") #load ROMS simulation function

#Set parameters for functions
dir <- "~/DisMAP project/Location, Location, Location/Location Workshop/ROMS/hadley" #directory where ROMS data is stored (on dropbox, email steph for access)

dat <- SimulateWorld_ROMS_FishDepFun_WAORports(dir=dir, nsamples = 100) #takes a few mins
names(dat)[names(dat) == 'sst'] <- 'temp' #matching roms names. Quick temporary fix.
#head(dat)

#----Explore Data----
#Create dataframe with historical/forecast data
dat_hist_1 <- dat[dat$year<=2010,]
dat_hist<-dat_hist_1[dat_hist_1$year>1984,] ## so dat_hist used to fit the model is 1985-2010
dat_fcast <- dat[dat$year>2010,] #forecast using 2011-2100

dat_hist$log_abundance <- log(dat_hist$abundance)


dat_hist_random<-dat_hist[dat_hist$random_sampled>0,]
dat_hist_Tar<-dat_hist[dat_hist$pref_sampled>0,]
dat_hist_Opt<-dat_hist[dat_hist$opt_sampled>0,]
dat_hist_Dist<-dat_hist[dat_hist$Dist_sampled>0,]
dat_hist_BY<-dat_hist[dat_hist$BY_sampled>0,]
dat_hist_ByDist<-dat_hist[dat_hist$ByDist_sampled>0,]
dat_hist_CA<-dat_hist[dat_hist$Closed_sampled>0,]

dat_fcast_random<-dat_fcast[dat_fcast$random_sampled>0,]
dat_fcast_Tar<-dat_fcast[dat_fcast$pref_sampled>0,]
dat_fcast_Opt<-dat_fcast[dat_fcast$opt_sampled>0,]
dat_fcast_Dist<-dat_fcast[dat_fcast$Dist_sampled>0,]
dat_fcast_BY<-dat_fcast[dat_fcast$BY_sampled>0,]
dat_fcast_ByDist<-dat_fcast[dat_fcast$ByDist_sampled>0,]
head(dat_hist)

##Testing it out with just a GAM first

  #Random sampling
gam_Ran_P <- gam(pres~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_random, family=binomial) #note - chl-surface is log of chl-surface
#plot(gam_Ran_P, pages=1)
gam_Ran_N <- gam(log_abundance~s(temp, bs="gp")+ s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_random[dat_hist_random$abundance>0,], family=gaussian)
#plot(gam_Ran_N, pages=1)

dat_hist$presxR <- predict(gam_Ran_P, dat_hist, type="response")
abundxR <- predict(gam_Ran_N, dat_hist, type="response") 
dat_hist$gam_Ran_Abund <- dat_hist$presxR * exp(abundxR)  #predicted catch

dat_fcast$presxR <- predict(gam_Ran_P, dat_fcast, type="response")
abundxR <- predict(gam_Ran_N, dat_fcast, type="response")
dat_fcast$gam_Ran_Abund <- dat_fcast$presxR * exp(abundxR)

  #Pref sampling
gam_Tar_P <- gam(pres~s(temp, bs="gp")+ s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_Tar, family=binomial)
#plot(gam_Tar_P, pages=1)
gam_Tar_N <- gam(log_abundance~s(temp, bs="gp")+ s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_Tar[dat_hist_Tar$abundance>0,], family=gaussian)
#plot(gam_Tar_N, pages=1)

dat_hist$presxT <- predict(gam_Tar_P, dat_hist, type="response")
abundxT <- predict(gam_Tar_N, dat_hist, type="response") 
dat_hist$gam_Tar_Abund <- dat_hist$presxT * exp(abundxT)

dat_fcast$presxT <- predict(gam_Tar_P, dat_fcast, type="response")
abundxT <- predict(gam_Tar_N, dat_fcast, type="response")
dat_fcast$gam_Tar_Abund <- dat_fcast$presxT * exp(abundxT)

#Opt sampling
gam_Opt_P <- gam(pres~s(temp, bs="gp")+ s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_Opt, family=binomial)
#plot(gam_Tar_P, pages=1)
gam_Opt_N <- gam(log_abundance~s(temp, bs="gp")+ s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_Opt[dat_hist_Opt$abundance>0,], family=gaussian)
#plot(gam_Tar_N, pages=1)

dat_hist$presxO <- predict(gam_Opt_P, dat_hist, type="response")
abundxO <- predict(gam_Opt_N, dat_hist, type="response") 
dat_hist$gam_Opt_Abund <- dat_hist$presxO * exp(abundxO)

dat_fcast$presxO <- predict(gam_Opt_P, dat_fcast, type="response")
abundxO <- predict(gam_Opt_N, dat_fcast, type="response")
dat_fcast$gam_Opt_Abund <- dat_fcast$presxO * exp(abundxO)

  #Dist sampling
gam_Dist_P <- gam(pres~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_Dist, family=binomial)
#plot(gam_E_P, pages=1)
gam_Dist_N <- gam(log_abundance~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_Dist[dat_hist_Dist$abundance>0,], family=gaussian)
#plot(gam_E_N, pages=1)

dat_hist$presxD <- predict(gam_Dist_P, dat_hist, type="response")
abundxD <- predict(gam_Dist_N, dat_hist, type="response") 
dat_hist$gam_Dist_Abund <- dat_hist$presxD * exp(abundxD)

dat_fcast$presxD <- predict(gam_Dist_P, dat_fcast, type="response")
abundxD <- predict(gam_Dist_N, dat_fcast, type="response")
dat_fcast$gam_Dist_Abund <- dat_fcast$presxD * exp(abundxD)

 #BY sampling
gam_BY_P <- gam(pres~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_BY, family=binomial)
#plot(gam_E_P, pages=1)
gam_BY_N <- gam(log_abundance~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_BY[dat_hist_BY$abundance>0,], family=gaussian)
#plot(gam_E_N, pages=1)

dat_hist$presxB <- predict(gam_BY_P, dat_hist, type="response")
abundxB <- predict(gam_BY_N, dat_hist, type="response") 
dat_hist$gam_BY_Abund <- dat_hist$presxB * exp(abundxB)

dat_fcast$presxB <- predict(gam_BY_P, dat_fcast, type="response")
abundxB <- predict(gam_BY_N, dat_fcast, type="response")
dat_fcast$gam_BY_Abund <- dat_fcast$presxB * exp(abundxB)  

  #ByDist sampling
gam_BD_P <- gam(pres~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_ByDist, family=binomial)
#plot(gam_E_P, pages=1)
gam_BD_N <- gam(log_abundance~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_ByDist[dat_hist_ByDist$abundance>0,], family=gaussian)
#plot(gam_E_N, pages=1)

dat_hist$presxBD <- predict(gam_BD_P, dat_hist, type="response")
abundxBD <- predict(gam_BD_N, dat_hist, type="response") 
dat_hist$gam_BD_Abund <- dat_hist$presxBD * exp(abundxBD)

dat_fcast$presxBD <- predict(gam_BD_P, dat_fcast, type="response")
abundxBD <- predict(gam_BD_N, dat_fcast, type="response")
dat_fcast$gam_BD_Abund <- dat_fcast$presxBD * exp(abundxBD)

  #Closed Area Sampling
gam_CA_P <- gam(pres~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_CA, family=binomial)
#plot(gam_E_P, pages=1)
gam_CA_N <- gam(log_abundance~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_CA[dat_hist_CA$abundance>0,], family=gaussian)
#plot(gam_E_N, pages=1)

dat_hist$presxCA <- predict(gam_CA_P, dat_hist, type="response")
abundxCA <- predict(gam_CA_N, dat_hist, type="response") 
dat_hist$gam_CA_Abund <- dat_hist$presxCA * exp(abundxCA)

dat_fcast$presxCA <- predict(gam_CA_P, dat_fcast, type="response")
abundxCA <- predict(gam_CA_N, dat_fcast, type="response")
dat_fcast$gam_CA_Abund <- dat_fcast$presxCA * exp(abundxCA)


###### Boosted Regression Trees (BRTs) ######

#random sampling
brt_R_P <- gbm.step(data=dat_hist_random,
                    gbm.x = 14:16,
                    gbm.y = 'pres',
                    family = "bernoulli",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)

brt_R_N <- gbm.step(data=dat_hist_random[dat_hist_random$abundance>0,], 
                    gbm.x = 14:16,
                    gbm.y = 'log_abundance',
                    family = "gaussian",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)

presRx <- predict(brt_R_P, dat_hist, n.trees=brt_R_P$gbm.call$best.trees, type="response")
abundRx <- exp(predict(brt_R_N, dat_hist, n.trees=brt_R_N$gbm.call$best.trees, type="response"))
dat_hist$brt_Ran <- presRx * abundRx

presRx <- predict(brt_R_P, dat_fcast, n.trees=brt_R_P$gbm.call$best.trees, type="response")
abundRx <- exp(predict(brt_R_N, dat_fcast, n.trees=brt_R_N$gbm.call$best.trees, type="response"))
dat_fcast$brt_Ran <- presRx * abundRx

#Target sampling
brt_T_P <- gbm.step(data=dat_hist_Tar,
                    gbm.x = 14:16,
                    gbm.y = 'pres',
                    family = "bernoulli",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)

brt_T_N <- gbm.step(data=dat_hist_Tar[dat_hist_Tar$abundance>0,], 
                    gbm.x = 14:16,
                    gbm.y = 'log_abundance',
                    family = "gaussian",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)

presTx <- predict(brt_T_P, dat_hist, n.trees=brt_T_P$gbm.call$best.trees, type="response")
abundTx <- exp(predict(brt_T_N, dat_hist, n.trees=brt_T_N$gbm.call$best.trees, type="response"))
dat_hist$brt_Tar <- presTx * abundTx

presTx <- predict(brt_T_P, dat_fcast, n.trees=brt_T_P$gbm.call$best.trees, type="response")
abundTx <- exp(predict(brt_T_N, dat_fcast, n.trees=brt_T_N$gbm.call$best.trees, type="response"))
dat_fcast$brt_Tar <- presTx * abundTx

##Fishing in Optimal habitat suitability (stronger pref for high quality target spp habitat suit)
brt_O_P <- gbm.step(data=dat_hist_Opt,
                    gbm.x = 14:16,
                    gbm.y = 'pres',
                    family = "bernoulli",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)

brt_O_N <- gbm.step(data=dat_hist_Opt[dat_hist_Opt$abundance>0,], 
                    gbm.x = 14:16,
                    gbm.y = 'log_abundance',
                    family = "gaussian",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)

presOx <- predict(brt_O_P, dat_hist, n.trees=brt_O_P$gbm.call$best.trees, type="response")
abundOx <- exp(predict(brt_O_N, dat_hist, n.trees=brt_O_N$gbm.call$best.trees, type="response"))
dat_hist$brt_Opt <- presOx * abundOx

presOx <- predict(brt_O_P, dat_fcast, n.trees=brt_O_P$gbm.call$best.trees, type="response")
abundOx <- exp(predict(brt_O_N, dat_fcast, n.trees=brt_O_N$gbm.call$best.trees, type="response"))
dat_fcast$brt_Opt <- presOx * abundOx

##Distance sampling + Opt Target Spp Habitat Suitability 
brt_D_P <- gbm.step(data=dat_hist_Dist,
                    gbm.x = 14:16,
                    gbm.y = 'pres',
                    family = "bernoulli",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)

brt_D_N <- gbm.step(data=dat_hist_Dist[dat_hist_Dist$abundance>0,], 
                    gbm.x = 14:16,
                    gbm.y = 'log_abundance',
                    family = "gaussian",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)

presDx <- predict(brt_D_P, dat_hist, n.trees=brt_D_P$gbm.call$best.trees, type="response")
abundDx <- exp(predict(brt_D_N, dat_hist, n.trees=brt_D_N$gbm.call$best.trees, type="response"))
dat_hist$brt_Dist <- presDx * abundDx

presDx <- predict(brt_D_P, dat_fcast, n.trees=brt_D_P$gbm.call$best.trees, type="response")
abundDx <- exp(predict(brt_D_N, dat_fcast, n.trees=brt_D_N$gbm.call$best.trees, type="response"))
dat_fcast$brt_Dist <- presDx * abundDx

## Bycatch + Opt Target Sampling
brt_B_P <- gbm.step(data=dat_hist_BY,
                    gbm.x = 14:16,
                    gbm.y = 'pres',
                    family = "bernoulli",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)

brt_B_N <- gbm.step(data=dat_hist_BY[dat_hist_BY$abundance>0,], 
                    gbm.x = 14:16,
                    gbm.y = 'log_abundance',
                    family = "gaussian",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)

presBx <- predict(brt_B_P, dat_hist, n.trees=brt_B_P$gbm.call$best.trees, type="response")
abundBx <- exp(predict(brt_B_N, dat_hist, n.trees=brt_B_N$gbm.call$best.trees, type="response"))
dat_hist$brt_BY <- presBx * abundBx

presBx <- predict(brt_B_P, dat_fcast, n.trees=brt_B_P$gbm.call$best.trees, type="response")
abundBx <- exp(predict(brt_B_N, dat_fcast, n.trees=brt_B_N$gbm.call$best.trees, type="response"))
dat_fcast$brt_BY<- presBx * abundBx

##Bycatch + Distance + Opt Species Habitat Suitability 
brt_BD_P <- gbm.step(data=dat_hist_ByDist,
                    gbm.x = 14:16,
                    gbm.y = 'pres',
                    family = "bernoulli",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)

brt_BD_N <- gbm.step(data=dat_hist_ByDist[dat_hist_ByDist$abundance>0,], 
                    gbm.x = 14:16,
                    gbm.y = 'log_abundance',
                    family = "gaussian",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)

presBDx <- predict(brt_BD_P, dat_hist, n.trees=brt_BD_P$gbm.call$best.trees, type="response")
abundBDx <- exp(predict(brt_BD_N, dat_hist, n.trees=brt_BD_N$gbm.call$best.trees, type="response"))
dat_hist$brt_BD <- presBDx * abundBDx

presBDx <- predict(brt_BD_P, dat_fcast, n.trees=brt_BD_P$gbm.call$best.trees, type="response")
abundBDx <- exp(predict(brt_BD_N, dat_fcast, n.trees=brt_BD_N$gbm.call$best.trees, type="response"))
dat_fcast$brt_BD<- presBDx * abundBDx

##Closed Areas + Opt Target Species
brt_CA_P <- gbm.step(data=dat_hist_CA,
                     gbm.x = 14:16,
                     gbm.y = 'pres',
                     family = "bernoulli",
                     tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                     plot.main=FALSE, verbose = FALSE)

brt_CA_N <- gbm.step(data=dat_hist_CA[dat_hist_CA$abundance>0,], 
                     gbm.x = 14:16,
                     gbm.y = 'log_abundance',
                     family = "gaussian",
                     tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                     plot.main=FALSE, verbose = FALSE)

presCAx <- predict(brt_CA_P, dat_hist, n.trees=brt_CA_P$gbm.call$best.trees, type="response")
abundCAx <- exp(predict(brt_CA_N, dat_hist, n.trees=brt_CA_N$gbm.call$best.trees, type="response"))
dat_hist$brt_CA <- presCAx * abundCAx

presCAx <- predict(brt_CA_P, dat_fcast, n.trees=brt_CA_P$gbm.call$best.trees, type="response")
abundCAx <- exp(predict(brt_CA_N, dat_fcast, n.trees=brt_CA_N$gbm.call$best.trees, type="response"))
dat_fcast$brt_CA<- presCAx * abundCAx


## Compare true Abundance to predictions - GAMS

#combine data from histori and future into one dataframe
dat_fcast$log_abundance <- log(dat_fcast$abundance)

dat_all<-rbind(dat_hist, dat_fcast)

par(mfrow=c(1,1))
plot(aggregate(abundance~year,dat_all,FUN="sum"),type="l",  lwd=2,ylab="Biomass")
lines(aggregate(gam_Ran_Abund~year,dat_all,FUN="sum"),type="l",lwd=2,ylab="Biomass",col="blue")
lines(aggregate(gam_Tar_Abund~year,dat_all,FUN="sum"),type="l",lwd=2,ylab="Biomass",col="green")
lines(aggregate(gam_Opt_Abund~year,dat_all,FUN="sum"),type="l",lwd=2,ylab="Biomass",col="light blue")
lines(aggregate(gam_Dist_Abund~year,dat_all,FUN="sum"),type="l",lwd=2,ylab="Biomass",col="red")
lines(aggregate(gam_BY_Abund~year,dat_all,FUN="sum"),type="l",lwd=2,ylab="Biomass",col="yellow")
lines(aggregate(gam_BD_Abund~year,dat_all,FUN="sum"),type="l",lwd=2,ylab="Biomass",col="purple")
lines(aggregate(gam_CA_Abund~year,dat_all,FUN="sum"),type="l",lwd=2,ylab="Biomass",col="orange")
abline(v = 2010, col="black", lwd=3, lty=2)

##Compare True COG to predictions
cog_hist_lat <- as.data.frame(matrix(NA,nrow=nrow(dat_hist),ncol=9))
colnames(cog_hist_lat) <- c("year","truth", "gam_Ran_Abund", "gam_Tar_Abund","gam_Opt_Abund", "gam_Dist_Abund", "gam_BY_Abund", "gam_BD_Abund", "gam_CA_Abund")
counter=1
for (y in unique(dat_all$year)){
  cog_hist_lat[counter,1] <- y
  cog_hist_lat[counter,2] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$abundance[dat_all$year==y])
  cog_hist_lat[counter,3] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$gam_Ran_Abund[dat_all$year==y])
  cog_hist_lat[counter,4] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$gam_Tar_Abund[dat_all$year==y])
  cog_hist_lat[counter,5] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$gam_Opt_Abund[dat_all$year==y])
  cog_hist_lat[counter,6] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$gam_Dist_Abund[dat_all$year==y])
  cog_hist_lat[counter,7] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$gam_BY_Abund[dat_all$year==y])
  cog_hist_lat[counter,8] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$gam_BD_Abund[dat_all$year==y])
  cog_hist_lat[counter,9] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$gam_CA_Abund[dat_all$year==y])
  counter = counter + 1
}
par(mfrow=c(1,1))
plot(cog_hist_lat$year,cog_hist_lat$truth, type='b', ylim=c(35, 44))
lines(cog_hist_lat$year,cog_hist_lat$gam_Ran_Abund, type='b', col="blue")
lines(cog_hist_lat$year,cog_hist_lat$gam_Tar_Abund, type='b', col="green")
lines(cog_hist_lat$year,cog_hist_lat$gam_Opt_Abund, type='b', col="light blue")
lines(cog_hist_lat$year,cog_hist_lat$gam_Dist_Abund, type='b', col="red")
lines(cog_hist_lat$year,cog_hist_lat$gam_BY_Abund, type='b', col="yellow")
lines(cog_hist_lat$year,cog_hist_lat$gam_BD_Abund, type='b', col="purple")
lines(cog_hist_lat$year,cog_hist_lat$gam_CA_Abund, type='b', col="orange")
abline(v = 2010, col="black", lwd=3, lty=2)

### Distributio Map: Presence (Binomial)
par(mfrow=c(3,7), mar=c(2,3,3,2))
for (yy in c(1985, 2030, 2100)) {
  year_x <- yy
  dat_x <- dat_all[dat_all$year==year_x,]
  
  #dat_x$presxR <- predict(Truth, dat_x, type="response")
  Truth <- rasterFromXYZ(dat_x[,c("lon","lat","suitability_t")])
  plot(Truth, asp=1, main=paste0("Truth, year ", year_x))
  
  #dat_x$presxT <- predict(gam_Tar_P, dat_x, type="response")
  r_pred_Tar <- rasterFromXYZ(dat_x[,c("lon","lat","presxT")])
  plot(r_pred_Tar, asp=1, main=paste0("Predicted, Target, year ", year_x))
  
  #dat_x$presxO <- predict(gam_Opt_P, dat_x, type="response")
  r_pred_Opt <- rasterFromXYZ(dat_x[,c("lon","lat","presxO")])
  plot(r_pred_Opt, asp=1, main=paste0("Predicted, Optimal, year ", year_x))
  
  #dat_x$presxD <- predict(gam_Dist_P, dat_x, type="response")
  r_pred_Dist <- rasterFromXYZ(dat_x[,c("lon","lat","presxD")])
  plot(r_pred_Dist, asp=1, main=paste0("Predicted, Distance, year ", year_x))
  
  #dat_x$presxB <- predict(gam_BY_P, dat_x, type="response")
  r_pred_BY <- rasterFromXYZ(dat_x[,c("lon","lat","presxB")])
  plot(r_pred_BY, asp=1, main=paste0("Predicted, Bycatch, year ", year_x))
  
  #dat_x$presxBD <- predict(gam_BD_P, dat_x, type="response")
  r_pred_BD <- rasterFromXYZ(dat_x[,c("lon","lat","presxBD")])
  plot(r_pred_BD, asp=1, main=paste0("Predicted, All Biases, year ", year_x))
  
  #dat_x$presxCA <- predict(gam_CA_P, dat_x, type="response")
  r_pred_CA <- rasterFromXYZ(dat_x[,c("lon","lat","presxCA")])
  plot(r_pred_CA, asp=1, main=paste0("Predicted, Closed Area, year ", year_x))
  
}


####Map of Distributions:Abundance
par(mfrow=c(3,6), mar=c(2,3,3,2))
for (yy in c(1985, 2030, 2100)) {
  year_x <- yy
  dat_x <- dat_all[dat_all$year==year_x,]
  
  #dat_x$presxR <- predict(gam_Ran_P, dat_x, type="response")
  r_pred_Ran <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Ran_Abund")])
  plot(r_pred_Ran, asp=1, main=paste0("Predicted, Random, year ", year_x))
  
  #dat_x$presxT <- predict(gam_Tar_P, dat_x, type="response")
  r_pred_Tar <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Tar_Abund")])
  plot(r_pred_Tar, asp=1, main=paste0("Predicted, Target, year ", year_x))
  
  #dat_x$presxT <- predict(gam_Tar_P, dat_x, type="response")
  r_pred_Opt <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Opt_Abund")])
  plot(r_pred_Opt, asp=1, main=paste0("Predicted, Optimal, year ", year_x))
  
  #dat_x$presxD <- predict(gam_Dist_P, dat_x, type="response")
  r_pred_Dist <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_Abund")])
  plot(r_pred_Dist, asp=1, main=paste0("Predicted, Distance, year ", year_x))
  
  #dat_x$presxB <- predict(gam_BY_P, dat_x, type="response")
  r_pred_BY <- rasterFromXYZ(dat_x[,c("lon","lat","gam_BY_Abund")])
  plot(r_pred_BY, asp=1, main=paste0("Predicted, Bycatch, year ", year_x))
  
  #dat_x$presxBD <- predict(gam_BD_P, dat_x, type="response")
  r_pred_BD <- rasterFromXYZ(dat_x[,c("lon","lat","gam_BD_Abund")])
  plot(r_pred_BD, asp=1, main=paste0("Predicted, All Biases, year ", year_x))
  
  #dat_x$presxCA <- predict(gam_CA_P, dat_x, type="response")
  r_pred_CA <- rasterFromXYZ(dat_x[,c("lon","lat","gam_CA_Abund")])
  plot(r_pred_CA, asp=1, main=paste0("Predicted, Closed Area, year ", year_x))
  
}

## Compare true abundance to predictions - BRTs

par(mfrow=c(1,1))
plot(aggregate(abundance~year,dat_all,FUN="sum"),type="l",  lwd=2,ylab="Biomass")
lines(aggregate(brt_Ran~year,dat_all,FUN="sum"),type="l",lwd=2,ylab="Biomass",col="blue")
lines(aggregate(brt_Tar~year,dat_all,FUN="sum"),type="l",lwd=2,ylab="Biomass",col="green")
lines(aggregate(brt_Opt~year,dat_all,FUN="sum"),type="l",lwd=2,ylab="Biomass",col="light blue")
lines(aggregate(brt_Dist~year,dat_all,FUN="sum"),type="l",lwd=2,ylab="Biomass",col="red")
lines(aggregate(brt_BY~year,dat_all,FUN="sum"),type="l",lwd=2,ylab="Biomass",col="yellow")
lines(aggregate(brt_BD~year,dat_all,FUN="sum"),type="l",lwd=2,ylab="Biomass",col="purple")
lines(aggregate(brt_CA~year,dat_all,FUN="sum"),type="l",lwd=2,ylab="Biomass",col="orange")
abline(v = 2010, col="black", lwd=3, lty=2)

##Compare True COG to predictions
cog_hist_lat <- as.data.frame(matrix(NA,nrow=nrow(dat_hist),ncol=9))
colnames(cog_hist_lat) <- c("year","truth", "brt_Ran", "brt_Tar","brt_Opt", "brt_Dist", "brt_BY", "brt_BD", "brt_CA")
counter=1
for (y in unique(dat_all$year)){
  cog_hist_lat[counter,1] <- y
  cog_hist_lat[counter,2] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$abundance[dat_all$year==y])
  cog_hist_lat[counter,3] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$brt_Ran[dat_all$year==y])
  cog_hist_lat[counter,4] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$brt_Tar[dat_all$year==y])
  cog_hist_lat[counter,5] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$brt_Opt[dat_all$year==y])
  cog_hist_lat[counter,6] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$brt_Dist[dat_all$year==y])
  cog_hist_lat[counter,7] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$brt_BY[dat_all$year==y])
  cog_hist_lat[counter,8] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$brt_BD[dat_all$year==y])
  cog_hist_lat[counter,9] <- weighted.mean(dat_all$lat[dat_all$year==y],w=dat_all$brt_CA[dat_all$year==y])
  counter = counter + 1
}
par(mfrow=c(1,1))
plot(cog_hist_lat$year,cog_hist_lat$truth, type='b', ylim=c(35, 44))
lines(cog_hist_lat$year,cog_hist_lat$brt_Ran, type='b', col="blue")
lines(cog_hist_lat$year,cog_hist_lat$brt_Tar, type='b', col="green")
lines(cog_hist_lat$year,cog_hist_lat$brt_Opt, type='b', col="light blue")
lines(cog_hist_lat$year,cog_hist_lat$brt_Dist, type='b', col="red")
lines(cog_hist_lat$year,cog_hist_lat$brt_BY, type='b', col="yellow")
lines(cog_hist_lat$year,cog_hist_lat$brt_BD, type='b', col="purple")
lines(cog_hist_lat$year,cog_hist_lat$brt_CA, type='b', col="orange")
abline(v = 2010, col="black", lwd=3, lty=2)

####Map of Distributions:Abundance
par(mfrow=c(3,7), mar=c(2,3,3,2))
for (yy in c(1985, 2010, 2090)) {
  year_x <- yy
  dat_x <- dat_all[dat_all$year==year_x,]
  
  #Truth 
  Truth <- rasterFromXYZ(dat_x[,c("lon","lat","abundance")])
  plot(Truth, asp=1, main=paste0("Truth, year ", year_x))
  
  #dat_x$presxR <- predict(gam_Ran_P, dat_x, type="response")
  r_pred_Ran <- rasterFromXYZ(dat_x[,c("lon","lat","brt_Ran")])
  plot(r_pred_Ran, asp=1, main=paste0("Predicted, Random, year ", year_x))
  
  #dat_x$presxT <- predict(gam_Tar_P, dat_x, type="response")
  r_pred_Tar <- rasterFromXYZ(dat_x[,c("lon","lat","brt_Tar")])
  plot(r_pred_Tar, asp=1, main=paste0("Predicted, Target, year ", year_x))
  
  #dat_x$presxT <- predict(gam_Tar_P, dat_x, type="response")
  r_pred_Opt <- rasterFromXYZ(dat_x[,c("lon","lat","brt_Opt")])
  plot(r_pred_Opt, asp=1, main=paste0("Predicted, Optimal, year ", year_x))

  #dat_x$presxD <- predict(gam_Dist_P, dat_x, type="response")
  r_pred_Dist <- rasterFromXYZ(dat_x[,c("lon","lat","brt_Dist")])
  plot(r_pred_Dist, asp=1, main=paste0("Predicted, Distance, year ", year_x))
  
  #dat_x$presxB <- predict(gam_BY_P, dat_x, type="response")
  r_pred_BY <- rasterFromXYZ(dat_x[,c("lon","lat","brt_BY")])
  plot(r_pred_BY, asp=1, main=paste0("Predicted, Bycatch, year ", year_x))
  
  #dat_x$presxBD <- predict(gam_BD_P, dat_x, type="response")
  r_pred_BD <- rasterFromXYZ(dat_x[,c("lon","lat","brt_BD")])
  plot(r_pred_BD, asp=1, main=paste0("Predicted, All Biases, year ", year_x))
  
  #dat_x$presxCA <- predict(gam_CA_P, dat_x, type="response")
  r_pred_CA <- rasterFromXYZ(dat_x[,c("lon","lat","brt_CA")])
  plot(r_pred_CA, asp=1, main=paste0("Predicted, Closed Area, year ", year_x))
  
}

##### Performance Metrics (e.g. RMSE) ######
write.csv(dat_all, 'EM_FullModel_TrueAbund.csv',row.names = FALSE)

source("~/DisMAP project/Location, Location, Location/Location Workshop/PerformanceMeasures_RSME.R") #load ROMS simulation function

RMSE<-sdm_rmse(dat_hist=dat_hist, dat_fcast=dat_fcast)

###EXPLORING THE DATA###
##average temperature sampled with each sampling scenario
Temp_Ran<-aggregate(temp~year,dat_hist_random,FUN="mean")
Temp_Tar<-aggregate(temp~year,dat_hist_Tar,FUN="mean")
Temp_Opt<-aggregate(temp~year,dat_hist_Opt,FUN="mean")
Temp_Dist<-aggregate(temp~year,dat_hist_Dist,FUN="mean")
Temp_BY<-aggregate(temp~year,dat_hist_BY,FUN="mean")
Temp_BD<-aggregate(temp~year,dat_hist_ByDist,FUN="mean")
Temp_CA<-aggregate(temp~year,dat_hist_CA,FUN="mean")

dat_hist_suit<-dat_hist[dat_hist$pres>0,]

plot(aggregate(temp~year, dat_hist_suit, FUN="max"), type="l",lwd=2,ylab="avg temp of samples", ylim=c(15,19))
lines(aggregate(temp~year,dat_hist_random,FUN="max"),type="l",lwd=2,ylab="avg temp of samples", col="grey")
lines(aggregate(temp~year,dat_hist_Tar,FUN="max"),type="l",lwd=2,ylab="Biomass",col="blue")
lines(aggregate(temp~year,dat_hist_Opt,FUN="max"),type="l",lwd=2,ylab="Biomass",col="green")
lines(aggregate(temp~year,dat_hist_Dist,FUN="max"),type="l",lwd=2,ylab="Biomass",col="light blue")
lines(aggregate(temp~year,dat_hist_BY,FUN="max"),type="l",lwd=2,ylab="Biomass",col="red")
lines(aggregate(temp~year, dat_hist_ByDist,FUN="max"),type="l",lwd=2,ylab="Biomass",col="yellow")
lines(aggregate(temp~year,dat_hist_CA,FUN="max"),type="l",lwd=2,ylab="Biomass",col="purple")

## Total Number of zeros with each sampling design
library(tidyverse)
dat_hist_random%>%
  gather(x, value, 4)%>%
  group_by(x)%>%
  tally(value > 0)
dat_hist_Tar%>%
  gather(x, value, 4)%>%
  group_by(x)%>%
  tally(value > 0)
dat_hist_Dist%>%
  gather(x, value, 4)%>%
  group_by(x)%>%
  tally(value > 0)
dat_hist_BY%>%
  gather(x, value, 4)%>%
  group_by(x)%>%
  tally(value > 0)
dat_hist_ByDist%>%
  gather(x, value, 4)%>%
  group_by(x)%>%
  tally(value > 0)

