
setwd("~/DisMAP project/Location, Location, Location/Location Workshop")

#----Load Library & Function----
library(mgcv)
library(ggplot2)
library(viridis)
library(raster)
library(gbm)
library(dismo)
source("~/DisMAP project/Location, Location, Location/Location Workshop/SimulatedWorld_ROMS_FishDep_Final.R") #load ROMS simulation function

#Set parameters for functions
dir <- "~/DisMAP project/Location, Location, Location/Location Workshop/ROMS/hadley" #directory where ROMS data is stored (on dropbox, email steph for access)

#dat <- SimulateWorld_ROMS_FishDepFun_Final(dir=dir, nsamples = 100) #takes a few mins
names(dat)[names(dat) == 'sst'] <- 'temp' #matching roms names. Quick temporary fix.
#head(dat)
dat<-read.csv('FisheryDependent_OM_Simulation_Final.csv', header=T, sep=",")

#----Explore Data----
#Create dataframe with historical/forecast data
dat_hist_1 <- dat[dat$year<=2010,]
dat_hist<-dat_hist_1[dat_hist_1$year>1984,] ## so dat_hist used to fit the model is 1985-2010
dat_fcast <- dat[dat$year>2010,] #forecast using 2011-2100

dat_hist$log_abundance <- log(dat_hist$abundance)

## still need to update code after this point....12/14/20

dat_hist_random<-dat_hist[dat_hist$random_sampled>0,]
dat_hist_Tar<-dat_hist[dat_hist$pref_sampled>0,]
dat_hist_Opt<-dat_hist[dat_hist$opt_sampled>0,]
dat_hist_Dist_npo<-dat_hist[dat_hist$dist_sampled_npo>0,]
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

##### FULL COVARIATE MODELS ######

##GAM w/ Chl

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

#   #Pref sampling
# gam_Tar_P <- gam(pres~s(temp, bs="gp")+ s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_Tar, family=binomial)
# #plot(gam_Tar_P, pages=1)
# gam_Tar_N <- gam(log_abundance~s(temp, bs="gp")+ s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_Tar[dat_hist_Tar$abundance>0,], family=gaussian)
# #plot(gam_Tar_N, pages=1)
# 
# dat_hist$presxT <- predict(gam_Tar_P, dat_hist, type="response")
# abundxT <- predict(gam_Tar_N, dat_hist, type="response") 
# dat_hist$gam_Tar_Abund <- dat_hist$presxT * exp(abundxT)
# 
# dat_fcast$presxT <- predict(gam_Tar_P, dat_fcast, type="response")
# abundxT <- predict(gam_Tar_N, dat_fcast, type="response")
# dat_fcast$gam_Tar_Abund <- dat_fcast$presxT * exp(abundxT)

#Opt sampling
gam_Opt_P <- gam(pres~s(temp, bs="gp")+ s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_Opt, family=binomial)
#plot(gam_Opt_P, pages=1)
gam_Opt_N <- gam(log_abundance~s(temp, bs="gp")+ s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_Opt[dat_hist_Opt$abundance>0,], family=gaussian)
#plot(gam_Opt_N, pages=1)

dat_hist$presxO <- predict(gam_Opt_P, dat_hist, type="response")
abundxO <- predict(gam_Opt_N, dat_hist, type="response") 
dat_hist$gam_Opt_Abund <- dat_hist$presxO * exp(abundxO)

dat_fcast$presxO <- predict(gam_Opt_P, dat_fcast, type="response")
abundxO <- predict(gam_Opt_N, dat_fcast, type="response")
dat_fcast$gam_Opt_Abund <- dat_fcast$presxO * exp(abundxO)

  #Dist sampling
gam_Dist_P <- gam(pres~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_Dist, family=binomial)
#plot(gam_Dist_P, pages=1)
gam_Dist_N <- gam(log_abundance~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_Dist[dat_hist_Dist$abundance>0,], family=gaussian)
#plot(gam_Dist_N, pages=1)

dat_hist$presxD <- predict(gam_Dist_P, dat_hist, type="response")
abundxD <- predict(gam_Dist_N, dat_hist, type="response") 
dat_hist$gam_Dist_Abund <- dat_hist$presxD * exp(abundxD)

dat_fcast$presxD <- predict(gam_Dist_P, dat_fcast, type="response")
abundxD <- predict(gam_Dist_N, dat_fcast, type="response")
dat_fcast$gam_Dist_Abund <- dat_fcast$presxD * exp(abundxD)

 #BY sampling
gam_BY_P <- gam(pres~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_BY, family=binomial)
#plot(gam_BY_P, pages=1)
gam_BY_N <- gam(log_abundance~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_BY[dat_hist_BY$abundance>0,], family=gaussian)
#plot(gam_BY_N, pages=1)

dat_hist$presxB <- predict(gam_BY_P, dat_hist, type="response")
abundxB <- predict(gam_BY_N, dat_hist, type="response") 
dat_hist$gam_BY_Abund <- dat_hist$presxB * exp(abundxB)

dat_fcast$presxB <- predict(gam_BY_P, dat_fcast, type="response")
abundxB <- predict(gam_BY_N, dat_fcast, type="response")
dat_fcast$gam_BY_Abund <- dat_fcast$presxB * exp(abundxB)  

#   #ByDist sampling
# gam_BD_P <- gam(pres~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_ByDist, family=binomial)
# #plot(gam_BD_P, pages=1)
# gam_BD_N <- gam(log_abundance~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_ByDist[dat_hist_ByDist$abundance>0,], family=gaussian)
# #plot(gam_BD_N, pages=1)
# 
# dat_hist$presxBD <- predict(gam_BD_P, dat_hist, type="response")
# abundxBD <- predict(gam_BD_N, dat_hist, type="response") 
# dat_hist$gam_BD_Abund <- dat_hist$presxBD * exp(abundxBD)
# 
# dat_fcast$presxBD <- predict(gam_BD_P, dat_fcast, type="response")
# abundxBD <- predict(gam_BD_N, dat_fcast, type="response")
# dat_fcast$gam_BD_Abund <- dat_fcast$presxBD * exp(abundxBD)

  #Closed Area Sampling
gam_CA_P <- gam(pres~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_CA, family=binomial)
#plot(gam_CA_P, pages=1)
gam_CA_N <- gam(log_abundance~s(temp, bs="gp") + s(mld, bs="gp") + s(chl_surface, bs="gp"), data=dat_hist_CA[dat_hist_CA$abundance>0,], family=gaussian)
#plot(gam_CA_N, pages=1)

dat_hist$presxCA <- predict(gam_CA_P, dat_hist, type="response")
abundxCA <- predict(gam_CA_N, dat_hist, type="response") 
dat_hist$gam_CA_Abund <- dat_hist$presxCA * exp(abundxCA)

dat_fcast$presxCA <- predict(gam_CA_P, dat_fcast, type="response")
abundxCA <- predict(gam_CA_N, dat_fcast, type="response")
dat_fcast$gam_CA_Abund <- dat_fcast$presxCA * exp(abundxCA)


###### Boosted Regression Trees (BRTs) ######

#random sampling
brt_R_P <- gbm.step(data=dat_hist_random,
                    gbm.x = c(25, 27, 28),
                    gbm.y = 'pres',
                    family = "bernoulli",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_R_P, write.title=F, main="brt_Ran_P", plot.layout = c(2,2))

brt_R_N <- gbm.step(data=dat_hist_random[dat_hist_random$abundance>0,], 
                    gbm.x = c(25, 27, 28),
                    gbm.y = 'log_abundance',
                    family = "gaussian",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_R_N, write.title=F, main="brt_Ran_N", plot.layout = c(2,2))


presRx <- predict(brt_R_P, dat_hist, n.trees=brt_R_P$gbm.call$best.trees, type="response")
abundRx <- exp(predict(brt_R_N, dat_hist, n.trees=brt_R_N$gbm.call$best.trees, type="response"))
dat_hist$brt_Ran <- presRx * abundRx

presRx <- predict(brt_R_P, dat_fcast, n.trees=brt_R_P$gbm.call$best.trees, type="response")
abundRx <- exp(predict(brt_R_N, dat_fcast, n.trees=brt_R_N$gbm.call$best.trees, type="response"))
dat_fcast$brt_Ran <- presRx * abundRx

#Target sampling
# brt_T_P <- gbm.step(data=dat_hist_Tar,
#                     gbm.x = c(14,16,17),
#                     gbm.y = 'pres',
#                     family = "bernoulli",
#                     tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
#                     plot.main=FALSE, verbose = FALSE)
# 
# brt_T_N <- gbm.step(data=dat_hist_Tar[dat_hist_Tar$abundance>0,], 
#                     gbm.x = c(14,16,17),
#                     gbm.y = 'log_abundance',
#                     family = "gaussian",
#                     tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
#                     plot.main=FALSE, verbose = FALSE)
# 
# presTx <- predict(brt_T_P, dat_hist, n.trees=brt_T_P$gbm.call$best.trees, type="response")
# abundTx <- exp(predict(brt_T_N, dat_hist, n.trees=brt_T_N$gbm.call$best.trees, type="response"))
# dat_hist$brt_Tar <- presTx * abundTx
# 
# presTx <- predict(brt_T_P, dat_fcast, n.trees=brt_T_P$gbm.call$best.trees, type="response")
# abundTx <- exp(predict(brt_T_N, dat_fcast, n.trees=brt_T_N$gbm.call$best.trees, type="response"))
# dat_fcast$brt_Tar <- presTx * abundTx

##Fishing in Optimal habitat suitability (stronger pref for high quality target spp habitat suit)
brt_O_P <- gbm.step(data=dat_hist_Opt,
                    gbm.x = c(14,16,17),
                    gbm.y = 'pres',
                    family = "bernoulli",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_O_P, write.title=F, main="brt_Opt_P", plot.layout = c(2,2))

brt_O_N <- gbm.step(data=dat_hist_Opt[dat_hist_Opt$abundance>0,], 
                    gbm.x = c(14,16,17),
                    gbm.y = 'log_abundance',
                    family = "gaussian",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_O_N, write.title=F, main="brt_Opt_N", plot.layout = c(2,2))

presOx <- predict(brt_O_P, dat_hist, n.trees=brt_O_P$gbm.call$best.trees, type="response")
abundOx <- exp(predict(brt_O_N, dat_hist, n.trees=brt_O_N$gbm.call$best.trees, type="response"))
dat_hist$brt_Opt <- presOx * abundOx

presOx <- predict(brt_O_P, dat_fcast, n.trees=brt_O_P$gbm.call$best.trees, type="response")
abundOx <- exp(predict(brt_O_N, dat_fcast, n.trees=brt_O_N$gbm.call$best.trees, type="response"))
dat_fcast$brt_Opt <- presOx * abundOx

##Distance sampling + Opt Target Spp Habitat Suitability 
brt_D_P <- gbm.step(data=dat_hist_Dist,
                    gbm.x = c(14,16,17),
                    gbm.y = 'pres',
                    family = "bernoulli",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_D_P, write.title=F, main="brt_Dist_P", plot.layout = c(2,2))

brt_D_N <- gbm.step(data=dat_hist_Dist[dat_hist_Dist$abundance>0,], 
                    gbm.x = c(14,16,17),
                    gbm.y = 'log_abundance',
                    family = "gaussian",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_D_N, write.title=F, main="brt_Dist_N", plot.layout = c(2,2))


presDx <- predict(brt_D_P, dat_hist, n.trees=brt_D_P$gbm.call$best.trees, type="response")
abundDx <- exp(predict(brt_D_N, dat_hist, n.trees=brt_D_N$gbm.call$best.trees, type="response"))
dat_hist$brt_Dist <- presDx * abundDx

presDx <- predict(brt_D_P, dat_fcast, n.trees=brt_D_P$gbm.call$best.trees, type="response")
abundDx <- exp(predict(brt_D_N, dat_fcast, n.trees=brt_D_N$gbm.call$best.trees, type="response"))
dat_fcast$brt_Dist <- presDx * abundDx

## Bycatch + Opt Target Sampling
brt_B_P <- gbm.step(data=dat_hist_BY,
                    gbm.x = c(14,16,17),
                    gbm.y = 'pres',
                    family = "bernoulli",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_B_P, write.title=F, main="brt_BY_P", plot.layout = c(2,2))

brt_B_N <- gbm.step(data=dat_hist_BY[dat_hist_BY$abundance>0,], 
                    gbm.x = c(14,16,17),
                    gbm.y = 'log_abundance',
                    family = "gaussian",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_B_N, write.title=F, main="brt_BY_N", plot.layout = c(2,2))

presBx <- predict(brt_B_P, dat_hist, n.trees=brt_B_P$gbm.call$best.trees, type="response")
abundBx <- exp(predict(brt_B_N, dat_hist, n.trees=brt_B_N$gbm.call$best.trees, type="response"))
dat_hist$brt_BY <- presBx * abundBx

presBx <- predict(brt_B_P, dat_fcast, n.trees=brt_B_P$gbm.call$best.trees, type="response")
abundBx <- exp(predict(brt_B_N, dat_fcast, n.trees=brt_B_N$gbm.call$best.trees, type="response"))
dat_fcast$brt_BY<- presBx * abundBx

##Bycatch + Distance + Opt Species Habitat Suitability 
# brt_BD_P <- gbm.step(data=dat_hist_ByDist,
#                     gbm.x = c(14,16,17),
#                     gbm.y = 'pres',
#                     family = "bernoulli",
#                     tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
#                     plot.main=FALSE, verbose = FALSE)
# 
# brt_BD_N <- gbm.step(data=dat_hist_ByDist[dat_hist_ByDist$abundance>0,], 
#                     gbm.x = c(14,16,17),
#                     gbm.y = 'log_abundance',
#                     family = "gaussian",
#                     tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
#                     plot.main=FALSE, verbose = FALSE)
# 
# presBDx <- predict(brt_BD_P, dat_hist, n.trees=brt_BD_P$gbm.call$best.trees, type="response")
# abundBDx <- exp(predict(brt_BD_N, dat_hist, n.trees=brt_BD_N$gbm.call$best.trees, type="response"))
# dat_hist$brt_BD <- presBDx * abundBDx
# 
# presBDx <- predict(brt_BD_P, dat_fcast, n.trees=brt_BD_P$gbm.call$best.trees, type="response")
# abundBDx <- exp(predict(brt_BD_N, dat_fcast, n.trees=brt_BD_N$gbm.call$best.trees, type="response"))
# dat_fcast$brt_BD<- presBDx * abundBDx

##Closed Areas + Opt Target Species
brt_CA_P <- gbm.step(data=dat_hist_CA,
                     gbm.x = c(14,16,17),
                     gbm.y = 'pres',
                     family = "bernoulli",
                     tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                     plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_CA_P, write.title=F, main="brt_CA_P", plot.layout = c(2,2))

brt_CA_N <- gbm.step(data=dat_hist_CA[dat_hist_CA$abundance>0,], 
                     gbm.x = c(14,16,17),
                     gbm.y = 'log_abundance',
                     family = "gaussian",
                     tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                     plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_CA_N, write.title=F, main="brt_CA_N", plot.layout = c(2,2))

presCAx <- predict(brt_CA_P, dat_hist, n.trees=brt_CA_P$gbm.call$best.trees, type="response")
abundCAx <- exp(predict(brt_CA_N, dat_hist, n.trees=brt_CA_N$gbm.call$best.trees, type="response"))
dat_hist$brt_CA <- presCAx * abundCAx

presCAx <- predict(brt_CA_P, dat_fcast, n.trees=brt_CA_P$gbm.call$best.trees, type="response")
abundCAx <- exp(predict(brt_CA_N, dat_fcast, n.trees=brt_CA_N$gbm.call$best.trees, type="response"))
dat_fcast$brt_CA<- presCAx * abundCAx


######### MISSING COVARIATE MODELS __ Models withoug chl-a #########
#Random sampling
gam_Ran_P_nochl <- gam(pres~s(temp, bs="gp") + s(mld, bs="gp"), data=dat_hist_random, family=binomial)
#plot(gam_Ran_P_nochl, pages=1)
gam_Ran_N_nochl <- gam(log_abundance~s(temp, bs="gp")+ s(mld, bs="gp"), data=dat_hist_random[dat_hist_random$abundance>0,], family=gaussian)
#plot(gam_Ran_N_nochl, pages=1)

dat_hist$presxR_nochl <- predict(gam_Ran_P_nochl, dat_hist, type="response")
abundxR_nochl <- predict(gam_Ran_N_nochl, dat_hist, type="response") 
dat_hist$gam_Ran_Abund_nochl <- dat_hist$presxR_nochl * exp(abundxR_nochl)  #predicted catch

dat_fcast$presxR_nochl <- predict(gam_Ran_P_nochl, dat_fcast, type="response")
abundxR_nochl <- predict(gam_Ran_N_nochl, dat_fcast, type="response")
dat_fcast$gam_Ran_Abund_nochl <- dat_fcast$presxR_nochl * exp(abundxR_nochl)

# #Pref sampling
# gam_Tar_P_nochl <- gam(pres~s(temp, bs="gp")+ s(mld, bs="gp"), data=dat_hist_Tar, family=binomial)
# #plot(gam_Tar_P, pages=1)
# gam_Tar_N_nochl <- gam(log_abundance~s(temp, bs="gp")+ s(mld, bs="gp"), data=dat_hist_Tar[dat_hist_Tar$abundance>0,], family=gaussian)
# #plot(gam_Tar_N, pages=1)
# 
# dat_hist$presxT_nochl <- predict(gam_Tar_P_nochl, dat_hist, type="response")
# abundxT_nochl <- predict(gam_Tar_N_nochl, dat_hist, type="response") 
# dat_hist$gam_Tar_Abund_nochl <- dat_hist$presxT_nochl * exp(abundxT_nochl)
# 
# dat_fcast$presxT_nochl <- predict(gam_Tar_P_nochl, dat_fcast, type="response")
# abundxT_nochl <- predict(gam_Tar_N_nochl, dat_fcast, type="response")
# dat_fcast$gam_Tar_Abund_nochl <- dat_fcast$presxT_nochl * exp(abundxT_nochl)

#Opt sampling
gam_Opt_P_nochl <- gam(pres~s(temp, bs="gp")+ s(mld, bs="gp"), data=dat_hist_Opt, family=binomial)
#plot(gam_Opt_P_nochl, pages=1)
gam_Opt_N_nochl <- gam(log_abundance~s(temp, bs="gp")+ s(mld, bs="gp"), data=dat_hist_Opt[dat_hist_Opt$abundance>0,], family=gaussian)
#plot(gam_Opt_N_nochl, pages=1)

dat_hist$presxO_nochl <- predict(gam_Opt_P_nochl, dat_hist, type="response")
abundxO_nochl <- predict(gam_Opt_N_nochl, dat_hist, type="response") 
dat_hist$gam_Opt_Abund_nochl <- dat_hist$presxO_nochl * exp(abundxO_nochl)

dat_fcast$presxO_nochl <- predict(gam_Opt_P_nochl, dat_fcast, type="response")
abundxO_nochl <- predict(gam_Opt_N_nochl, dat_fcast, type="response")
dat_fcast$gam_Opt_Abund_nochl <- dat_fcast$presxO_nochl * exp(abundxO_nochl)

#Dist sampling
gam_Dist_P_nochl <- gam(pres~s(temp, bs="gp") + s(mld, bs="gp"), data=dat_hist_Dist, family=binomial)
#plot(gam_Dist_P_nochl, pages=1)
gam_Dist_N_nochl <- gam(log_abundance~s(temp, bs="gp") + s(mld, bs="gp"), data=dat_hist_Dist[dat_hist_Dist$abundance>0,], family=gaussian)
#plot(gam_Dist_N_nochl, pages=1)

dat_hist$presxD_nochl <- predict(gam_Dist_P_nochl, dat_hist, type="response")
abundxD_nochl <- predict(gam_Dist_N_nochl, dat_hist, type="response") 
dat_hist$gam_Dist_Abund_nochl <- dat_hist$presxD_nochl * exp(abundxD_nochl)

dat_fcast$presxD_nochl <- predict(gam_Dist_P_nochl, dat_fcast, type="response")
abundxD_nochl <- predict(gam_Dist_N_nochl, dat_fcast, type="response")
dat_fcast$gam_Dist_Abund_nochl <- dat_fcast$presxD_nochl * exp(abundxD_nochl)

#BY sampling
gam_BY_P_nochl <- gam(pres~s(temp, bs="gp") + s(mld, bs="gp"), data=dat_hist_BY, family=binomial)
#plot(gam_BY_P_nochl, pages=1)
gam_BY_N_nochl <- gam(log_abundance~s(temp, bs="gp") + s(mld, bs="gp"), data=dat_hist_BY[dat_hist_BY$abundance>0,], family=gaussian)
#plot(gam_BY_N_nochl, pages=1)

dat_hist$presxB_nochl <- predict(gam_BY_P_nochl, dat_hist, type="response")
abundxB_nochl <- predict(gam_BY_N_nochl, dat_hist, type="response") 
dat_hist$gam_BY_Abund_nochl <- dat_hist$presxB_nochl * exp(abundxB_nochl)

dat_fcast$presxB_nochl <- predict(gam_BY_P_nochl, dat_fcast, type="response")
abundxB_nochl <- predict(gam_BY_N_nochl, dat_fcast, type="response")
dat_fcast$gam_BY_Abund_nochl <- dat_fcast$presxB_nochl * exp(abundxB_nochl)  

# #ByDist sampling
# gam_BD_P_nochl <- gam(pres~s(temp, bs="gp") + s(mld, bs="gp"), data=dat_hist_ByDist, family=binomial)
# #plot(gam_E_P, pages=1)
# gam_BD_N_nochl <- gam(log_abundance~s(temp, bs="gp") + s(mld, bs="gp"), data=dat_hist_ByDist[dat_hist_ByDist$abundance>0,], family=gaussian)
# #plot(gam_E_N, pages=1)
# 
# dat_hist$presxBD_nochl <- predict(gam_BD_P_nochl, dat_hist, type="response")
# abundxBD_nochl <- predict(gam_BD_N_nochl, dat_hist, type="response") 
# dat_hist$gam_BD_Abund_nochl <- dat_hist$presxBD_nochl * exp(abundxBD_nochl)
# 
# dat_fcast$presxBD_nochl <- predict(gam_BD_P_nochl, dat_fcast, type="response")
# abundxBD_nochl <- predict(gam_BD_N_nochl, dat_fcast, type="response")
# dat_fcast$gam_BD_Abund_nochl <- dat_fcast$presxBD_nochl * exp(abundxBD_nochl)

#Closed Area Sampling
gam_CA_P_nochl <- gam(pres~s(temp, bs="gp") + s(mld, bs="gp"), data=dat_hist_CA, family=binomial)
#plot(gam_CA_P_nochl, pages=1)
gam_CA_N_nochl <- gam(log_abundance~s(temp, bs="gp") + s(mld, bs="gp"), data=dat_hist_CA[dat_hist_CA$abundance>0,], family=gaussian)
#plot(gam_CA_N_nochl, pages=1)

dat_hist$presxCA_nochl <- predict(gam_CA_P_nochl, dat_hist, type="response")
abundxCA_nochl <- predict(gam_CA_N_nochl, dat_hist, type="response") 
dat_hist$gam_CA_Abund_nochl <- dat_hist$presxCA_nochl * exp(abundxCA_nochl)

dat_fcast$presxCA_nochl <- predict(gam_CA_P_nochl, dat_fcast, type="response")
abundxCA_nochl <- predict(gam_CA_N_nochl, dat_fcast, type="response")
dat_fcast$gam_CA_Abund_nochl <- dat_fcast$presxCA_nochl * exp(abundxCA_nochl)

#BRTs - no Chlorophyll
###### Boosted Regression Trees (BRTs) ######

#random sampling
brt_R_P_nochl <- gbm.step(data=dat_hist_random,
                    gbm.x = c(14,16),
                    gbm.y = 'pres',
                    family = "bernoulli",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_R_P_nochl, write.title=F, main="brt_Ran_P_nochl", plot.layout = c(2,2))

brt_R_N_nochl <- gbm.step(data=dat_hist_random[dat_hist_random$abundance>0,], 
                    gbm.x = c(14,16),
                    gbm.y = 'log_abundance',
                    family = "gaussian",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_R_N_nochl, write.title=F, main="brt_Ran_N_nochl", plot.layout = c(2,2))

presRx_nochl <- predict(brt_R_P_nochl, dat_hist, n.trees=brt_R_P_nochl$gbm.call$best.trees, type="response")
abundRx_nochl <- exp(predict(brt_R_N_nochl, dat_hist, n.trees=brt_R_N_nochl$gbm.call$best.trees, type="response"))
dat_hist$brt_Ran_nochl <- presRx_nochl * abundRx_nochl

presRx_nochl <- predict(brt_R_P_nochl, dat_fcast, n.trees=brt_R_P_nochl$gbm.call$best.trees, type="response")
abundRx_nochl <- exp(predict(brt_R_N_nochl, dat_fcast, n.trees=brt_R_N_nochl$gbm.call$best.trees, type="response"))
dat_fcast$brt_Ran_nochl <- presRx_nochl * abundRx_nochl

#Target sampling
# brt_T_P_nochl <- gbm.step(data=dat_hist_Tar,
#                     gbm.x = c(14,16),
#                     gbm.y = 'pres',
#                     family = "bernoulli",
#                     tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
#                     plot.main=FALSE, verbose = FALSE)
# 
# brt_T_N_nochl <- gbm.step(data=dat_hist_Tar[dat_hist_Tar$abundance>0,], 
#                     gbm.x = c(14,16),
#                     gbm.y = 'log_abundance',
#                     family = "gaussian",
#                     tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
#                     plot.main=FALSE, verbose = FALSE)
# 
# presTx_nochl <- predict(brt_T_P_nochl, dat_hist, n.trees=brt_T_P_nochl$gbm.call$best.trees, type="response")
# abundTx_nochl <- exp(predict(brt_T_N_nochl, dat_hist, n.trees=brt_T_N_nochl$gbm.call$best.trees, type="response"))
# dat_hist$brt_Tar_nochl <- presTx_nochl * abundTx_nochl
# 
# presTx_nochl <- predict(brt_T_P_nochl, dat_fcast, n.trees=brt_T_P_nochl$gbm.call$best.trees, type="response")
# abundTx_nochl <- exp(predict(brt_T_N_nochl, dat_fcast, n.trees=brt_T_N_nochl$gbm.call$best.trees, type="response"))
# dat_fcast$brt_Tar_nochl <- presTx * abundTx_nochl

##Fishing in Optimal habitat suitability (stronger pref for high quality target spp habitat suit)
brt_O_P_nochl <- gbm.step(data=dat_hist_Opt,
                    gbm.x = c(14,16),
                    gbm.y = 'pres',
                    family = "bernoulli",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_O_P_nochl, write.title=F, main="brt_Opt_P_nochl", plot.layout = c(2,2))

brt_O_N_nochl <- gbm.step(data=dat_hist_Opt[dat_hist_Opt$abundance>0,], 
                    gbm.x = c(14,16),
                    gbm.y = 'log_abundance',
                    family = "gaussian",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_R_N_nochl, write.title=F, main="brt_Opt_N_nochl", plot.layout = c(2,2))

presOx_nochl <- predict(brt_O_P_nochl, dat_hist, n.trees=brt_O_P_nochl$gbm.call$best.trees, type="response")
abundOx_nochl <- exp(predict(brt_O_N_nochl, dat_hist, n.trees=brt_O_N_nochl$gbm.call$best.trees, type="response"))
dat_hist$brt_Opt_nochl <- presOx_nochl * abundOx_nochl

presOx_nochl <- predict(brt_O_P_nochl, dat_fcast, n.trees=brt_O_P_nochl$gbm.call$best.trees, type="response")
abundOx_nochl <- exp(predict(brt_O_N_nochl, dat_fcast, n.trees=brt_O_N_nochl$gbm.call$best.trees, type="response"))
dat_fcast$brt_Opt_nochl <- presOx_nochl * abundOx_nochl

##Distance sampling + Opt Target Spp Habitat Suitability 
brt_D_P_nochl <- gbm.step(data=dat_hist_Dist_npo,
                    gbm.x = c(25, 27),
                    gbm.y = 'pres',
                    family = "bernoulli",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_D_P_nochl, write.title=F, main="brt_Dist_P_nochl", plot.layout = c(2,2))

brt_D_N_nochl <- gbm.step(data=dat_hist_Dist_npo[dat_hist_Dist_npo$abundance>0,], 
                    gbm.x = c(25, 27),
                    gbm.y = 'log_abundance',
                    family = "gaussian",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_D_N_nochl, write.title=F, main="brt_Dist_N_nochl", plot.layout = c(2,2))

presDx_nochl <- predict(brt_D_P_nochl, dat_hist, n.trees=brt_D_P_nochl$gbm.call$best.trees, type="response")
abundDx_nochl <- exp(predict(brt_D_N_nochl, dat_hist, n.trees=brt_D_N_nochl$gbm.call$best.trees, type="response"))
dat_hist$brt_Dist_npo_nochl <- presDx_nochl * abundDx_nochl

presDx_nochl <- predict(brt_D_P_nochl, dat_fcast, n.trees=brt_D_P_nochl$gbm.call$best.trees, type="response")
abundDx_nochl <- exp(predict(brt_D_N_nochl, dat_fcast, n.trees=brt_D_N_nochl$gbm.call$best.trees, type="response"))
dat_fcast$brt_Dist_npo_nochl <- presDx_nochl * abundDx_nochl

## Bycatch + Opt Target Sampling
brt_B_P_nochl <- gbm.step(data=dat_hist_BY,
                    gbm.x = c(14,16),
                    gbm.y = 'pres',
                    family = "bernoulli",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_B_P_nochl, write.title=F, main="brt_BY_P_nochl", plot.layout = c(2,2))

brt_B_N_nochl <- gbm.step(data=dat_hist_BY[dat_hist_BY$abundance>0,], 
                    gbm.x = c(14,16),
                    gbm.y = 'log_abundance',
                    family = "gaussian",
                    tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                    plot.main=FALSE, verbose = FALSE)
gbm.plot(brt_B_N_nochl, write.title=F, main="brt_BY_N_nochl", plot.layout = c(2,2))


presBx_nochl <- predict(brt_B_P_nochl, dat_hist, n.trees=brt_B_P_nochl$gbm.call$best.trees, type="response")
abundBx_nochl <- exp(predict(brt_B_N_nochl, dat_hist, n.trees=brt_B_N_nochl$gbm.call$best.trees, type="response"))
dat_hist$brt_BY_nochl <- presBx_nochl * abundBx_nochl

presBx_nochl <- predict(brt_B_P_nochl, dat_fcast, n.trees=brt_B_P_nochl$gbm.call$best.trees, type="response")
abundBx_nochl <- exp(predict(brt_B_N_nochl, dat_fcast, n.trees=brt_B_N_nochl$gbm.call$best.trees, type="response"))
dat_fcast$brt_BY_nochl<- presBx_nochl * abundBx_nochl

##Bycatch + Distance + Opt Species Habitat Suitability 
# brt_BD_P_nochl <- gbm.step(data=dat_hist_ByDist,
#                      gbm.x = c(14,16),
#                      gbm.y = 'pres',
#                      family = "bernoulli",
#                      tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
#                      plot.main=FALSE, verbose = FALSE)
# gbm.plot(brt_CA_P_nochl, write.title=F, main="brt_CA_P_nochl", plot.layout = c(2,2))
# 
# brt_BD_N_nochl <- gbm.step(data=dat_hist_ByDist[dat_hist_ByDist$abundance>0,], 
#                      gbm.x = c(14,16),
#                      gbm.y = 'log_abundance',
#                      family = "gaussian",
#                      tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
#                      plot.main=FALSE, verbose = FALSE)
# gbm.plot(brt_CA_N_nochl, write.title=F, main="brt_CA_N_nochl", plot.layout = c(2,2))
# 
# presBDx_nochl <- predict(brt_BD_P_nochl, dat_hist, n.trees=brt_BD_P_nochl$gbm.call$best.trees, type="response")
# abundBDx_nochl <- exp(predict(brt_BD_N_nochl, dat_hist, n.trees=brt_BD_N_nochl$gbm.call$best.trees, type="response"))
# dat_hist$brt_BD_nochl <- presBDx_nochl * abundBDx_nochl
# 
# presBDx_nochl <- predict(brt_BD_P_nochl, dat_fcast, n.trees=brt_BD_P_nochl$gbm.call$best.trees, type="response")
# abundBDx_nochl <- exp(predict(brt_BD_N_nochl, dat_fcast, n.trees=brt_BD_N_nochl$gbm.call$best.trees, type="response"))
# dat_fcast$brt_BD_nochl<- presBDx_nochl * abundBDx_nochl

##Closed Areas + Opt Target Species
brt_CA_P_nochl <- gbm.step(data=dat_hist_CA,
                     gbm.x = c(14,16),
                     gbm.y = 'pres',
                     family = "bernoulli",
                     tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                     plot.main=FALSE, verbose = FALSE)

brt_CA_N_nochl <- gbm.step(data=dat_hist_CA[dat_hist_CA$abundance>0,], 
                     gbm.x = c(14,16),
                     gbm.y = 'log_abundance',
                     family = "gaussian",
                     tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                     plot.main=FALSE, verbose = FALSE)

presCAx_nochl <- predict(brt_CA_P_nochl, dat_hist, n.trees=brt_CA_P_nochl$gbm.call$best.trees, type="response")
abundCAx_nochl <- exp(predict(brt_CA_N_nochl, dat_hist, n.trees=brt_CA_N_nochl$gbm.call$best.trees, type="response"))
dat_hist$brt_CA_nochl <- presCAx_nochl * abundCAx_nochl

presCAx_nochl <- predict(brt_CA_P_nochl, dat_fcast, n.trees=brt_CA_P_nochl$gbm.call$best.trees, type="response")
abundCAx_nochl <- exp(predict(brt_CA_N_nochl, dat_fcast, n.trees=brt_CA_N_nochl$gbm.call$best.trees, type="response"))
dat_fcast$brt_CA_nochl<- presCAx_nochl * abundCAx_nochl


#################### PERFORMANCE MEASURES #############################

#### This section runs RMSE for abundance, COG, and correlation, and spatial RMSE##
##RMSE - Abundance
source("~/DisMAP project/Location, Location, Location/Location Workshop/PerformanceMeasures_RSME.R") 
RMSE<-sdm_rmse(dat_hist=dat_hist, dat_fcast=dat_fcast)
RMSE
par(mfrow=c(1,1), mar=c(8,5,2,2))
rmseall<-data.frame(RMSE$rmse_abund_all_year)
barplot(rmseall$rmse_hist~rmseall$model, main="Average RMSE Historic", ylab="RMSE", xlab="",cex.names=0.7, las=2)
barplot(rmseall$rmse_fcast~rmseall$model, main="Average RMSE Forecast", ylab="RMSE", xlab="",cex.names=0.7, las=2)

#RMSE - COG
source("~/DisMAP project/Location, Location, Location/Location Workshop/RMSE_CenterGravity.R") 
RMSE_cog<-sdm_cog(dat_hist=dat_hist, dat_fcast=dat_fcast)
RMSE_cog
rmsecog_y<-data.frame(RMSE_cog$cog_lat_fcast_y)

#Correlation 
source("~/DisMAP project/Location, Location, Location/Location Workshop/Performance_COR.R") 
Cor_sdm<-sdm_cor(dat_hist=dat_hist, dat_fcast=dat_fcast)

#spatial RMSE
source("~/DisMAP project/Location, Location, Location/Location Workshop/spatial_RMSE.R") 
SRMSE<-sdm_sqErr(dat_hist=dat_hist, dat_fcast=dat_fcast)
SEH<-data.frame(SRMSE$spatial_err_all_hist)
# group by cell(i.e., lat, lon) and find RMSE for each grid cell
SRMSE_hist<-SEH%>%
group_by(lat, lon) %>% 
 summarise(meanErr=mean(brt_Ran))

SRMSE_hist_map<-ggplot(data = SRMSE_hist, aes(x=lon,y=lat))+
  geom_tile(aes(fill=meanErr))+
  theme_classic() +  labs(y="", x="") +
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradient2() +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48))


####Additional ways to plot RMSE through time
par(mfrow=c(3,5), mar=c(1,4,2,1))
#RMSE Abundance
#Random
plot(gam_Ran_Abund~year, rmseall, type="l", lwd=2, col="blue", ylim=c(5, 11), ylab="RMSE - Abundace", main="Random")
lines(gam_Ran_Abund_nochl~year, rmseall, type='l', lwd=2, col="light blue")
lines(brt_Ran~year, rmseall, type='l', lwd=2, col="red")
lines(brt_Ran_nochl~year, rmseall, type='l', lwd=2, col="pink")
#optimum
plot(gam_Opt_Abund~year, rmseall, type="l", lwd=2, col="blue", ylim=c(5, 11), ylab="",main="Target")
lines(gam_Opt_Abund_nochl~year, rmseall, type='l', lwd=2, col="light blue")
lines(brt_Opt~year, rmseall, type='l', lwd=2, col="red")
lines(brt_Opt_nochl~year, rmseall, type='l', lwd=2, col="pink")
#Distance
plot(gam_Dist_Abund~year, rmseall, type="l", lwd=2, col="blue", ylim=c(5, 30),ylab="", main="Distance")
lines(gam_Dist_Abund_nochl~year, rmseall, type='l', lwd=2, col="light blue")
lines(brt_Dist~year, rmseall, type='l', lwd=2, col="red")
lines(brt_Dist_nochl~year, rmseall, type='l', lwd=2, col="pink")
#Bycatch
plot(gam_BY_Abund~year, rmseall, type="l", lwd=2, col="blue", ylim=c(5, 11), ylab="", main="Bycatch")
lines(gam_BY_Abund_nochl~year, rmseall, type='l', lwd=2, col="light blue")
lines(brt_BY~year, rmseall, type='l', lwd=2, col="red")
lines(brt_BY_nochl~year, rmseall, type='l', lwd=2, col="pink")
#Closed Area
plot(gam_CA_Abund~year, rmseall, type="l", lwd=2, col="blue", ylim=c(5, 11), ylab="", main="Closed Area")
lines(gam_CA_Abund_nochl~year, rmseall, type='l', lwd=2, col="light blue")
lines(brt_CA~year, rmseall, type='l', lwd=2, col="red")
lines(brt_CA_nochl~year, rmseall, type='l', lwd=2, col="pink")

#RMSE Center of Gravity
#Random
plot(gam_Ran_Abund~year, rmsecog_y, type="l", lwd=2, col="blue", ylab="COG: Obs - Pred")
lines(gam_Ran_Abund_nochl~year, rmsecog_y, type='l', lwd=2, col="light blue")
lines(brt_Ran~year, rmsecog_y, type='l', lwd=2, col="red")
lines(brt_Ran_nochl~year, rmsecog_y, type='l', lwd=2, col="pink")
#optimum
plot(gam_Opt_Abund~year, rmsecog_y, type="l", lwd=2, col="blue", ylab="")
lines(gam_Opt_Abund_nochl~year, rmsecog_y, type='l', lwd=2, col="light blue")
lines(brt_Opt~year, rmsecog_y, type='l', lwd=2, col="red")
lines(brt_Opt_nochl~year, rmsecog_y, type='l', lwd=2, col="pink")
#Distance
plot(gam_Dist_Abund~year, rmsecog_y, type="l", lwd=2, col="blue", ylab="", ylim=c(-2, 8))
lines(gam_Dist_Abund_nochl~year, rmsecog_y, type='l', lwd=2, col="light blue")
lines(brt_Dist~year, rmsecog_y, type='l', lwd=2, col="red")
lines(brt_Dist_nochl~year, rmsecog_y, type='l', lwd=2, col="pink")
#Bycatch
plot(gam_BY_Abund~year, rmsecog_y, type="l", lwd=2, col="blue", ylab="")
lines(gam_BY_Abund_nochl~year, rmsecog_y, type='l', lwd=2, col="light blue")
lines(brt_BY~year, rmsecog_y, type='l', lwd=2, col="red")
lines(brt_BY_nochl~year, rmsecog_y, type='l', lwd=2, col="pink")
#Closed Area
plot(gam_CA_Abund~year, rmsecog_y, type="l", lwd=2, col="blue", ylab="")
lines(gam_CA_Abund_nochl~year, rmsecog_y, type='l', lwd=2, col="light blue")
lines(brt_CA~year, rmsecog_y, type='l', lwd=2, col="red")
lines(brt_CA_nochl~year, rmsecog_y, type='l', lwd=2, col="pink")

###EXPLORING THE DATA###
##average temperature sampled with each sampling scenario
dat_suit<-dat_all[dat_all$pres>0,]

#max temp
plot(aggregate(temp~year, dat_suit, FUN="max"), type="l",lwd=2,col="gray",ylab="max temp of samples", ylim=c(8,25), main="temp sampled")
lines(aggregate(temp~year,dat_hist_random,FUN="max"),type="l",lwd=2, col="light blue")
#lines(aggregate(temp~year,dat_hist_Tar,FUN="max"),type="l",lwd=2,ylab="Biomass",col="green")
lines(aggregate(temp~year,dat_hist_Opt,FUN="max"),type="l",lwd=2, col="green")
lines(aggregate(temp~year,dat_hist_Dist,FUN="max"),type="l",lwd=2,col="pink")
lines(aggregate(temp~year,dat_hist_BY,FUN="max"),type="l",lwd=2,col="yellow")
#lines(aggregate(temp~year, dat_hist_ByDist,FUN="max"),type="l",lwd=2,ylab="Biomass",col="purple")
lines(aggregate(temp~year,dat_hist_CA,FUN="max"),type="l",lwd=2,col="orange")

#min temp
lines(aggregate(temp~year, dat_suit, FUN="min"), type="l",lwd=2,ylab="min temp of samples",col="gray", ylim=c(9,21), main="min temp sampled")
lines(aggregate(temp~year,dat_hist_random,FUN="min"),type="l",lwd=2,ylab="avg temp of samples", col="light blue")
#lines(aggregate(temp~year,dat_hist_Tar,FUN="min"),type="l",lwd=2,ylab="Biomass",col="green")
lines(aggregate(temp~year,dat_hist_Opt,FUN="min"),type="l",lwd=2,ylab="Biomass",col="green")
lines(aggregate(temp~year,dat_hist_Dist,FUN="min"),type="l",lwd=2,ylab="Biomass",col="pink")
lines(aggregate(temp~year,dat_hist_BY,FUN="min"),type="l",lwd=2,ylab="Biomass",col="yellow")
#lines(aggregate(temp~year, dat_hist_ByDist,FUN="min"),type="l",lwd=2,ylab="Biomass",col="purple")
lines(aggregate(temp~year,dat_hist_CA,FUN="min"),type="l",lwd=2,ylab="Biomass",col="orange")

#average Temp
lines(aggregate(temp~year, dat_suit, FUN="mean"), type="l",lwd=2,ylab="avg temp of samples", ylim=c(11,19), main="Average temp sampled")
lines(aggregate(temp~year,dat_hist_random,FUN="mean"),type="l",lwd=2,ylab="avg temp of samples", col="blue")
#lines(aggregate(temp~year,dat_hist_Tar,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="green")
lines(aggregate(temp~year,dat_hist_Opt,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="dark green")
lines(aggregate(temp~year,dat_hist_Dist,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="red")
lines(aggregate(temp~year,dat_hist_BY,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="yellow")
#lines(aggregate(temp~year, dat_hist_ByDist,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="purple")
lines(aggregate(temp~year,dat_hist_CA,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="orange")

#Chl-a
#max chla
plot(aggregate(chl_surface~year, dat_suit, FUN="max"), type="l",lwd=2,col="gray",ylab="max chl of samples", main="chla")
lines(aggregate(chl_surface~year,dat_hist_random,FUN="max"),type="l",lwd=2, col="blue")
#lines(aggregate(temp~year,dat_hist_Tar,FUN="max"),type="l",lwd=2,ylab="Biomass",col="green")
lines(aggregate(chl_surface~year,dat_hist_Opt,FUN="max"),type="l",lwd=2, col="green")
lines(aggregate(chl_surface~year,dat_hist_Dist,FUN="max"),type="l",lwd=2,col="red")
lines(aggregate(chl_surface~year,dat_hist_BY,FUN="max"),type="l",lwd=2,col="yellow")
#lines(aggregate(temp~year, dat_hist_ByDist,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="purple")
lines(aggregate(chl_surface~year,dat_hist_CA,FUN="max"),type="l",lwd=2,col="orange")

#MLD
plot(aggregate(mld~year, dat_suit, FUN="mean"), type="l",lwd=2,col="gray",ylab="mld", main="average mld")
lines(aggregate(mld~year,dat_hist_random,FUN="mean"),type="l",lwd=2, col="blue")
#lines(aggregate(temp~year,dat_hist_Tar,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="green")
lines(aggregate(mld~year,dat_hist_Opt,FUN="mean"),type="l",lwd=2, col="green")
lines(aggregate(mld~year,dat_hist_Dist,FUN="mean"),type="l",lwd=2,col="red")
lines(aggregate(mld~year,dat_hist_BY,FUN="mean"),type="l",lwd=2,col="yellow")
#lines(aggregate(temp~year, dat_hist_ByDist,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="purple")
lines(aggregate(mld~year,dat_hist_CA,FUN="mean"),type="l",lwd=2,col="orange")

#average Lat
par(mfrow=c(1,1))
plot(aggregate(lat~year,dat_hist_random,FUN="mean"),type="l",lwd=2, col="blue", ylim=c(33,43))
#lines(aggregate(temp~year,dat_hist_Tar,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="green")
lines(aggregate(lat~year,dat_hist_Opt,FUN="mean"),type="l",lwd=2, col="green")
lines(aggregate(lat~year,dat_hist_Dist,FUN="mean"),type="l",lwd=2,col="red")
lines(aggregate(lat~year,dat_hist_BY,FUN="mean"),type="l",lwd=2,col="yellow")
#lines(aggregate(temp~year, dat_hist_ByDist,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="purple")
lines(aggregate(lat~year,dat_hist_CA,FUN="mean"),type="l",lwd=2,col="orange")

###AVERAGE TEMP of predictions...
dat_Ran<-dat_all[dat_all$gam_Ran_Abund_nochl>1,]
dat_Opt<-dat_all[dat_all$gam_Opt_Abund_nochl>1,]
dat_Dist<-dat_all[dat_all$gam_Dist_Abund_nochl>1,]
plot(aggregate(temp~year,dat_Ran,FUN="mean"),type="l",lwd=2, col="blue", ylim=c(13, 19))
lines(aggregate(temp~year, dat_suit, FUN="mean"), type="l", lwd=2, col="gray")
lines(aggregate(temp~year, dat_all, FUN="mean"), type="l", lwd=2, col="black")
lines(aggregate(temp~year,dat_Opt,FUN="mean"),type="l",lwd=2, col="green")
lines(aggregate(temp~year,dat_Dist,FUN="mean"),type="l",lwd=2, col="red")
#for distance scenario the model predicts the species to be present at higher temperatures


##panel plot of ROMS env variables




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

