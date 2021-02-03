
setwd("~/DisMAP project/Location, Location, Location/Location Workshop")

#################################
##       Load Libraries        ##
#################################

library(mgcv)
library(ggplot2)
library(viridis)
library(raster)
library(gbm)
library(dismo)

#################################
##       Data Simulation       ##
#################################

source("SimulatedWorld_ROMS_FishDep_Final.R") #load ROMS simulation function
dir <- "~/DisMAP project/Location, Location, Location/Location Workshop/ROMS/hadley" #directory where ROMS data is stored (on dropbox, email steph for access)

dat <- SimulateWorld_ROMS_FishDepFun_Final(dir=dir, nsamples = 100) #takes a few mins
names(dat)[names(dat) == 'sst'] <- 'temp' #matching roms names. Quick temporary fix.

############################
#    Prepare Data Sets     #
############################

#Create dataframe with historical/forecast data
dat_hist_1 <- dat[dat$year<=2010,]
dat_hist<-dat_hist_1[dat_hist_1$year>1984,] ## so dat_hist used to fit the model is 1985-2010
dat_fcast <- dat[dat$year>2010,] #forecast using 2011-2100

dat_hist$log_abundance <- log(dat_hist$abundance)

dat_hist_random<-dat_hist[dat_hist$random_sampled>0,]
dat_hist_Tar_1<-dat_hist[dat_hist$pref_sampled_1>0,]
dat_hist_Tar_2<-dat_hist[dat_hist$pref_sampled_2>0,]
dat_hist_Tar_3<-dat_hist[dat_hist$pref_sampled_3>0,]
dat_hist_Tar_4<-dat_hist[dat_hist$pref_sampled_4>0,]
dat_hist_Tar_5<-dat_hist[dat_hist$pref_sampled_5>0,]
dat_hist_Dist_npo<-dat_hist[dat_hist$dist_sampled_npo>0,]
dat_hist_Dist_npn<-dat_hist[dat_hist$dist_sampled_npn>0,]
dat_hist_Dist_mpo<-dat_hist[dat_hist$dist_sampled_mpo>0,]
dat_hist_Dist_mpn<-dat_hist[dat_hist$dist_sampled_mpn>0,]
dat_hist_Dist_spo<-dat_hist[dat_hist$dist_sampled_spo>0,]
dat_hist_Dist_spn<-dat_hist[dat_hist$dist_sampled_spn>0,]
dat_hist_Dist_allo<-dat_hist[dat_hist$dist_sampled_allo>0,]
dat_hist_Dist_alln<-dat_hist[dat_hist$dist_sampled_alln>0,]
dat_hist_BY<-dat_hist[dat_hist$BY_sampled>0,]
dat_hist_CA_sm<-dat_hist[dat_hist$Closed_sampled_1>0,]
dat_hist_CA_med<-dat_hist[dat_hist$Closed_sampled_2>0,]
dat_hist_CA_lar<-dat_hist[dat_hist$Closed_sampled_3>0,]

# dat_fcast_random<-dat_fcast[dat_fcast$random_sampled>0,]
# dat_fcast_Tar_1<-dat_fcast[dat_fcast$pref_sampled_1>0,]
# dat_fcast_Tar_2<-dat_fcast[dat_fcast$pref_sampled_2>0,]
# dat_fcast_Tar_3<-dat_fcast[dat_fcast$pref_sampled_3>0,]
# dat_fcast_Tar_4<-dat_fcast[dat_fcast$pref_sampled_4>0,]
# dat_fcast_Tar_5<-dat_fcast[dat_fcast$pref_sampled_5>0,]
# dat_fcast_Dist_npo<-dat_fcast[dat_fcast$dist_sampled_npo>0,]
# dat_fcast_Dist_npn<-dat_fcast[dat_fcast$dist_sampled_npn>0,]
# dat_fcast_Dist_mpo<-dat_fcast[dat_fcast$dist_sampled_mpo>0,]
# dat_fcast_Dist_mpn<-dat_fcast[dat_fcast$dist_sampled_mpn>0,]
# dat_fcast_Dist_spo<-dat_fcast[dat_fcast$dist_sampled_spo>0,]
# dat_fcast_Dist_spn<-dat_fcast[dat_fcast$dist_sampled_spn>0,]
# dat_fcast_Dist_allo<-dat_fcast[dat_fcast$dist_sampled_allo>0,]
# dat_fcast_Dist_alln<-dat_fcast[dat_fcast$dist_sampled_alln>0,]
# dat_fcast_BY<-dat_fcast[dat_fcast$BY_sampled>0,]
# dat_fcast_CA_sm<-dat_fcast[dat_fcast$Closed_sampled_1>0,]
# dat_fcast_CA_med<-dat_fcast[dat_fcast$Closed_sampled_2>0,]
# dat_fcast_CA_lar<-dat_fcast[dat_hist$Closed_sampled_3>0,]


###########################
#   MODEL FITTING         #
###########################
#sampling scenario options: ran, tar_0.5, tar_0.6, tar_0.7, tar_0.8, tar_0.9, npo, npn, mpo,
# mpn, spo, spn, allo, alln, BY, CA_sm, CA_med, CA_lar
#total of 18 different sampling scenarios/rules X 2 model algorithms x 6 covariate configurations
#total of 216 possible models


#NOTE: if you are running individual models (not all of them) then will need to
# go into each fitting code and "#" the "PlOTS" section out and turn on the relevant 
#lines for the specific models to plot just those models being run for it to work 
  
 ### GAMS #####
  #### GAMS - Full Model #####
    sampling <- c("ran", "tar_0.5", "tar_0.6", "tar_0.7", "tar_0.8", "tar_0.9",
                  "npo", "npn", "mpo", "mpn", "spo", "spn", "allo", "alln", "BY", 
                  "CA_sm", "CA_med", "CA_lar")
    source("Fitting_GAMs.R") 
    
  #### GAMS - Missing Covariate Model #####
    sampling <- c("ran", "tar_0.5", "tar_0.6", "tar_0.7", "tar_0.8", "tar_0.9",
                  "npo", "npn", "mpo", "mpn", "spo", "spn", "allo", "alln", "BY", 
                  "CA_sm", "CA_med", "CA_lar")
    source("Fitting_GAM_noChl.R") 
    
  #### GAMS - Full Model time-space tensor #####
    sampling <- c("ran", "tar_0.5", "tar_0.6", "tar_0.7", "tar_0.8", "tar_0.9",
                  "npo", "npn", "mpo", "mpn", "spo", "spn", "allo", "alln", "BY", 
                  "CA_sm", "CA_med", "CA_lar")
    source("GAM_SpaceTime_Config3.R") 

  #### GAMS - Partial Model time-space tensor #####
    sampling <- c("ran", "tar_0.5", "tar_0.6", "tar_0.7", "tar_0.8", "tar_0.9",
                  "npo", "npn", "mpo", "mpn", "spo", "spn", "allo", "alln", "BY", 
                  "CA_sm", "CA_med", "CA_lar")
    source("GAM_SpaceTime_nochl_Config3.R") 

  #### GAMS - Full Model Space Smoother #####
  sampling <- c("ran", "tar_0.5", "tar_0.6", "tar_0.7", "tar_0.8", "tar_0.9",
                "npo", "npn", "mpo", "mpn", "spo", "spn", "allo", "alln", "BY", 
                "CA_sm", "CA_med", "CA_lar")
  source("GAM_Space_Smoother.R") 

  #### GAMS - Partial Model time-space tensor #####
  sampling <- c("ran", "tar_0.5", "tar_0.6", "tar_0.7", "tar_0.8", "tar_0.9",
                "npo", "npn", "mpo", "mpn", "spo", "spn", "allo", "alln", "BY", 
                "CA_sm", "CA_med", "CA_lar")
  source("GAM_nochl_Space_Smoother.R") 

 ###BRTS ####  
  #### Boosted Regression Trees (BRTs) - Full Model #####
    sampling <- c("ran", "tar_0.5", "tar_0.6", "tar_0.7", "tar_0.8", "tar_0.9",
                  "npo", "npn", "mpo", "mpn", "spo", "spn", "allo", "alln", "BY", 
                  "CA_sm", "CA_med", "CA_lar")
    source("Fitting_BRTs.R") 

  #### Boosted Regression Trees (BRTs) - Missing Covariate Model #####
    sampling <- c("ran", "tar_0.5", "tar_0.6", "tar_0.7", "tar_0.8", "tar_0.9",
                  "npo", "npn", "mpo", "mpn", "spo", "spn", "allo", "alln", "BY", 
                  "CA_sm", "CA_med", "CA_lar")
    source("Fitting_BRT_noChl.R") 

  #### Boosted Regression Trees (BRTs) - Full Covariate Model + Spacetime #####
  sampling <- c("ran", "tar_0.5", "tar_0.6", "tar_0.7", "tar_0.8", "tar_0.9",
                "npo", "npn", "mpo", "mpn", "spo", "spn", "allo", "alln", "BY", 
                "CA_sm", "CA_med", "CA_lar")
  source("BRT_spacetime.R") 
  
  #### Boosted Regression Trees (BRTs) - Missing Covariate Model + Spacetime#####
  sampling <- c("ran", "tar_0.5", "tar_0.6", "tar_0.7", "tar_0.8", "tar_0.9",
                "npo", "npn", "mpo", "mpn", "spo", "spn", "allo", "alln", "BY", 
                "CA_sm", "CA_med", "CA_lar")
  source("BRT_SpaceTime_nochl.R") 

  #### Boosted Regression Trees (BRTs) - Full Covariate Model + SPACE #####
  sampling <- c("ran", "tar_0.5", "tar_0.6", "tar_0.7", "tar_0.8", "tar_0.9",
                "npo", "npn", "mpo", "mpn", "spo", "spn", "allo", "alln", "BY", 
                "CA_sm", "CA_med", "CA_lar")
  source("BRT_Space.R") 
  
  #### Boosted Regression Trees (BRTs) - Missing Covariate Model + SPACE #####
  sampling <- c("ran", "tar_0.5", "tar_0.6", "tar_0.7", "tar_0.8", "tar_0.9",
                "npo", "npn", "mpo", "mpn", "spo", "spn", "allo", "alln", "BY", 
                "CA_sm", "CA_med", "CA_lar")
  source("BRT_Space_nochl.R") 


  #### save the Rdata to load later ####
  saveRDS(dat_hist, "dat_hist_results_full_1_16_21.rds")  # * 'full' [full models] or 'temp' [temp-only models]
  saveRDS(dat_fcast, "dat_fcast_results_full_1_16_21.rds")
  
  all_mods <- c("gam_Ran", "gam_Tar_0.5", "gam_Tar_0.6", "gam_Tar_0.7","gam_Tar_0.8", 
                 "gam_Tar_0.9", "gam_Dist_npo", "gam_Dist_npn","gam_Dist_mpo", "gam_Dist_mpn",
                 "gam_Dist_spo", "gam_Dist_spn", "gam_Dist_allo", "gam_Dist_alln", 
                 "gam_BY", "gam_CA_sm", "gam_CA_med", "gam_CA_lar",
                 
                 "gam_Ran_nochl", "gam_Tar_0.5_nochl", "gam_Tar_0.6_nochl", "gam_Tar_0.7_nochl",
                 "gam_Tar_0.8_nochl", "gam_Tar_0.9_nochl", "gam_Dist_npo_nochl", "gam_Dist_npn_nochl",
                 "gam_Dist_mpo_nochl","gam_Dist_mpn_nochl","gam_Dist_spo_nochl","gam_Dist_spn_nochl",
                 "gam_Dist_allo_nochl", "gam_Dist_alln_nochl","gam_BY_nochl","gam_CA_sm_nochl", 
                 "gam_CA_med_nochl", "gam_CA_lar_nochl",
                 
                 "gam_Ran_te","gam_Tar_0.5_te", "gam_Tar_0.6_te", "gam_Tar_0.7_te","gam_Tar_0.8_te", 
                 "gam_Tar_0.9_te", "gam_Dist_npo_te", "gam_Dist_npn_te","gam_Dist_mpo_te", 
                 "gam_Dist_mpn_te","gam_Dist_spo_te", "gam_Dist_spn_te", "gam_Dist_allo_te", 
                 "gam_Dist_alln_te", "gam_BY_te", "gam_CA_sm_te", "gam_CA_med_te", "gam_CA_lar_te",
                 
                 "gam_Ran_nochl_te","gam_Tar_0.5_nochl_te", "gam_Tar_0.6_nochl_te", "gam_Tar_0.7_nochl_te",
                 "gam_Tar_0.8_nochl_te", "gam_Tar_0.9_nochl_te", "gam_Dist_npo_nochl_te", 
                 "gam_Dist_npn_nochl_te","gam_Dist_mpo_nochl_te","gam_Dist_mpn_nochl_te",
                 "gam_Dist_spo_nochl_te","gam_Dist_spn_nochl_te","gam_Dist_allo_nochl_te",
                 "gam_Dist_alln_nochl_te","gam_BY_nochl_te","gam_CA_sm_nochl_te", 
                 "gam_CA_med_nochl_te", "gam_CA_lar_nochl_te", 
                 
                 "gam_Ran_S","gam_Tar_0.5_S", "gam_Tar_0.6_S", "gam_Tar_0.7_S","gam_Tar_0.8_S", 
                 "gam_Tar_0.9_S", "gam_Dist_npo_S", "gam_Dist_npn_S","gam_Dist_mpo_S", 
                 "gam_Dist_mpn_S","gam_Dist_spo_S", "gam_Dist_spn_S", "gam_Dist_allo_S", 
                 "gam_Dist_alln_S", "gam_BY_S", "gam_CA_sm_S", "gam_CA_med_S", "gam_CA_lar_S",
                 
                 "gam_Ran_nochl_S","gam_Tar_0.5_nochl_S", "gam_Tar_0.6_nochl_S", "gam_Tar_0.7_nochl_S",
                 "gam_Tar_0.8_nochl_S", "gam_Tar_0.9_nochl_S", "gam_Dist_npo_nochl_S", 
                 "gam_Dist_npn_nochl_S","gam_Dist_mpo_nochl_S","gam_Dist_mpn_nochl_S",
                 "gam_Dist_spo_nochl_S","gam_Dist_spn_nochl_S","gam_Dist_allo_nochl_S",
                 "gam_Dist_alln_nochl_S","gam_BY_nochl_S","gam_CA_sm_nochl_S", 
                 "gam_CA_med_nochl_S", "gam_CA_lar_nochl_S")
  
  save(list = ls(pattern = paste0(all_mods, collapse="|")),
       file = "saved_models_full_1_16_21.RData")  # * 'full' or 'temp'; this will save all models of these names, even when outdated, if they're in the environment


###################################
##       Model Performance       ##
###################################
  # * load these, for full models
  dat_hist <- readRDS("~/DisMAP project/Location, Location, Location/Location Workshop/dat_hist_results_full_1_16_21.rds")
  dat_fcast <- readRDS("~/DisMAP project/Location, Location, Location/Location Workshop/dat_fcast_results_full_1_16_21.rds")
  load("~/DisMAP project/Location, Location, Location/Location Workshop/saved_models_full_1_16_21.RData")  #full models
  
#### This section runs RMSE for abundance, COG, and correlation, and spatial RMSE####
  dat_fcast_early <- dat_fcast[dat_fcast$year >= 2020 & dat_fcast$year < 2060, ]
  #dat_fcast <- dat_fcast[-c(123)]
  
  ##RMSE - Abundance
  source("PerformanceMeasures_RSME.R") 
  rmse_results<-sdm_rmse(dat_hist=dat_hist, dat_fcast=dat_fcast)
  rmse_results
  rmseall<-data.frame(rmse_results$rmse_abund_all)
  rmseall_hist<-order(rmseall$rmse_hist)
  #^ Hist: ordered best -> worst for historical
  rmse_results$rmse_abund_all[order(rmse_results$rmse_abund_all$rmse_fcast),]
  #^ Fcast: ordered best -> worst for forecast
  
  #write.csv(rmseall, 'RMSE_all.csv',row.names = FALSE)
  # rmse_Gam<-read.csv("~/DisMAP project/Location, Location, Location/Location Workshop/rmse_Gam.csv")
  # rmse_Gam_nochl<-read.csv("~/DisMAP project/Location, Location, Location/Location Workshop/rmse_Gam_nochl.csv")
  # rmse_Gam_te<-read.csv("~/DisMAP project/Location, Location, Location/Location Workshop/rmse_Gam_te.csv")
  # rmse_Gam_nochl_te<-read.csv("~/DisMAP project/Location, Location, Location/Location Workshop/rmse_Gam_nochl_te.csv")
  # 
  # par(mfrow=c(2,2))
  #  barplot(rmse_Gam$rmse_hist~rmse_Gam$model, main="Avg RMSE_hist GAM", ylab="RMSE", ylim=c(0,20), cex.names=0.7, las=2)
  #  barplot(rmse_Gam_nochl$rmse_hist~rmse_Gam_nochl$model, main="Avg RMSE_hist GAM nochl", ylab="RMSE", xlab="", ylim=c(0,20), cex.names=0.7, las=2)
  #  barplot(rmse_Gam_te$rmse_hist~rmse_Gam_te$model, main="Avg RMSE_hist GAM te", ylab="RMSE", xlab="", ylim=c(0,20), cex.names=0.7, las=2)
  #  barplot(rmse_Gam_nochl_te$rmse_hist~rmse_Gam_nochl_te$model, main="Avg RMSE_hist GAM nochl te", ylab="RMSE", xlab="", ylim=c(0,20), cex.names=0.7, las=2)
  #  
  ##Correlation - Abundance
  source("Performance_COR.R") 
  Cor_sdm<-sdm_cor(dat_hist=dat_hist, dat_fcast=dat_fcast)
  
  ##Center of Gravity 
  source("RMSE_CenterGravity.R") 
  RMSE_cog<-sdm_cog(dat_hist=dat_hist, dat_fcast=dat_fcast)
  RMSE_cog
  rmsecog_y<-data.frame(RMSE_cog$cog_lat_fcast_y)

  ##Species Distribution Maps ---- NOT sure this is the best way to map these
  source("SDM_plot_maps.R") 
  SDMmaps<-SDM_maps(dat_hist=dat_hist, dat_fcast=dat_fcast, YEARS=c("2000", "2040", "2100"))
 
  #spatial error
  source("spatial_error.R") 
  SRMSE<-sdm_spErr(dat_hist=dat_hist, dat_fcast=dat_fcast)
  
  
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
plot(aggregate(temp~year, dat_suit, FUN="mean"), type="l",lwd=2,ylab="avg temp of samples", ylim=c(11,19), main="Average temp sampled")
lines(aggregate(temp~year,dat_hist_random,FUN="mean"),type="l",lwd=2,ylab="avg temp of samples", col="blue")
lines(aggregate(temp~year,dat_hist_Tar_1,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="green")
lines(aggregate(temp~year,dat_hist_Tar_2,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="dark green")
lines(aggregate(temp~year,dat_hist_Tar_3,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="dark green")
lines(aggregate(temp~year,dat_hist_Tar_4,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="dark green")
lines(aggregate(temp~year,dat_hist_Tar_5,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="dark green")
lines(aggregate(temp~year,dat_hist_Dist_npo,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="red")
lines(aggregate(temp~year,dat_hist_Dist_npn,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="pink")
lines(aggregate(temp~year,dat_hist_Dist_mpo,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="red")
lines(aggregate(temp~year,dat_hist_Dist_mpn,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="red")
lines(aggregate(temp~year,dat_hist_Dist_spo,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="red")
lines(aggregate(temp~year,dat_hist_Dist_spn,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="red")
lines(aggregate(temp~year,dat_hist_Dist_allo,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="red")
lines(aggregate(temp~year,dat_hist_Dist_alln,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="red")
lines(aggregate(temp~year,dat_hist_BY,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="yellow")
lines(aggregate(temp~year,dat_hist_CA_sm,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="orange")
lines(aggregate(temp~year,dat_hist_CA_med,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="orange")
lines(aggregate(temp~year,dat_hist_CA_lar,FUN="mean"),type="l",lwd=2,ylab="Biomass",col="orange")

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

