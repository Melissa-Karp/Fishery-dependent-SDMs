#Function to calcualte spatial RMSE (RMSE for each cell over the time period)
#this function will create a data.frame with the difference between the true observed
#abundance at each grid cell and the predicted abundance for each model for the historic
#and forecast years

library(dplyr)
library(tidyverse)

sdm_sqErr<-function(dat_hist, dat_fcast) {

# calculate squared error for each observation
 sqerr<-function(p,o){
  ((p-o))
 }
 all_mods <-   c("gam_Ran_Abund", "gam_Tar_0.5_Abund", "gam_Tar_0.6_Abund", "gam_Tar_0.7_Abund", 
                 "gam_Tar_0.8_Abund", "gam_Tar_0.9_Abund", "gam_Dist_npo_Abund", "gam_Dist_npn_Abund",
                 "gam_Dist_mpo_Abund", "gam_Dist_mpn_Abund", "gam_Dist_spo_Abund", "gam_Dist_spn_Abund",
                 "gam_Dist_allo_Abund", "gam_Dist_alln_Abund", "gam_BY_Abund", "gam_CA_small_Abund", 
                 "gam_CA_med_Abund", "gam_CA_lar_Abund","brt_Ran","brt_Tar_0.5","brt_Tar_0.6",
                 "brt_Tar_0.7", "brt_Tar_0.8", "brt_Tar_0.9", "brt_Dist_npo", "brt_Dist_npn",
                 "brt_Dist_mpo","brt_Dist_mpn","brt_Dist_spo","brt_Dist_spn","brt_Dist_allo",
                 "brt_Dist_alln", "brt_BY", "brt_CA_small", "brt_CA_med", "brt_CA_lar",
                 "gam_Ran_Abund_nochl", "gam_Tar_0.5_Abund_nochl", "gam_Tar_0.6_Abund_nochl",
                 "gam_Tar_0.7_Abund_nochl", "gam_Tar_0.8_Abund_nochl", "gam_Tar_0.9_Abund_nochl",
                 "gam_Dist_npo_Abund_nochl", "gam_Dist_npn_Abund_nochl","gam_Dist_mpo_Abund_nochl",
                 "gam_Dist_mpn_Abund_nochl","gam_Dist_spo_Abund_nochl","gam_Dist_spn_Abund_nochl",
                 "gam_Dist_allo_Abund_nochl", "gam_Dist_alln_Abund_nochl","gam_BY_Abund_nochl",
                 "gam_CA_small_Abund_nochl", "gam_CA_med_Abund_nochl", "gam_CA_lar_Abund_nochl",
                 "brt_Ran_nochl", "brt_Tar_0.5_nochl", "brt_Tar_0.6_nochl", "brt_Tar_0.7_nochl",
                 "brt_Tar_0.8_nochl", "brt_Tar_0.9_nochl","brt_Dist_npo_nochl", "brt_Dist_npn_nochl",
                 "brt_Dist_mpo_nochl","brt_Dist_mpn_nochl","brt_Dist_spo_nochl","brt_Dist_spn_nochl",
                 "brt_Dist_allo_nochl", "brt_Dist_alln_nochl", "brt_BY_nochl", "brt_CA_small_nochl",
                 "brt_CA_med_nochl", "brt_CA_lar_nochl")  #all possible models
 
 mods <- names(dat_hist)[names(dat_hist) %in% all_mods]   #which of the possible models have been run
 
 #spatial error dataframe for historic period
 spatial_err_all_hist<- as.data.frame(matrix(0, nrow=nrow(dat_hist), ncol=length(mods)+3))
 names(spatial_err_all_hist) <- c("lat","lon","year", mods)
 spatial_err_all_hist$lat <- dat_hist$lat
 spatial_err_all_hist$lon <- dat_hist$lon
 spatial_err_all_hist$year<- dat_hist$year
 
 #spatial error dataframe for forecast period
 spatial_err_all_fcast<- as.data.frame(matrix(0, nrow=nrow(dat_fcast), ncol=length(mods)+3))
 names(spatial_err_all_fcast) <- c("lat","lon", "year", mods)
 spatial_err_all_fcast$lat <- dat_fcast$lat
 spatial_err_all_fcast$lon <- dat_fcast$lon
 spatial_err_all_fcast$year<- dat_fcast$year
 
 #Historic period Spatial error
 for (m in mods) {
   spatial_err_all_hist[,m] <- round(sqerr(dat_hist$abundance, dat_hist[,m]),3)
   }

 #Forecast period Spatial error
 for (m in mods) {
     spatial_err_all_fcast[,m] <- round(sqerr(dat_fcast$abundance, dat_fcast[,m]),3)
 }
 
 #group and summarize data tables

 return(list(spatial_err_all_fcast = spatial_err_all_fcast,
             spatial_err_all_hist = spatial_err_all_hist))
}

 # # group by cell(i.e., lat, lon) and find RMSE
 # group_by(lat, lon) %>% 
 #   summarise(rmse=sqrt(mean(sq_err))) --- maybe this code is used outside of the function
 
