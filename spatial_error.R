## Function to calcualte spatial error (error for each cell averaged over the time period)
## this function will create a data.frame with the difference between the true observed
## abundance at each grid cell and the predicted abundance for each model for the historic
## and forecast years

## Create 4 plots per model (2 hist, 2 fcast)

library(dplyr)
library(tidyverse)
library(gridExtra)
library(raster)

sdm_sqErr <- function(dat_hist, dat_fcast) {

# calculate error for each observation: pred-obs
 sqerr <- function(p,o) {
  ((p-o))
 }
 
 #list all possible models
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
 
 #find which of the possible models have been run
 mods <- names(dat_hist)[names(dat_hist) %in% all_mods]
 
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
   spatial_err_all_hist[,m] <- round(sqerr(dat_hist[,m]), dat_hist$abundance, 3)
   }

 #Forecast period Spatial error
 for (m in mods) {
     spatial_err_all_fcast[,m] <- round(sqerr(dat_fcast[,m]),dat_fcast$abundance, 3)
 }
 
 ## Plot maps of error
 SEH <- data.frame(SRMSE$spatial_err_all_hist)
 SEF <- data.frame(SRMSE$spatial_err_all_fcast)
 
 for (m in mods) {  #loop through each model (4 plots per model)
   
   print(paste0("Calculating ", m))
   
   # group by cell(i.e., lat, lon) and find RMSE for each grid cell
   SRMSE_hist <- SEH %>%  #hist
     group_by(lat, lon) %>% 
     summarise(meanErr=mean( .data[[m]] ), .groups='drop')  #this '.data' syntax is for standard evaluation using loop index 'm'
   
   SRMSE_fcast <- SEF %>%  #fcast
     group_by(lat, lon) %>% 
     summarise(meanErr=mean( .data[[m]] ), .groups='drop')
   
   # also do an aggregated version (coarser cell size)
   rast <- rasterFromXYZ(SRMSE_hist[,c("lon","lat","meanErr")])
   rast_agg <- aggregate(rast, fact=4, fun=mean)  #'fact' determines cell size ('4' means new cell width will be 4-times wider)
   vals <- as.vector(rast_agg)
   SRMSE_hist_agg <- as.data.frame(xyFromCell(rast_agg, cell=1:length(vals)))
   names(SRMSE_hist_agg) <- c("lon", "lat")
   SRMSE_hist_agg$meanErr <- vals
   
   rast <- rasterFromXYZ(SRMSE_fcast[,c("lon","lat","meanErr")])
   rast_agg <- aggregate(rast, fact=4, fun=mean)  #'fact' determines cell size ('4' means new cell width will be 4-times wider)
   vals <- as.vector(rast_agg)
   SRMSE_fcast_agg <- as.data.frame(xyFromCell(rast_agg, cell=1:length(vals)))
   names(SRMSE_fcast_agg) <- c("lon", "lat")
   SRMSE_fcast_agg$meanErr <- vals
   
   #Hist Raw resolution (0.1 degree)
   SRMSE_hist_map <- ggplot(data = SRMSE_hist, aes(x=lon,y=lat))+
     geom_tile(aes(fill=meanErr)) +
     ggtitle(paste0("Spatial Error (pred - obs); Hist - ",m)) +
     theme_classic() +  labs(y="", x="") +
     theme(legend.position="right",legend.title = element_blank())+
     theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
     scale_fill_gradient2() +
     annotation_map(map_data("world"), colour = "black", fill="grey50")+
     coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48))
   
   #Coarser resolution (currently 0.4 degree)
   SRMSE_hist_map_agg <- ggplot(data = SRMSE_hist_agg, aes(x=lon,y=lat))+
     geom_tile(aes(fill=meanErr)) +
     ggtitle(paste0("Spatial Error (pred - obs); Hist - ",m)) +
     theme_classic() +  labs(y="", x="") +
     theme(legend.position="right",legend.title = element_blank())+
     theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
     scale_fill_gradient2() +
     annotation_map(map_data("world"), colour = "black", fill="grey50")+
     coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48))
   
   #Fcast Raw resolution (0.1 degree)
   SRMSE_fcast_map <- ggplot(data = SRMSE_fcast, aes(x=lon,y=lat))+
     geom_tile(aes(fill=meanErr)) +
     ggtitle(paste0("Spatial Error (pred - obs); Fcast - ",m)) +
     theme_classic() +  labs(y="", x="") +
     theme(legend.position="right",legend.title = element_blank())+
     theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
     scale_fill_gradient2() +
     annotation_map(map_data("world"), colour = "black", fill="grey50")+
     coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48))
   
   #Fcast resolution (currently 0.4 degree)
   SRMSE_fcast_map_agg <- ggplot(data = SRMSE_fcast_agg, aes(x=lon,y=lat))+
     geom_tile(aes(fill=meanErr)) +
     ggtitle(paste0("Spatial Error (pred - obs); Fcast - ",m)) +
     theme_classic() +  labs(y="", x="") +
     theme(legend.position="right",legend.title = element_blank())+
     theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
     scale_fill_gradient2() +
     annotation_map(map_data("world"), colour = "black", fill="grey50")+
     coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48))
   
   #plot each model
   grid.arrange(SRMSE_hist_map, SRMSE_hist_map_agg, ncol=2)
   grid.arrange(SRMSE_fcast_map, SRMSE_fcast_map_agg, ncol=2)
   
 }
 
 ## Return the raw data
 return(list(spatial_err_all_fcast = spatial_err_all_fcast,
             spatial_err_all_hist = spatial_err_all_hist))
}


