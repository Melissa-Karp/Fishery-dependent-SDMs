sdm_spErr <- function(dat_hist, dat_fcast) {
  
  # calculate the relative error for each observation: RE = (pred-obs)/obs
  sperr <- function(p,o) {
    (p-o)
  }
  
  #list all possible models
  all_mods <-   c("gam_Ran", "gam_Tar_0.5", "gam_Tar_0.6", "gam_Tar_0.7","gam_Tar_0.8", 
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
                  "gam_CA_med_nochl_S", "gam_CA_lar_nochl_S",
                  
                  "brt_Ran","brt_Tar_0.5","brt_Tar_0.6", "brt_Tar_0.7", "brt_Tar_0.8", "brt_Tar_0.9", 
                  "brt_Dist_npo","brt_Dist_npn","brt_Dist_mpo","brt_Dist_mpn","brt_Dist_spo","brt_Dist_spn",
                  "brt_Dist_allo","brt_Dist_alln", "brt_BY", "brt_CA_sm", "brt_CA_med", "brt_CA_lar",
                  
                  "brt_Ran_nochl", "brt_Tar_0.5_nochl", "brt_Tar_0.6_nochl", "brt_Tar_0.7_nochl", 
                  "brt_Tar_0.8_nochl", "brt_Tar_0.9_nochl", "brt_Dist_npo_nochl", "brt_Dist_npn_nochl",
                  "brt_Dist_mpo_nochl","brt_Dist_mpn_nochl", "brt_Dist_spo_nochl","brt_Dist_spn_nochl",
                  "brt_Dist_allo_nochl", "brt_Dist_alln_nochl","brt_BY_nochl", "brt_CA_sm_nochl",
                  "brt_CA_med_nochl", "brt_CA_lar_nochl",
                  
                  "brt_Ran_te","brt_Tar_0.5_te","brt_Tar_0.6_te", "brt_Tar_0.7_te", "brt_Tar_0.8_te", "brt_Tar_0.9_te", 
                  "brt_Dist_npo_te","brt_Dist_npn_te","brt_Dist_mpo_te","brt_Dist_mpn_te","brt_Dist_spo_te",
                  "brt_Dist_spn_te","brt_Dist_allo_te","brt_Dist_alln_te", "brt_BY_te", "brt_CA_sm_te",
                  "brt_CA_med_te", "brt_CA_lar_te",
                  
                  "brt_Ran_nochl_te", "brt_Tar_0.5_nochl_te", "brt_Tar_0.6_nochl_te", "brt_Tar_0.7_nochl_te", 
                  "brt_Tar_0.8_nochl_te", "brt_Tar_0.9_nochl_te", "brt_Dist_npo_nochl_te", "brt_Dist_npn_nochl_te",
                  "brt_Dist_mpo_nochl_te","brt_Dist_mpn_nochl_te", "brt_Dist_spo_nochl_te","brt_Dist_spn_nochl_te",
                  "brt_Dist_allo_nochl_te", "brt_Dist_alln_nochl_te","brt_BY_nochl_te", "brt_CA_sm_nochl_te",
                  "brt_CA_med_nochl_te", "brt_CA_lar_nochl_te",
                  
                  "brt_Ran_S","brt_Tar_0.5_S","brt_Tar_0.6_S", "brt_Tar_0.7_S", "brt_Tar_0.8_S", "brt_Tar_0.9_S", 
                  "brt_Dist_npo_S","brt_Dist_npn_S","brt_Dist_mpo_S","brt_Dist_mpn_S","brt_Dist_spo_S",
                  "brt_Dist_spn_S","brt_Dist_allo_S","brt_Dist_alln_S", "brt_BY_S", "brt_CA_sm_S",
                  "brt_CA_med_S", "brt_CA_lar_S",
                  
                  "brt_Ran_nochl_S", "brt_Tar_0.5_nochl_S", "brt_Tar_0.6_nochl_S", "brt_Tar_0.7_nochl_S", 
                  "brt_Tar_0.8_nochl_S", "brt_Tar_0.9_nochl_S", "brt_Dist_npo_nochl_S", "brt_Dist_npn_nochl_S",
                  "brt_Dist_mpo_nochl_S","brt_Dist_mpn_nochl_S", "brt_Dist_spo_nochl_S","brt_Dist_spn_nochl_S",
                  "brt_Dist_allo_nochl_S", "brt_Dist_alln_nochl_S","brt_BY_nochl_S", "brt_CA_sm_nochl_S",
                  "brt_CA_med_nochl_S", "brt_CA_lar_nochl_S")   #all possible models
  
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
    spatial_err_all_hist[,m] <- round(sperr(dat_hist[,m], dat_hist$abundance), 3)
  }
  
  #Forecast period Spatial error
  for (m in mods) {
    spatial_err_all_fcast[,m] <- round(sperr(dat_fcast[,m],dat_fcast$abundance), 3)
  }
  
  ## Plot maps of error
  SEH <- data.frame(spatial_err_all_hist)
  SEF <- data.frame(spatial_err_all_fcast)
  
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


