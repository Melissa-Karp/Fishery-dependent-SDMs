
SDM_maps <- function(dat_hist, dat_fcast, YEARS){
  
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

mods <- names(dat_hist)[names(dat_hist) %in% all_mods]   #which of the possible models have been run


#combine data from historic and future into one dataframe
  # dat_fcast$log_abundance <- log(dat_fcast$abundance)
  # dat_hist<-dat_hist[-c(337)]
  # dat_fcast<-dat_fcast[-c(336)]
  
  dat_all<-rbind(dat_hist, dat_fcast)
  
 #SDM plots - all models
   par(mfrow=c(6,3), mar=c(2,2,2,1))
  for (yy in YEARS){
    for (m in mods) {
    
    year_x <- yy
    dat_x <- dat_all[dat_all$year==year_x,]
    
    r_obs <- rasterFromXYZ(dat_x[,c("lon","lat",m)])
    plot(r_obs, asp=0.8,main=paste0(m," - ", year_x))
    }
  }
}
test<-SDM_maps(dat_hist=dat_hist, dat_fcast=dat_fcast, YEARS=c(2000, 2040,2100))  
