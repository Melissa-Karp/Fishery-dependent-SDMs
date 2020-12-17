### Location^3 Estimation Model Code
## Create maps of prediction species distribution - observed and fitted SDMs

## Code by James Smith (april 2020) and modified by Melissa Karp December 2020 

SDM_maps <- function(dat_hist, dat_fcast, YEARS){
  
  #combine data from histori and future into one dataframe
  dat_fcast$log_abundance <- log(dat_fcast$abundance)
  dat_all<-rbind(dat_hist, dat_fcast)
  
  #Maps from GAMs - Full models
  par(mfrow=c(3,3), mar=c(2,3,3,4))
  for (yy in YEARS) {
    
    year_x <- yy
    dat_x <- dat_all[dat_all$year==year_x,]
    
    r_obs <- rasterFromXYZ(dat_x[,c("lon","lat","abundance")])
    plot(r_obs, asp=1, main=paste0("Observed, year ", year_x))
    
    r_Ran <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Ran")])
    plot(r_Ran, asp=1, main=paste0("Random Full, year ", year_x))
    
    r_Tar_1 <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Tar_0.5")])
    plot(r_Tar_1, asp=1, main=paste0("Tar 0.5 Full, year ", year_x))
    
    r_Tar_2 <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Tar_0.6")])
    plot(r_Tar_2, asp=1, main=paste0("Tar 0.6 Full, year ", year_x))
    
    r_Tar_3 <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Tar_0.7")])
    plot(r_Tar_3, asp=1, main=paste0("Tar 0.7 Full, year ", year_x))
    
    r_Tar_4 <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Tar_0.8")])
    plot(r_Tar_4, asp=1, main=paste0("Tar 0.8 Full, year ", year_x))
    
    r_Tar_5 <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Tar_0.9")])
    plot(r_Tar_5, asp=1, main=paste0("Tar 0.9 Full, year ", year_x))
    
    r_npo <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_npo")])
    plot(r_npo, asp=1, main=paste0("Dist NPO Full, year ", year_x))
    
    r_npn <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_npn")])
    plot(r_npn, asp=1, main=paste0("Dist NPN Full, year ", year_x))
    
    r_mpo <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_mpo")])
    plot(r_mpo, asp=1, main=paste0("Dist MPO Full, year ", year_x))
    
    r_mpn <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_mpn")])
    plot(r_mpn, asp=1, main=paste0("Dist MPN Full, year ", year_x))
    
    r_spo <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_spo")])
    plot(r_spo, asp=1, main=paste0("Dist SPO Full, year ", year_x))
    
    r_spn <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_spn")])
    plot(r_spn, asp=1, main=paste0("Dist SPN Full, year ", year_x))
    
    r_allo <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_allo")])
    plot(r_allo, asp=1, main=paste0("Dist MPO Full, year ", year_x))
    
    r_alln <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_alln")])
    plot(r_alln, asp=1, main=paste0("Dist MPN Full, year ", year_x))
    
    r_BY <- rasterFromXYZ(dat_x[,c("lon","lat","gam_BY")])
    plot(r_BY, asp=1, main=paste0("Bycatch Full, year ", year_x))
    
    r_sm<- rasterFromXYZ(dat_x[,c("lon","lat","gam_CA_sm")])
    plot(r_sm, asp=1, main=paste0("CA Sm Full, year ", year_x))
    
    r_med <- rasterFromXYZ(dat_x[,c("lon","lat","gam_CA_med")])
    plot(r_med, asp=1, main=paste0("CA Med Full, year ", year_x))
    
    r_lar <- rasterFromXYZ(dat_x[,c("lon","lat","gam_CA_lar")])
    plot(r_lar, asp=1, main=paste0("CA Lar Full, year ", year_x))
  }
  
  #Maps from GAMs - Missing Covariate models
  par(mfrow=c(3,3), mar=c(2,3,3,4))
  for (yy in YEARS) {
    
    year_x <- yy
    dat_x <- dat_all[dat_all$year==year_x,]
    
    r_obs <- rasterFromXYZ(dat_x[,c("lon","lat","abundance")])
    plot(r_obs, asp=1, main=paste0("Observed, year ", year_x))
    
    r_Ran <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Ran_nochl")])
    plot(r_Ran, asp=1, main=paste0("Random Partial, year ", year_x))
    
    r_Tar_1 <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Tar_0.5_nochl")])
    plot(r_Tar_1, asp=1, main=paste0("Tar 0.5 Partial, year ", year_x))
    
    r_Tar_2 <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Tar_0.6_nochl")])
    plot(r_Tar_2, asp=1, main=paste0("Tar 0.6 Partial, year ", year_x))
    
    r_Tar_3 <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Tar_0.7_nochl")])
    plot(r_Tar_3, asp=1, main=paste0("Tar 0.7 Partial, year ", year_x))
    
    r_Tar_4 <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Tar_0.8_nochl")])
    plot(r_Tar_4, asp=1, main=paste0("Tar 0.8 Partial, year ", year_x))
    
    r_Tar_5 <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Tar_0.9_nochl")])
    plot(r_Tar_5, asp=1, main=paste0("Tar 0.9 Partial, year ", year_x))
    
    r_npo <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_npo_nochl")])
    plot(r_npo, asp=1, main=paste0("Dist NPO Partial, year ", year_x))
    
    r_npn <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_npn_nochl")])
    plot(r_npn, asp=1, main=paste0("Dist NPN Partial, year ", year_x))
    
    r_mpo <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_mpo_nochl")])
    plot(r_mpo, asp=1, main=paste0("Dist MPO Partial, year ", year_x))
    
    r_mpn <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_mpn_nochl")])
    plot(r_mpn, asp=1, main=paste0("Dist MPN Partial, year ", year_x))
    
    r_spo <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_spo_nochl")])
    plot(r_spo, asp=1, main=paste0("Dist SPO Partial, year ", year_x))
    
    r_spn <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_spn_nochl")])
    plot(r_spn, asp=1, main=paste0("Dist SPN Partial, year ", year_x))
    
    r_allo <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_allo_nochl")])
    plot(r_allo, asp=1, main=paste0("Dist MPO Partial, year ", year_x))
    
    r_alln <- rasterFromXYZ(dat_x[,c("lon","lat","gam_Dist_alln_nochl")])
    plot(r_alln, asp=1, main=paste0("Dist MPN Partial, year ", year_x))
    
    r_BY <- rasterFromXYZ(dat_x[,c("lon","lat","gam_BY_nochl")])
    plot(r_BY, asp=1, main=paste0("Bycatch Partial, year ", year_x))
    
    r_sm<- rasterFromXYZ(dat_x[,c("lon","lat","gam_CA_sm_nochl")])
    plot(r_sm, asp=1, main=paste0("CA Sm Partial, year ", year_x))
    
    r_med <- rasterFromXYZ(dat_x[,c("lon","lat","gam_CA_med_nochl")])
    plot(r_med, asp=1, main=paste0("CA Med Partial, year ", year_x))
    
    r_lar <- rasterFromXYZ(dat_x[,c("lon","lat","gam_CA_lar_nochl")])
    plot(r_lar, asp=1, main=paste0("CA Lar Partial, year ", year_x))
  }
}