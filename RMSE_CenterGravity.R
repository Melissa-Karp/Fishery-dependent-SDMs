### Location^3 Estimation Model Code
## Function to measure latitudinal Centre Of Gravity for each SDM EM

## Code by Steph Brodie and James Smith - April 2020


sdm_cog <- function(dat_hist, dat_fcast) {
  #'dat_hist' are the 'historical' observations from the OM, and predicted abundance values from fitted models as additional columns
  #'dat_fcast' are the 'future' observations from the OM,  and predicted abundance values from fitted models as additional columns
  #'mods' are the types of models to report: c("gam", "glm", "brt", "mlp")
  
  # Perf. metrics are:
  #       - Plots of weighted annual mean latitude (COG)
  #       - RMSE of annual COG values
  
  
  RMSE = function(p, o){
    sqrt(mean((p - o)^2))
  }
  
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
                  "brt_CA_med_nochl_S", "brt_CA_lar_nochl_S") 
  
  mods <- names(dat_hist)[names(dat_hist) %in% all_mods]   #which of the possible models have been run
  
  rmse_cog <- data.frame(model=mods, rmse_hist=0, rmse_fcast=0)
  
  years_hist <- unique(dat_hist$year)
  years_fcast <- unique(dat_fcast$year)
  
  cog_lat_hist <- as.data.frame(matrix(0,nrow=length(years_hist),ncol=2+length(mods)))
  cog_lat_fcast <- as.data.frame(matrix(0,nrow=length(years_fcast),ncol=2+length(mods)))
  names(cog_lat_hist) <- c("year", "truth", mods)
  names(cog_lat_fcast) <- c("year", "truth", mods)
  cog_lat_hist$year <- years_hist
  cog_lat_fcast$year <- years_fcast
  
  cog_lat_fcast_d <- as.data.frame(matrix(0, nrow=length(years_fcast), ncol=length(mods)+1))
  names(cog_lat_fcast_d) <- c("year", mods)
  cog_lat_fcast_d$year <- unique(dat_fcast$year)
  
  
  # Calculate COG - hist period
  for (y in cog_lat_hist$year) {
    cog_lat_truth <- weighted.mean(dat_hist[dat_hist$year==y, "lat"],w=dat_hist[dat_hist$year==y, "abundance"])
    cog_lat_hist[cog_lat_hist$year==y, "truth"] <- cog_lat_truth
    for (m in mods) {
      cog_lat_mod <- weighted.mean(dat_hist[dat_hist$year==y, "lat"],w=dat_hist[dat_hist$year==y, m])
      cog_lat_hist[cog_lat_hist$year==y, m] <- cog_lat_mod
    }
  }
  
  # Calculate COG - fcast period
  for (y in cog_lat_fcast$year) {
    cog_lat_truth <- weighted.mean(dat_fcast[dat_fcast$year==y, "lat"],w=dat_fcast[dat_fcast$year==y, "abundance"])
    cog_lat_fcast[cog_lat_fcast$year==y, "truth"] <- cog_lat_truth
    for (m in mods) {
      cog_lat_mod <- weighted.mean(dat_fcast[dat_fcast$year==y, "lat"],w=dat_fcast[dat_fcast$year==y, m])
      cog_lat_fcast[cog_lat_fcast$year==y, m] <- cog_lat_mod
      
      #calculate difference between observed and predicted COG
      cog_lat_fcast_d[cog_lat_fcast_d$year == y, m] <- cog_lat_truth - cog_lat_mod
    }
  }
  
  # Plot COG - hist period
  par(mfrow=c(3,3))
  for (m in mods) {
    rmse_mod <- RMSE(cog_lat_hist$truth, cog_lat_hist[,m])
    rmse_cog$rmse_hist[rmse_cog$model==m] <- rmse_mod
    
    plot(cog_lat_hist$year,cog_lat_hist$truth, type='l', ylab="Latitude", xlim=c(min(years_hist),max(years_hist)),
         main=paste0("COG-lat, Hist, ", m, ", rmse=", round(rmse_mod,3)), lwd=2)
    lines(cog_lat_hist$year,cog_lat_hist[,m], col="blue", lwd=2)
  }
  
  # Plot COG - fcast period
  par(mfrow=c(3,3))
  for (m in mods) {
    rmse_mod <- RMSE(cog_lat_fcast$truth, cog_lat_fcast[,m])
    rmse_cog$rmse_fcast[rmse_cog$model==m] <- rmse_mod
    
    plot(cog_lat_fcast$year,cog_lat_fcast$truth, type='l', ylab="Latitude", xlim=c(min(years_fcast),max(years_fcast)),
         main=paste0("COG-lat, Fcast, ", m, ", rmse=", round(rmse_mod,3)), lwd=2)
    lines(cog_lat_fcast$year,cog_lat_fcast[,m], col="blue", lwd=2)
  }
  
  #plot obs-pred difference in COG during fcast period
  par(mfrow=c(3,3), mar=c(2.5,4,2.5,1))
  for (m in mods) {
    plot(cog_lat_fcast_d$year, cog_lat_fcast_d[,m], type="l", lwd=2, col="red",
         ylab="obs-pred cog", xlab="", main=paste(m),
         ylim=c(min(cog_lat_fcast_d[,m]),max(cog_lat_fcast_d[,m])))
    abline(h=0, lty=2, col="blue")
  }
  
  return(list(rmse_cog=rmse_cog,
              cog_lat_hist=cog_lat_hist,
              cog_lat_fcast=cog_lat_fcast,
              cog_lat_fcast_d=cog_lat_fcast_d))
}
