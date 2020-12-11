### Location^3 Estimation Model Code
## Function to take observed and predicted values from SDM EMs and measure their performance

## modifed from Code by James Smith - April 2020 for use in Fishery-Dependent data
#biases impacts on SDM performance project [Sept 2020]


sdm_rmse <- function(dat_hist, dat_fcast) {
  #'dat_hist' are the 'historical' observations from the OM, and predicted abundance values from fitted models as additional columns
  #'dat_fcast' are the 'future' observations from the OM,  and predicted abundance values from fitted models as additional columns
  
  # Perf. metrics are:
  #       - RMSE and Cor between entire subset of observed and predicted abundances, historical and future
  #       - Plots and RMSE of total annual abundance, historical and future
  # Maximum of 72 models (18 scenarios x 2 model algorithms x 2 covariate situations)
  
  
  RMSE = function(p, o){
    sqrt(mean((p - o)^2))
  }
  
  RRMSE = function(p, o){
    sqrt(mean((p - o)^2)) / mean(o)
  }
  
  all_mods <- c("gam_Ran_Abund", "gam_Tar_0.5_Abund", "gam_Tar_0.6_Abund", "gam_Tar_0.7_Abund", 
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
                "brt_CA_med_nochl", "brt_CA_lar_nochl")
  mods <- names(dat_hist)[names(dat_hist) %in% all_mods]   #which of the possible models have been run
  
  rmse_abund_all <- data.frame(model=mods, rmse_hist=0, rmse_fcast=0)
  rmse_abund_annual <- data.frame(model=mods, rmse_hist=0, rmse_fcast=0)
  
  rmse_abund_all_y <- as.data.frame(matrix(0, nrow=length(unique(dat_fcast$year)), ncol=length(mods)+1))
  names(rmse_abund_all_y) <- c("year", mods)
  rmse_abund_all_y$year <- unique(dat_fcast$year)
  
  rmse_abund_allh_y <- as.data.frame(matrix(0, nrow=length(unique(dat_hist$year)), ncol=length(mods)+1))
  names(rmse_abund_allh_y) <- c("year", mods)
  rmse_abund_allh_y$year <- unique(dat_hist$year)
  
  # Calculate RMSE, and plot total abundance trajectories
  
  #Hist
  par(mfrow=c(3,3), mar=c(2.5,4,2.5,1))
  
  for (m in mods) {
    rmse_abund_all$rmse_hist[rmse_abund_all$model==m]  <- round(RMSE(dat_hist$abundance, dat_hist[,m]),3)
    rmse_abund_all$rmse_fcast[rmse_abund_all$model==m]  <- round(RMSE(dat_fcast$abundance, dat_fcast[,m]),3)
    
    # hist
    ann_obs <- aggregate(abundance~year, dat=dat_hist, FUN="sum")
    ann_mod <- aggregate(formula(paste(m,"~ year")), dat=dat_hist, FUN="sum")
    rmse_abund_annual$rmse_hist[rmse_abund_all$model==m] <- RMSE(ann_obs$abundance, ann_mod[,m])
    
    plot(ann_obs, type="l", lwd=2, xlab="",
         main=paste0("Hist, ", m, ", rmse=",round(RMSE(dat_hist$abundance, dat_hist[,m]),2),
                     ", cor=", round(cor(dat_hist$abundance, dat_hist[,m]),2)))
    lines(ann_mod, lwd=2, col="red")
    
    for (yy in rmse_abund_allh_y$year) {  #rmse of all observations per year
      dat_hist_yy <- dat_hist[dat_hist$year == yy,]
      rmse_abund_allh_y[rmse_abund_allh_y$year == yy,m] <- round(RMSE(dat_hist_yy$abundance, dat_hist_yy[,m]),3)
    }
  }
  
  #Fcast
  par(mfrow=c(3,3), mar=c(2.5,4,2.5,1))
  
  for (m in mods) {
    ann_obs <- aggregate(abundance~year, dat=dat_fcast, FUN="sum")
    ann_mod <- aggregate(formula(paste(m,"~ year")), dat=dat_fcast, FUN="sum")
    rmse_abund_annual$rmse_fcast[rmse_abund_all$model==m] <- RMSE(ann_obs$abundance, ann_mod[,m])
    
    plot(ann_obs, type="l", lwd=2, xlab="",
         main=paste0("Fcast, ", m, ", rmse=",round(RMSE(dat_fcast$abundance, dat_fcast[,m]),2),
                     ", cor=", round(cor(dat_fcast$abundance, dat_fcast[,m]),2)))
    lines(ann_mod, lwd=2, col="red")
    
    for (yy in rmse_abund_all_y$year) {  #rmse of all 400 observations per year
      dat_fcast_yy <- dat_fcast[dat_fcast$year == yy,]
      rmse_abund_all_y[rmse_abund_all_y$year == yy,m] <- round(RMSE(dat_fcast_yy$abundance, dat_fcast_yy[,m]),3)
    }
  }
  
  par(mfrow=c(3,3), mar=c(2.5,4,2.5,1))
  
  for (m in mods) {  #plot degregdation of RMSE over time (if any)
    plot(rmse_abund_all_y$year, rmse_abund_all_y[,m], type="l", lwd=2, col="red",
         ylab="rmse", xlab="", main=paste(m), ylim=c(min(rmse_abund_all_y[,m]),max(rmse_abund_all_y[,m])))
    abline(h=mean(rmse_abund_allh_y[,m]), lty=2, col="blue")  #mean annual rmse of historical period
  }
  
  
  return(list(rmse_abund_all = rmse_abund_all,
              rmse_abund_annual = rmse_abund_annual,
              rmse_abund_all_year = rmse_abund_all_y,
              rmse_abund_allh_y = rmse_abund_allh_y))
  
}
