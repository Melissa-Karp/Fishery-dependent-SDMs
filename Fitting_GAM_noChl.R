### Fishery-dependent SDM (L^3) Estimation Model Code
## Function to run GAMs on data from OM -- with missing chla covariate & some w/ spatial...
## Function returns only the fitted and predicted values
#modified my M.Karp 12/15/20

##################################
#   Partial_models without Chla  #
##################################

#
##### Random sampling ####
#
    if("ran" %in% sampling){
      gam_Ran_P_nochl <- gam(pres~s(temp) + s(mld), data=dat_hist_random, family=binomial) 
      #plot(gam_Ran_P_nochl, pages=1)
      gam_Ran_N_nochl <- gam(log_abundance~s(temp)+ s(mld), data=dat_hist_random[dat_hist_random$abundance>0,], family=gaussian)
      #plot(gam_Ran_N_nochl, pages=1)
      
      dat_hist$presxR_nochl <- predict(gam_Ran_P_nochl, dat_hist, type="response")
      abundxR_nochl <- predict(gam_Ran_N_nochl, dat_hist, type="response") 
      dat_hist$gam_Ran_nochl <- dat_hist$presxR_nochl * exp(abundxR_nochl)  #predicted catch
      
      dat_fcast$presxR_nochl <- predict(gam_Ran_P_nochl, dat_fcast, type="response")
      abundxR_nochl <- predict(gam_Ran_N_nochl, dat_fcast, type="response")
      dat_fcast$gam_Ran_nochl <- dat_fcast$presxR_nochl * exp(abundxR_nochl)
    }
    
#
##### Preferential Sampling #######
#
    #Pref sampling - 0.5
    if("tar_0.5" %in% sampling) {
      gam_Tar_P_1_nochl <- gam(pres~s(temp)+ s(mld), data=dat_hist_Tar_1, family=binomial)
      #plot(gam_Tar_P_1_nochl, pages=1)
      gam_Tar_N_1_nochl <- gam(log_abundance~s(temp)+ s(mld), data=dat_hist_Tar_1[dat_hist_Tar_1$abundance>0,], family=gaussian)
      #plot(gam_Tar_N_1_nochl, pages=1)
      
      dat_hist$presxT_1_nochl <- predict(gam_Tar_P_1_nochl, dat_hist, type="response")
      abundxT_1_nochl <- predict(gam_Tar_N_1_nochl, dat_hist, type="response")
      dat_hist$gam_Tar_0.5_nochl <- dat_hist$presxT_1_nochl * exp(abundxT_1_nochl)
      
      dat_fcast$presxT_1_nochl <- predict(gam_Tar_P_1_nochl, dat_fcast, type="response")
      abundxT_1_nochl <- predict(gam_Tar_N_1_nochl, dat_fcast, type="response")
      dat_fcast$gam_Tar_0.5_nochl <- dat_fcast$presxT_1_nochl * exp(abundxT_1_nochl)
    }
    #Pref sampling -0.6
    if("tar_0.6" %in% sampling) {
      gam_Tar_P_2_nochl <- gam(pres~s(temp)+ s(mld), data=dat_hist_Tar_2, family=binomial)
      #plot(gam_Tar_P_2_nochl, pages=1)
      gam_Tar_N_2_nochl <- gam(log_abundance~s(temp)+ s(mld), data=dat_hist_Tar_2[dat_hist_Tar_2$abundance>0,], family=gaussian)
      #plot(gam_Tar_N_2_nochl, pages=1)
      
      dat_hist$presxT_2_nochl <- predict(gam_Tar_P_2_nochl, dat_hist, type="response")
      abundxT_2_nochl <- predict(gam_Tar_N_2_nochl, dat_hist, type="response")
      dat_hist$gam_Tar_0.6_nochl <- dat_hist$presxT_2_nochl * exp(abundxT_2_nochl)
      
      dat_fcast$presxT_2_nochl <- predict(gam_Tar_P_2_nochl, dat_fcast, type="response")
      abundxT_2_nochl <- predict(gam_Tar_N_2_nochl, dat_fcast, type="response")
      dat_fcast$gam_Tar_0.6_nochl <- dat_fcast$presxT_2_nochl * exp(abundxT_2_nochl)
    }
    #Pref sampling -0.7
    if("tar_0.7" %in% sampling) {
      gam_Tar_P_3_nochl <- gam(pres~s(temp)+ s(mld), data=dat_hist_Tar_3, family=binomial)
      #plot(gam_Tar_P_3_nochl, pages=1)
      gam_Tar_N_3_nochl <- gam(log_abundance~s(temp)+ s(mld), data=dat_hist_Tar_3[dat_hist_Tar_3$abundance>0,], family=gaussian)
      #plot(gam_Tar_N_3_nochl, pages=1)
      
      dat_hist$presxT_3_nochl <- predict(gam_Tar_P_3_nochl, dat_hist, type="response")
      abundxT_3_nochl <- predict(gam_Tar_N_3_nochl, dat_hist, type="response")
      dat_hist$gam_Tar_0.7_nochl <- dat_hist$presxT_3_nochl * exp(abundxT_3_nochl)
      
      dat_fcast$presxT_3_nochl <- predict(gam_Tar_P_3_nochl, dat_fcast, type="response")
      abundxT_3_nochl <- predict(gam_Tar_N_3_nochl, dat_fcast, type="response")
      dat_fcast$gam_Tar_0.7_nochl <- dat_fcast$presxT_3_nochl * exp(abundxT_3_nochl)
    }
    #Pref sampling -0.8
    if("tar_0.8" %in% sampling) {
      gam_Tar_P_4_nochl <- gam(pres~s(temp)+ s(mld), data=dat_hist_Tar_4, family=binomial)
      #plot(gam_Tar_P_4_nochl, pages=1)
      gam_Tar_N_4_nochl <- gam(log_abundance~s(temp)+ s(mld), data=dat_hist_Tar_4[dat_hist_Tar_4$abundance>0,], family=gaussian)
      #plot(gam_Tar_N_4_nochl, pages=1)
      
      dat_hist$presxT_4_nochl <- predict(gam_Tar_P_4_nochl, dat_hist, type="response")
      abundxT_4_nochl <- predict(gam_Tar_N_4_nochl, dat_hist, type="response")
      dat_hist$gam_Tar_0.8_nochl <- dat_hist$presxT_4_nochl * exp(abundxT_4_nochl)
      
      dat_fcast$presxT_4_nochl <- predict(gam_Tar_P_4_nochl, dat_fcast, type="response")
      abundxT_4_nochl <- predict(gam_Tar_N_4_nochl, dat_fcast, type="response")
      dat_fcast$gam_Tar_0.8_nochl <- dat_fcast$presxT_4_nochl * exp(abundxT_4_nochl)
    }
    #Pref sampling -0.9
    if("tar_0.9" %in% sampling) {
      gam_Tar_P_5_nochl <- gam(pres~s(temp)+ s(mld), data=dat_hist_Tar_5, family=binomial)
      #plot(gam_Tar_P_5_nochl, pages=1)
      gam_Tar_N_5_nochl <- gam(log_abundance~s(temp)+ s(mld), data=dat_hist_Tar_5[dat_hist_Tar_5$abundance>0,], family=gaussian)
      #plot(gam_Tar_N_5_nochl, pages=1)
      
      dat_hist$presxT_5_nochl <- predict(gam_Tar_P_5_nochl, dat_hist, type="response")
      abundxT_5_nochl <- predict(gam_Tar_N_5_nochl, dat_hist, type="response")
      dat_hist$gam_Tar_0.9_nochl <- dat_hist$presxT_5_nochl* exp(abundxT_5_nochl)
      
      dat_fcast$presxT_5_nochl <- predict(gam_Tar_P_5_nochl, dat_fcast, type="response")
      abundxT_5_nochl <- predict(gam_Tar_N_5_nochl, dat_fcast, type="response")
      dat_fcast$gam_Tar_0.9_nochl <- dat_fcast$presxT_5_nochl * exp(abundxT_5_nochl)
    }
    
#
##### Distance from Port Sampling #######
#
    #Dist sampling - NPO
    if("npo" %in% sampling){
      gam_Dist_P_npo_nochl <- gam(pres~s(temp) + s(mld), data=dat_hist_Dist_npo, family=binomial)
      #plot(gam_Dist_P_npo_nochl, pages=1)
      gam_Dist_N_npo_nochl <- gam(log_abundance~s(temp) + s(mld), data=dat_hist_Dist_npo[dat_hist_Dist_npo$abundance>0,], family=gaussian)
      #plot(gam_Dist_N_npo_nochl, pages=1)
      
      dat_hist$presxD_npo_nochl <- predict(gam_Dist_P_npo_nochl, dat_hist, type="response")
      abundxD_npo_nochl <- predict(gam_Dist_N_npo_nochl, dat_hist, type="response")
      dat_hist$gam_Dist_npo_nochl <- dat_hist$presxD_npo_nochl * exp(abundxD_npo_nochl)
      
      dat_fcast$presxD_npo_nochl <- predict(gam_Dist_P_npo_nochl, dat_fcast, type="response")
      abundxD_npo_nochl <- predict(gam_Dist_N_npo_nochl, dat_fcast, type="response")
      dat_fcast$gam_Dist_npo_nochl <- dat_fcast$presxD_npo_nochl * exp(abundxD_npo_nochl)
    }
    #Dist sampling - NPN
    if("npn" %in% sampling){
      gam_Dist_P_npn_nochl <- gam(pres~s(temp) + s(mld), data=dat_hist_Dist_npn, family=binomial)
      #plot(gam_Dist_P_npn_nochl, pages=1)
      gam_Dist_N_npn_nochl <- gam(log_abundance~s(temp) + s(mld), data=dat_hist_Dist_npn[dat_hist_Dist_npn$abundance>0,], family=gaussian)
      #plot(gam_Dist_N_npn_nochl, pages=1)
      
      dat_hist$presxD_npn_nochl <- predict(gam_Dist_P_npn_nochl, dat_hist, type="response")
      abundxD_npn_nochl <- predict(gam_Dist_N_npn_nochl, dat_hist, type="response")
      dat_hist$gam_Dist_npn_nochl <- dat_hist$presxD_npn_nochl * exp(abundxD_npn_nochl)
      
      dat_fcast$presxD_npn_nochl <- predict(gam_Dist_P_npn_nochl, dat_fcast, type="response")
      abundxD_npn_nochl <- predict(gam_Dist_N_npn_nochl, dat_fcast, type="response")
      dat_fcast$gam_Dist_npn_nochl <- dat_fcast$presxD_npn_nochl * exp(abundxD_npn_nochl)
    }
    #Dist sampling - MPO
    if("mpo" %in% sampling){
      gam_Dist_P_mpo_nochl <- gam(pres~s(temp) + s(mld), data=dat_hist_Dist_mpo, family=binomial)
      #plot(gam_Dist_P_mpo_nochl, pages=1)
      gam_Dist_N_mpo_nochl <- gam(log_abundance~s(temp) + s(mld), data=dat_hist_Dist_mpo[dat_hist_Dist_mpo$abundance>0,], family=gaussian)
      #plot(gam_Dist_N_mpo_nochl, pages=1)
      
      dat_hist$presxD_mpo_nochl <- predict(gam_Dist_P_mpo_nochl, dat_hist, type="response")
      abundxD_mpo_nochl <- predict(gam_Dist_N_mpo_nochl, dat_hist, type="response")
      dat_hist$gam_Dist_mpo_nochl <- dat_hist$presxD_mpo_nochl * exp(abundxD_mpo_nochl)
      
      dat_fcast$presxD_mpo_nochl <- predict(gam_Dist_P_mpo_nochl, dat_fcast, type="response")
      abundxD_mpo_nochl <- predict(gam_Dist_N_mpo_nochl, dat_fcast, type="response")
      dat_fcast$gam_Dist_mpo_nochl <- dat_fcast$presxD_mpo_nochl * exp(abundxD_mpo_nochl)
    }
    #Dist sampling - MPN
    if("mpn" %in% sampling){
      gam_Dist_P_mpn_nochl <- gam(pres~s(temp) + s(mld), data=dat_hist_Dist_mpn, family=binomial)
      #plot(gam_Dist_P_mpn_nochl, pages=1)
      gam_Dist_N_mpn_nochl <- gam(log_abundance~s(temp) + s(mld), data=dat_hist_Dist_mpn[dat_hist_Dist_mpn$abundance>0,], family=gaussian)
      #plot(gam_Dist_N_mpn_nochl, pages=1)
      
      dat_hist$presxD_mpn_nochl <- predict(gam_Dist_P_mpn_nochl, dat_hist, type="response")
      abundxD_mpn_nochl <- predict(gam_Dist_N_mpn_nochl, dat_hist, type="response")
      dat_hist$gam_Dist_mpn_nochl <- dat_hist$presxD_mpn_nochl * exp(abundxD_mpn_nochl)
      
      dat_fcast$presxD_mpn_nochl <- predict(gam_Dist_P_mpn_nochl, dat_fcast, type="response")
      abundxD_mpn_nochl <- predict(gam_Dist_N_mpn_nochl, dat_fcast, type="response")
      dat_fcast$gam_Dist_mpn_nochl <- dat_fcast$presxD_mpn_nochl * exp(abundxD_mpn_nochl)
    }
    #Dist sampling - SPO
    if("spo" %in% sampling){
      gam_Dist_P_spo_nochl <- gam(pres~s(temp) + s(mld), data=dat_hist_Dist_spo, family=binomial)
      #plot(gam_Dist_P_spo_nochl, pages=1)
      gam_Dist_N_spo_nochl <- gam(log_abundance~s(temp) + s(mld), data=dat_hist_Dist_spo[dat_hist_Dist_spo$abundance>0,], family=gaussian)
      #plot(gam_Dist_N_spo_nochl, pages=1)
      
      dat_hist$presxD_spo_nochl <- predict(gam_Dist_P_spo_nochl, dat_hist, type="response")
      abundxD_spo_nochl <- predict(gam_Dist_N_spo_nochl, dat_hist, type="response")
      dat_hist$gam_Dist_spo_nochl <- dat_hist$presxD_spo_nochl * exp(abundxD_spo_nochl)
      
      dat_fcast$presxD_spo_nochl <- predict(gam_Dist_P_spo_nochl, dat_fcast, type="response")
      abundxD_spo_nochl <- predict(gam_Dist_N_spo_nochl, dat_fcast, type="response")
      dat_fcast$gam_Dist_spo_nochl <- dat_fcast$presxD_spo_nochl * exp(abundxD_spo_nochl)
    }
    #Dist sampling - SPN
    if("spn" %in% sampling){
      gam_Dist_P_spn_nochl <- gam(pres~s(temp) + s(mld), data=dat_hist_Dist_spn, family=binomial)
      #plot(gam_Dist_P_spn_nochl, pages=1)
      gam_Dist_N_spn_nochl <- gam(log_abundance~s(temp) + s(mld), data=dat_hist_Dist_spn[dat_hist_Dist_spn$abundance>0,], family=gaussian)
      #plot(gam_Dist_N_spn_nochl, pages=1)
      
      dat_hist$presxD_spn_nochl <- predict(gam_Dist_P_spn_nochl, dat_hist, type="response")
      abundxD_spn_nochl <- predict(gam_Dist_N_spn_nochl, dat_hist, type="response")
      dat_hist$gam_Dist_spn_nochl <- dat_hist$presxD_spn_nochl * exp(abundxD_spn_nochl)
      
      dat_fcast$presxD_spn_nochl <- predict(gam_Dist_P_spn_nochl, dat_fcast, type="response")
      abundxD_spn_nochl <- predict(gam_Dist_N_spn_nochl, dat_fcast, type="response")
      dat_fcast$gam_Dist_spn_nochl <- dat_fcast$presxD_spn_nochl * exp(abundxD_spn_nochl)
    }
    #Dist sampling - ALLO
    if("allo" %in% sampling){
      gam_Dist_P_allo_nochl <- gam(pres~s(temp) + s(mld), data=dat_hist_Dist_allo, family=binomial)
      #plot(gam_Dist_P_allo_nochl, pages=1)
      gam_Dist_N_allo_nochl <- gam(log_abundance~s(temp) + s(mld), data=dat_hist_Dist_allo[dat_hist_Dist_allo$abundance>0,], family=gaussian)
      #plot(gam_Dist_N_allo_nochl, pages=1)
      
      dat_hist$presxD_allo_nochl <- predict(gam_Dist_P_allo_nochl, dat_hist, type="response")
      abundxD_allo_nochl <- predict(gam_Dist_N_allo_nochl, dat_hist, type="response")
      dat_hist$gam_Dist_allo_nochl <- dat_hist$presxD_allo_nochl * exp(abundxD_allo_nochl)
      
      dat_fcast$presxD_allo_nochl <- predict(gam_Dist_P_allo_nochl, dat_fcast, type="response")
      abundxD_allo_nochl <- predict(gam_Dist_N_allo_nochl, dat_fcast, type="response")
      dat_fcast$gam_Dist_allo_nochl <- dat_fcast$presxD_allo_nochl * exp(abundxD_allo_nochl)
    }
    #Dist sampling - Alln
    if("alln" %in% sampling){
      gam_Dist_P_alln_nochl <- gam(pres~s(temp) + s(mld), data=dat_hist_Dist_alln, family=binomial)
      #plot(gam_Dist_P_alln_nochl, pages=1)
      gam_Dist_N_alln_nochl<- gam(log_abundance~s(temp) + s(mld), data=dat_hist_Dist_alln[dat_hist_Dist_alln$abundance>0,], family=gaussian)
      #plot(gam_Dist_N_alln_nochl, pages=1)
      
      dat_hist$presxD_alln_nochl <- predict(gam_Dist_P_alln_nochl, dat_hist, type="response")
      abundxD_alln_nochl <- predict(gam_Dist_N_alln_nochl, dat_hist, type="response")
      dat_hist$gam_Dist_alln_nochl <- dat_hist$presxD_alln_nochl * exp(abundxD_alln_nochl)
      
      dat_fcast$presxD_alln_nochl <- predict(gam_Dist_P_alln, dat_fcast, type="response")
      abundxD_alln_nochl <- predict(gam_Dist_N_alln, dat_fcast, type="response")
      dat_fcast$gam_Dist_alln_nochl <- dat_fcast$presxD_alln_nochl * exp(abundxD_alln_nochl)
    }
    
#
######BY sampling
#
    if("BY" %in% sampling){
      gam_BY_P_nochl <- gam(pres~s(temp) + s(mld), data=dat_hist_BY, family=binomial)
      #plot(gam_BY_P_nochl, pages=1)
      gam_BY_N_nochl<- gam(log_abundance~s(temp) + s(mld), data=dat_hist_BY[dat_hist_BY$abundance>0,], family=gaussian)
      #plot(gam_BY_N_nochl, pages=1)
      
      dat_hist$presxB_nochl <- predict(gam_BY_P_nochl, dat_hist, type="response")
      abundxB_nochl <- predict(gam_BY_N_nochl, dat_hist, type="response")
      dat_hist$gam_BY_nochl <- dat_hist$presxB_nochl * exp(abundxB_nochl)
      
      dat_fcast$presxB_nochl <- predict(gam_BY_P_nochl, dat_fcast, type="response")
      abundxB_nochl <- predict(gam_BY_N_nochl, dat_fcast, type="response")
      dat_fcast$gam_BY_nochl <- dat_fcast$presxB_nochl * exp(abundxB_nochl)
    }
    
#
#####Closed Area Sampling
#
    #Closed Area - Small
    if("CA_sm" %in% sampling) {
      gam_CA_P_sm_nochl <- gam(pres~s(temp) + s(mld), data=dat_hist_CA_sm, family=binomial)
      #plot(gam_CA_P_sm_nochl, pages=1)
      gam_CA_N_sm_nochl <- gam(log_abundance~s(temp) + s(mld), data=dat_hist_CA_sm[dat_hist_CA_sm$abundance>0,], family=gaussian)
      #plot(gam_CA_N_sm_nochl, pages=1)
      
      dat_hist$presxCA_sm_nochl <- predict(gam_CA_P_sm_nochl, dat_hist, type="response")
      abundxCA_sm_nochl <- predict(gam_CA_N_sm_nochl, dat_hist, type="response")
      dat_hist$gam_CA_sm_nochl <- dat_hist$presxCA_sm_nochl * exp(abundxCA_sm_nochl)
      
      dat_fcast$presxCA_sm_nochl <- predict(gam_CA_P_sm_nochl, dat_fcast, type="response")
      abundxCA_sm_nochl <- predict(gam_CA_N_sm_nochl, dat_fcast, type="response")
      dat_fcast$gam_CA_sm_nochl <- dat_fcast$presxCA_sm_nochl * exp(abundxCA_sm_nochl)
    }
    #Closed Area - Medium
    if("CA_med" %in% sampling) {
      gam_CA_P_med_nochl <- gam(pres~s(temp) + s(mld), data=dat_hist_CA_med, family=binomial)
      #plot(gam_CA_P_med_nochl, pages=1)
      gam_CA_N_med_nochl <- gam(log_abundance~s(temp) + s(mld), data=dat_hist_CA_med[dat_hist_CA_med$abundance>0,], family=gaussian)
      #plot(gam_CA_N_med_nochl, pages=1)
      
      dat_hist$presxCA_med_nochl <- predict(gam_CA_P_med_nochl, dat_hist, type="response")
      abundxCA_med_nochl <- predict(gam_CA_N_med_nochl, dat_hist, type="response")
      dat_hist$gam_CA_med_nochl <- dat_hist$presxCA_med_nochl * exp(abundxCA_med_nochl)
      
      dat_fcast$presxCA_med_nochl <- predict(gam_CA_P_med_nochl, dat_fcast, type="response")
      abundxCA_med_nochl <- predict(gam_CA_N_med_nochl, dat_fcast, type="response")
      dat_fcast$gam_CA_med_nochl <- dat_fcast$presxCA_med_nochl * exp(abundxCA_med_nochl)
    }
    #Closed Area - Large
    if("CA_lar" %in% sampling) {
      gam_CA_P_lar_nochl <- gam(pres~s(temp) + s(mld), data=dat_hist_CA_lar, family=binomial)
      #plot(gam_CA_P_lar_nochl, pages=1)
      gam_CA_N_lar_nochl <- gam(log_abundance~s(temp) + s(mld), data=dat_hist_CA_lar[dat_hist_CA_lar$abundance>0,], family=gaussian)
      #plot(gam_CA_N_lar_nochl, pages=1)
      
      dat_hist$presxCA_lar_nochl <- predict(gam_CA_P_lar, dat_hist, type="response")
      abundxCA_lar_nochl <- predict(gam_CA_N_lar, dat_hist, type="response")
      dat_hist$gam_CA_lar_nochl <- dat_hist$presxCA_lar_nochl * exp(abundxCA_lar_nochl)
      
      dat_fcast$presxCA_lar_nochl <- predict(gam_CA_P_lar_nochl, dat_fcast, type="response")
      abundxCA_lar_nochl <- predict(gam_CA_N_lar_nochl, dat_fcast, type="response")
      dat_fcast$gam_CA_lar_nochl <- dat_fcast$presxCA_lar_nochl * exp(abundxCA_lar_nochl)
    }

#### PLOTS ####
#NOTE: if you are running individual models (not all of them) then will need to
# "#" the below out and turn on the relevant lines for the specific models
# being run above for it to work 
#Env Cov 1
  #Pres-Abs
  par(mfrow=c(3,3))
  plot(gam_Ran_P_nochl, select=1, main="gam_Ran, Pres-Abs", scale=0)
  plot(gam_Tar_P_1_nochl, select=1, main="gam_Tar_0.5, Pres-Abs", scale=0)
  plot(gam_Tar_P_2_nochl, select=1, main="gam_Tar_0.6, Pres-Abs", scale=0)
  plot(gam_Tar_P_3_nochl, select=1, main="gam_Tar_0.7, Pres-Abs", scale=0)
  plot(gam_Tar_P_4_nochl, select=1, main="gam_Tar_0.8, Pres-Abs", scale=0)
  plot(gam_Tar_P_5_nochl, select=1, main="gam_Tar_0.9, Pres-Abs", scale=0)
  plot(gam_Dist_P_npo_nochl, select=1, main="gam_Dist_npo, Pres-Abs", scale=0)
  plot(gam_Dist_P_npn_nochl, select=1, main="gam_Dist_npn, Pres-Abs", scale=0)
  plot(gam_Dist_P_spo_nochl, select=1, main="gam_Dist_spo, Pres-Abs", scale=0)
  plot(gam_Dist_P_mpo_nochl, select=1, main="gam_Dist_mpo, Pres-Abs", scale=0)
  plot(gam_Dist_P_mpn_nochl, select=1, main="gam_Dist_mpn, Pres-Abs", scale=0)
  plot(gam_Dist_P_allo_nochl, select=1, main="gam_Dist_allo, Pres-Abs", scale=0)
  plot(gam_Dist_P_alln_nochl, select=1, main="gam_Dist_alln, Pres-Abs", scale=0)
  plot(gam_BY_P_nochl, select=1, main="gam_BY, Pres-Abs", scale=0)
  plot(gam_CA_P_sm_nochl, select=1, main="gam_CA_sm, Pres-Abs", scale=0)
  plot(gam_CA_P_med_nochl, select=1, main="gam_CA_med, Pres-Abs", scale=0)
  plot(gam_CA_P_lar_nochl, select=1, main="gam_CA_lar, Pres-Abs", scale=0)
  
  #Abund
  par(mfrow=c(3,3))
  plot(gam_Ran_N_nochl, select=1, main="gam_Ran Nochl, Abund", scale=0)
  plot(gam_Tar_N_1_nochl, select=1, main="gam_Tar_0.5 Nochl, Abund", scale=0)
  plot(gam_Tar_N_2_nochl, select=1, main="gam_Tar_0.6 Nochl, Abund", scale=0)
  plot(gam_Tar_N_3_nochl, select=1, main="gam_Tar_0.7 Nochl, Abund", scale=0)
  plot(gam_Tar_N_4_nochl, select=1, main="gam_Tar_0.8 Nochl, Abund", scale=0)
  plot(gam_Tar_N_5_nochl, select=1, main="gam_Tar_0.9 Nochl, Abund", scale=0)
  plot(gam_Dist_N_npo_nochl, select=1, main="gam_Dist_npo Nochl, Abund", scale=0)
  plot(gam_Dist_N_npn_nochl, select=1, main="gam_Dist_npn Nochl, Abund", scale=0)
  plot(gam_Dist_N_spo_nochl, select=1, main="gam_Dist_spo Nochl, Abund", scale=0)
  plot(gam_Dist_N_mpo_nochl, select=1, main="gam_Dist_mpo Nochl, Abund", scale=0)
  plot(gam_Dist_N_mpn_nochl, select=1, main="gam_Dist_mpn Nochl, Abund", scale=0)
  plot(gam_Dist_N_allo_nochl, select=1, main="gam_Dist_allo Nochl, Abund", scale=0)
  plot(gam_Dist_N_alln_nochl, select=1, main="gam_Dist_alln Nochl, Abund", scale=0)
  plot(gam_BY_N_nochl, select=1, main="gam_BY NoChl, Abund", scale=0)
  plot(gam_CA_N_sm_nochl, select=1, main="gam_CA_sm Nochl, Abund", scale=0)
  plot(gam_CA_N_med_nochl, select=1, main="gam_CA_med Nochl, Abund", scale=0)
  plot(gam_CA_N_lar_nochl, select=1, main="gam_CA_lar Nochl, Abund", scale=0)

#Env Cov 2
  #Pres-Abs
  par(mfrow=c(3,3))
  plot(gam_Ran_P_nochl, select=2, main="gam_Ran, Pres-Abs", scale=0)
  plot(gam_Tar_P_1_nochl, select=2, main="gam_Tar_0.5, Pres-Abs", scale=0)
  plot(gam_Tar_P_2_nochl, select=2, main="gam_Tar_0.6, Pres-Abs", scale=0)
  plot(gam_Tar_P_3_nochl, select=2, main="gam_Tar_0.7, Pres-Abs", scale=0)
  plot(gam_Tar_P_4_nochl, select=2, main="gam_Tar_0.8, Pres-Abs", scale=0)
  plot(gam_Tar_P_5_nochl, select=2, main="gam_Tar_0.9, Pres-Abs", scale=0)
  plot(gam_Dist_P_npo_nochl, select=2, main="gam_Dist_npo, Pres-Abs", scale=0)
  plot(gam_Dist_P_npn_nochl, select=2, main="gam_Dist_npn, Pres-Abs", scale=0)
  plot(gam_Dist_P_spo_nochl, select=2, main="gam_Dist_spo, Pres-Abs", scale=0)
  plot(gam_Dist_P_mpo_nochl, select=2, main="gam_Dist_mpo, Pres-Abs", scale=0)
  plot(gam_Dist_P_mpn_nochl, select=2, main="gam_Dist_mpn, Pres-Abs", scale=0)
  plot(gam_Dist_P_allo_nochl, select=2, main="gam_Dist_allo, Pres-Abs", scale=0)
  plot(gam_Dist_P_alln_nochl, select=2, main="gam_Dist_alln, Pres-Abs", scale=0)
  plot(gam_BY_P_nochl, select=2, main="gam_BY, Pres-Abs", scale=0)
  plot(gam_CA_P_sm_nochl, select=2, main="gam_CA_sm, Pres-Abs", scale=0)
  plot(gam_CA_P_med_nochl, select=2, main="gam_CA_med, Pres-Abs", scale=0)
  plot(gam_CA_P_lar_nochl, select=2, main="gam_CA_lar, Pres-Abs", scale=0)
  
  #Abund
  par(mfrow=c(3,3))
  plot(gam_Ran_N_nochl, select=2, main="gam_Ran Nochl, Abund", scale=0)
  plot(gam_Tar_N_1_nochl, select=2, main="gam_Tar_0.5 Nochl, Abund", scale=0)
  plot(gam_Tar_N_2_nochl, select=2, main="gam_Tar_0.6 Nochl, Abund", scale=0)
  plot(gam_Tar_N_3_nochl, select=2, main="gam_Tar_0.7 Nochl, Abund", scale=0)
  plot(gam_Tar_N_4_nochl, select=2, main="gam_Tar_0.8 Nochl, Abund", scale=0)
  plot(gam_Tar_N_5_nochl, select=2, main="gam_Tar_0.9 Nochl, Abund", scale=0)
  plot(gam_Dist_N_npo_nochl, select=2, main="gam_Dist_npo Nochl, Abund", scale=0)
  plot(gam_Dist_N_npn_nochl, select=2, main="gam_Dist_npn Nochl, Abund", scale=0)
  plot(gam_Dist_N_spo_nochl, select=2, main="gam_Dist_spo Nochl, Abund", scale=0)
  plot(gam_Dist_N_mpo_nochl, select=2, main="gam_Dist_mpo Nochl, Abund", scale=0)
  plot(gam_Dist_N_mpn_nochl, select=2, main="gam_Dist_mpn Nochl, Abund", scale=0)
  plot(gam_Dist_N_allo_nochl, select=2, main="gam_Dist_allo Nochl, Abund", scale=0)
  plot(gam_Dist_N_alln_nochl, select=2, main="gam_Dist_alln Nochl, Abund", scale=0)
  plot(gam_BY_N_nochl, select=2, main="gam_BY NoChl, Abund", scale=0)
  plot(gam_CA_N_sm_nochl, select=2, main="gam_CA_sm Nochl, Abund", scale=0)
  plot(gam_CA_N_med_nochl, select=2, main="gam_CA_med Nochl, Abund", scale=0)
  plot(gam_CA_N_lar_nochl, select=2, main="gam_CA_lar Nochl, Abund", scale=0)
  
