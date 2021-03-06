### Fishery-dependent SDM (L^3) Estimation Model Code
## Function to run GAMs on data from OM -- with missing chla covariate & some w/ spatial...
## Function returns only the fitted and predicted values
#modified my M.Karp 1/15/21

############################################################
#   Partial_models without Chla and with Space Smoother #
############################################################

#
##### Random sampling ####
#
if("ran" %in% sampling){
  print("Fitting GAM-Ran_nochl_S")
  gam_Ran_P_nochl_S <- gam(pres~s(temp) + s(mld) + s(lat,lon), data=dat_hist_random, family=binomial) 
  #plot(gam_Ran_P_nochl_S, pages=1)
  gam_Ran_N_nochl_S <- gam(log_abundance~s(temp)+ s(mld) + s(lat,lon), data=dat_hist_random[dat_hist_random$abundance>0,], family=gaussian)
  #plot(gam_Ran_N_nochl_S, pages=1)
  
  dat_hist$presxR_nochl_S <- predict(gam_Ran_P_nochl_S, dat_hist, type="response")
  abundxR_nochl_S <- predict(gam_Ran_N_nochl_S, dat_hist, type="response") 
  dat_hist$gam_Ran_nochl_S <- dat_hist$presxR_nochl_S * exp(abundxR_nochl_S)  #predicted catch
  
  dat_fcast$presxR_nochl_S <- predict(gam_Ran_P_nochl_S, dat_fcast, type="response")
  abundxR_nochl_S <- predict(gam_Ran_N_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_Ran_nochl_S <- dat_fcast$presxR_nochl_S * exp(abundxR_nochl_S)
}

#
##### Preferential Sampling #######
#
#Pref sampling - 0.5
if("tar_0.5" %in% sampling) {
  print("Fitting GAM-Pref0.5_nochl_S")
  gam_Tar_P_1_nochl_S <- gam(pres~s(temp)+ s(mld) + s(lat,lon), data=dat_hist_Tar_1, family=binomial)
  #plot(gam_Tar_P_1_nochl_S, pages=1)
  gam_Tar_N_1_nochl_S <- gam(log_abundance~s(temp)+ s(mld) + s(lat,lon), data=dat_hist_Tar_1[dat_hist_Tar_1$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_1_nochl_S, pages=1)
  
  dat_hist$presxT_1_nochl_S <- predict(gam_Tar_P_1_nochl_S, dat_hist, type="response")
  abundxT_1_nochl_S <- predict(gam_Tar_N_1_nochl_S, dat_hist, type="response")
  dat_hist$gam_Tar_0.5_nochl_S <- dat_hist$presxT_1_nochl_S * exp(abundxT_1_nochl_S)
  
  dat_fcast$presxT_1_nochl_S <- predict(gam_Tar_P_1_nochl_S, dat_fcast, type="response")
  abundxT_1_nochl_S <- predict(gam_Tar_N_1_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.5_nochl_S <- dat_fcast$presxT_1_nochl_S * exp(abundxT_1_nochl_S)
}
#Pref sampling -0.6
if("tar_0.6" %in% sampling) {
  print("Fitting GAM-Pref0.6_nochl_S")
  gam_Tar_P_2_nochl_S <- gam(pres~s(temp)+ s(mld) + s(lat,lon), data=dat_hist_Tar_2, family=binomial)
  #plot(gam_Tar_P_2_nochl_S, pages=1)
  gam_Tar_N_2_nochl_S <- gam(log_abundance~s(temp)+ s(mld) + s(lat,lon), data=dat_hist_Tar_2[dat_hist_Tar_2$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_2_nochl_S, pages=1)
  
  dat_hist$presxT_2_nochl_S <- predict(gam_Tar_P_2_nochl_S, dat_hist, type="response")
  abundxT_2_nochl_S <- predict(gam_Tar_N_2_nochl_S, dat_hist, type="response")
  dat_hist$gam_Tar_0.6_nochl_S <- dat_hist$presxT_2_nochl_S * exp(abundxT_2_nochl_S)
  
  dat_fcast$presxT_2_nochl_S <- predict(gam_Tar_P_2_nochl_S, dat_fcast, type="response")
  abundxT_2_nochl_S <- predict(gam_Tar_N_2_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.6_nochl_S <- dat_fcast$presxT_2_nochl_S * exp(abundxT_2_nochl_S)
}
#Pref sampling -0.7
if("tar_0.7" %in% sampling) {
  print("Fitting GAM-Pref0.7_nochl_S")
  gam_Tar_P_3_nochl_S <- gam(pres~s(temp)+ s(mld) + s(lat,lon), data=dat_hist_Tar_3, family=binomial)
  #plot(gam_Tar_P_3_nochl_S, pages=1)
  gam_Tar_N_3_nochl_S <- gam(log_abundance~s(temp)+ s(mld) + s(lat,lon), data=dat_hist_Tar_3[dat_hist_Tar_3$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_3_nochl_S, pages=1)
  
  dat_hist$presxT_3_nochl_S <- predict(gam_Tar_P_3_nochl_S, dat_hist, type="response")
  abundxT_3_nochl_S <- predict(gam_Tar_N_3_nochl_S, dat_hist, type="response")
  dat_hist$gam_Tar_0.7_nochl_S <- dat_hist$presxT_3_nochl_S * exp(abundxT_3_nochl_S)
  
  dat_fcast$presxT_3_nochl_S <- predict(gam_Tar_P_3_nochl_S, dat_fcast, type="response")
  abundxT_3_nochl_S <- predict(gam_Tar_N_3_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.7_nochl_S <- dat_fcast$presxT_3_nochl_S * exp(abundxT_3_nochl_S)
}
#Pref sampling -0.8
if("tar_0.8" %in% sampling) {
  print("Fitting GAM-Pref0.8_nochl_S")
  gam_Tar_P_4_nochl_S <- gam(pres~s(temp)+ s(mld) + s(lat,lon), data=dat_hist_Tar_4, family=binomial)
  #plot(gam_Tar_P_4_nochl_S, pages=1)
  gam_Tar_N_4_nochl_S <- gam(log_abundance~s(temp)+ s(mld) + s(lat,lon), data=dat_hist_Tar_4[dat_hist_Tar_4$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_4_nochl_S, pages=1)
  
  dat_hist$presxT_4_nochl_S <- predict(gam_Tar_P_4_nochl_S, dat_hist, type="response")
  abundxT_4_nochl_S <- predict(gam_Tar_N_4_nochl_S, dat_hist, type="response")
  dat_hist$gam_Tar_0.8_nochl_S <- dat_hist$presxT_4_nochl_S * exp(abundxT_4_nochl_S)
  
  dat_fcast$presxT_4_nochl_S <- predict(gam_Tar_P_4_nochl_S, dat_fcast, type="response")
  abundxT_4_nochl_S <- predict(gam_Tar_N_4_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.8_nochl_S <- dat_fcast$presxT_4_nochl_S * exp(abundxT_4_nochl_S)
}
#Pref sampling -0.9
if("tar_0.9" %in% sampling) {
  print("Fitting GAM-Pref0.9_nochl_S")
  gam_Tar_P_5_nochl_S <- gam(pres~s(temp)+ s(mld) + s(lat,lon), data=dat_hist_Tar_5, family=binomial)
  #plot(gam_Tar_P_5_nochl_S, pages=1)
  gam_Tar_N_5_nochl_S <- gam(log_abundance~s(temp)+ s(mld) + s(lat,lon), data=dat_hist_Tar_5[dat_hist_Tar_5$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_5_nochl_S, pages=1)
  
  dat_hist$presxT_5_nochl_S <- predict(gam_Tar_P_5_nochl_S, dat_hist, type="response")
  abundxT_5_nochl_S <- predict(gam_Tar_N_5_nochl_S, dat_hist, type="response")
  dat_hist$gam_Tar_0.9_nochl_S <- dat_hist$presxT_5_nochl_S* exp(abundxT_5_nochl_S)
  
  dat_fcast$presxT_5_nochl_S <- predict(gam_Tar_P_5_nochl_S, dat_fcast, type="response")
  abundxT_5_nochl_S <- predict(gam_Tar_N_5_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.9_nochl_S <- dat_fcast$presxT_5_nochl_S * exp(abundxT_5_nochl_S)
}

#
##### Distance from Port Sampling #######
#
#Dist sampling - NPO
if("npo" %in% sampling){
  print("Fitting GAM-NPO_nochl_S")
  gam_Dist_P_npo_nochl_S <- gam(pres~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_npo, family=binomial)
  #plot(gam_Dist_P_npo_nochl_S, pages=1)
  gam_Dist_N_npo_nochl_S <- gam(log_abundance~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_npo[dat_hist_Dist_npo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_npo_nochl_S, pages=1)
  
  dat_hist$presxD_npo_nochl_S <- predict(gam_Dist_P_npo_nochl_S, dat_hist, type="response")
  abundxD_npo_nochl_S <- predict(gam_Dist_N_npo_nochl_S, dat_hist, type="response")
  dat_hist$gam_Dist_npo_nochl_S <- dat_hist$presxD_npo_nochl_S * exp(abundxD_npo_nochl_S)
  
  dat_fcast$presxD_npo_nochl_S <- predict(gam_Dist_P_npo_nochl_S, dat_fcast, type="response")
  abundxD_npo_nochl_S <- predict(gam_Dist_N_npo_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_Dist_npo_nochl_S <- dat_fcast$presxD_npo_nochl_S * exp(abundxD_npo_nochl_S)
}
#Dist sampling - NPN
if("npn" %in% sampling){
  print("Fitting GAM-NPN_nochl_S")
  gam_Dist_P_npn_nochl_S <- gam(pres~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_npn, family=binomial)
  #plot(gam_Dist_P_npn_nochl_S, pages=1)
  gam_Dist_N_npn_nochl_S <- gam(log_abundance~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_npn[dat_hist_Dist_npn$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_npn_nochl_S, pages=1)
  
  dat_hist$presxD_npn_nochl_S <- predict(gam_Dist_P_npn_nochl_S, dat_hist, type="response")
  abundxD_npn_nochl_S <- predict(gam_Dist_N_npn_nochl_S, dat_hist, type="response")
  dat_hist$gam_Dist_npn_nochl_S <- dat_hist$presxD_npn_nochl_S * exp(abundxD_npn_nochl_S)
  
  dat_fcast$presxD_npn_nochl_S <- predict(gam_Dist_P_npn_nochl_S, dat_fcast, type="response")
  abundxD_npn_nochl_S <- predict(gam_Dist_N_npn_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_Dist_npn_nochl_S <- dat_fcast$presxD_npn_nochl_S * exp(abundxD_npn_nochl_S)
}

#Dist sampling - MPO
if("mpo" %in% sampling){
  print("Fitting GAM-MPO_nochl_S")
  gam_Dist_P_mpo_nochl_S <- gam(pres~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_mpo, family=binomial)
  #plot(gam_Dist_P_mpo_nochl_S, pages=1)
  gam_Dist_N_mpo_nochl_S <- gam(log_abundance~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_mpo[dat_hist_Dist_mpo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_mpo_nochl_S, pages=1)
  
  dat_hist$presxD_mpo_nochl_S <- predict(gam_Dist_P_mpo_nochl_S, dat_hist, type="response")
  abundxD_mpo_nochl_S <- predict(gam_Dist_N_mpo_nochl_S, dat_hist, type="response")
  dat_hist$gam_Dist_mpo_nochl_S <- dat_hist$presxD_mpo_nochl_S * exp(abundxD_mpo_nochl_S)
  
  dat_fcast$presxD_mpo_nochl_S <- predict(gam_Dist_P_mpo_nochl_S, dat_fcast, type="response")
  abundxD_mpo_nochl_S <- predict(gam_Dist_N_mpo_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_Dist_mpo_nochl_S <- dat_fcast$presxD_mpo_nochl_S * exp(abundxD_mpo_nochl_S)
}
#Dist sampling - MPN
if("mpn" %in% sampling){
  print("Fitting GAM-NPN_nochl_S")
  gam_Dist_P_mpn_nochl_S <- gam(pres~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_mpn, family=binomial)
  #plot(gam_Dist_P_mpn_nochl_S, pages=1)
  gam_Dist_N_mpn_nochl_S <- gam(log_abundance~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_mpn[dat_hist_Dist_mpn$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_mpn_nochl_S, pages=1)
  
  dat_hist$presxD_mpn_nochl_S <- predict(gam_Dist_P_mpn_nochl_S, dat_hist, type="response")
  abundxD_mpn_nochl_S <- predict(gam_Dist_N_mpn_nochl_S, dat_hist, type="response")
  dat_hist$gam_Dist_mpn_nochl_S <- dat_hist$presxD_mpn_nochl_S * exp(abundxD_mpn_nochl_S)
  
  dat_fcast$presxD_mpn_nochl_S <- predict(gam_Dist_P_mpn_nochl_S, dat_fcast, type="response")
  abundxD_mpn_nochl_S <- predict(gam_Dist_N_mpn_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_Dist_mpn_nochl_S <- dat_fcast$presxD_mpn_nochl_S * exp(abundxD_mpn_nochl_S)
}
#Dist sampling - SPO
if("spo" %in% sampling){
  print("Fitting GAM-SPO_nochl_S")
  gam_Dist_P_spo_nochl_S <- gam(pres~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_spo, family=binomial)
  #plot(gam_Dist_P_spo_nochl_S, pages=1)
  gam_Dist_N_spo_nochl_S <- gam(log_abundance~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_spo[dat_hist_Dist_spo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_spo_nochl_S, pages=1)
  
  dat_hist$presxD_spo_nochl_S <- predict(gam_Dist_P_spo_nochl_S, dat_hist, type="response")
  abundxD_spo_nochl_S <- predict(gam_Dist_N_spo_nochl_S, dat_hist, type="response")
  dat_hist$gam_Dist_spo_nochl_S <- dat_hist$presxD_spo_nochl_S * exp(abundxD_spo_nochl_S)
  
  dat_fcast$presxD_spo_nochl_S <- predict(gam_Dist_P_spo_nochl_S, dat_fcast, type="response")
  abundxD_spo_nochl_S <- predict(gam_Dist_N_spo_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_Dist_spo_nochl_S <- dat_fcast$presxD_spo_nochl_S * exp(abundxD_spo_nochl_S)
}
#Dist sampling - SPN
if("spn" %in% sampling){
  print("Fitting GAM-SPN_nochl_S")
  gam_Dist_P_spn_nochl_S <- gam(pres~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_spn, family=binomial)
  #plot(gam_Dist_P_spn_nochl_S, pages=1)
  gam_Dist_N_spn_nochl_S <- gam(log_abundance~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_spn[dat_hist_Dist_spn$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_spn_nochl_S, pages=1)
  
  dat_hist$presxD_spn_nochl_S <- predict(gam_Dist_P_spn_nochl_S, dat_hist, type="response")
  abundxD_spn_nochl_S <- predict(gam_Dist_N_spn_nochl_S, dat_hist, type="response")
  dat_hist$gam_Dist_spn_nochl_S <- dat_hist$presxD_spn_nochl_S * exp(abundxD_spn_nochl_S)
  
  dat_fcast$presxD_spn_nochl_S <- predict(gam_Dist_P_spn_nochl_S, dat_fcast, type="response")
  abundxD_spn_nochl_S <- predict(gam_Dist_N_spn_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_Dist_spn_nochl_S <- dat_fcast$presxD_spn_nochl_S * exp(abundxD_spn_nochl_S)
}
#Dist sampling - ALLO
if("allo" %in% sampling){
  print("Fitting GAM-ALLO_nochl_S")
  gam_Dist_P_allo_nochl_S <- gam(pres~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_allo, family=binomial)
  #plot(gam_Dist_P_allo_nochl_S, pages=1)
  gam_Dist_N_allo_nochl_S <- gam(log_abundance~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_allo[dat_hist_Dist_allo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_allo_nochl_S, pages=1)
  
  dat_hist$presxD_allo_nochl_S <- predict(gam_Dist_P_allo_nochl_S, dat_hist, type="response")
  abundxD_allo_nochl_S <- predict(gam_Dist_N_allo_nochl_S, dat_hist, type="response")
  dat_hist$gam_Dist_allo_nochl_S <- dat_hist$presxD_allo_nochl_S * exp(abundxD_allo_nochl_S)
  
  dat_fcast$presxD_allo_nochl_S <- predict(gam_Dist_P_allo_nochl_S, dat_fcast, type="response")
  abundxD_allo_nochl_S <- predict(gam_Dist_N_allo_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_Dist_allo_nochl_S <- dat_fcast$presxD_allo_nochl_S * exp(abundxD_allo_nochl_S)
}
#Dist sampling - Alln
if("alln" %in% sampling){
  print("Fitting GAM-ALLN_nochl_S")
  gam_Dist_P_alln_nochl_S <- gam(pres~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_alln, family=binomial)
  #plot(gam_Dist_P_alln_nochte, pages=1)
  gam_Dist_N_alln_nochl_S <- gam(log_abundance~s(temp) + s(mld) + s(lat,lon), data=dat_hist_Dist_alln[dat_hist_Dist_alln$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_alln_nochl_S, pages=1)
  
  dat_hist$presxD_alln_nochl_S <- predict(gam_Dist_P_alln_nochl_S, dat_hist, type="response")
  abundxD_alln_nochl_S <- predict(gam_Dist_N_alln_nochl_S, dat_hist, type="response")
  dat_hist$gam_Dist_alln_nochl_S <- dat_hist$presxD_alln_nochl_S * exp(abundxD_alln_nochl_S)
  
  dat_fcast$presxD_alln_nochl_S <- predict(gam_Dist_P_alln_nochl_S, dat_fcast, type="response")
  abundxD_alln_nochl_S <- predict(gam_Dist_N_alln_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_Dist_alln_nochl_S <- dat_fcast$presxD_alln_nochl_S * exp(abundxD_alln_nochl_S)
}

#
######BY sampling
#
if("BY" %in% sampling){
  print("Fitting GAM-BY_nochl_S")
  gam_BY_P_nochl_S <- gam(pres~s(temp) + s(mld)+ s(lat,lon), data=dat_hist_BY, family=binomial)
  #plot(gam_BY_P_nochl_S, pages=1)
  gam_BY_N_nochl_S<- gam(log_abundance~s(temp) + s(mld)+ s(lat,lon), data=dat_hist_BY[dat_hist_BY$abundance>0,], family=gaussian)
  #plot(gam_BY_N_nochl_S, pages=1)
  
  dat_hist$presxB_nochl_S <- predict(gam_BY_P_nochl_S, dat_hist, type="response")
  abundxB_nochl_S <- predict(gam_BY_N_nochl_S, dat_hist, type="response")
  dat_hist$gam_BY_nochl_S <- dat_hist$presxB_nochl_S * exp(abundxB_nochl_S)
  
  dat_fcast$presxB_nochl_S <- predict(gam_BY_P_nochl_S, dat_fcast, type="response")
  abundxB_nochl_S <- predict(gam_BY_N_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_BY_nochl_S <- dat_fcast$presxB_nochl_S * exp(abundxB_nochl_S)
}

#
#####Closed Area Sampling
#
#Closed Area - Small
if("CA_sm" %in% sampling) {
  print("Fitting GAM-CASM_nochl_S")
  gam_CA_P_sm_nochl_S <- gam(pres~s(temp) + s(mld) + s(lat,lon), data=dat_hist_CA_sm, family=binomial)
  #plot(gam_CA_P_sm_nochl_S, pages=1)
  gam_CA_N_sm_nochl_S <- gam(log_abundance~s(temp) + s(mld) + s(lat,lon), data=dat_hist_CA_sm[dat_hist_CA_sm$abundance>0,], family=gaussian)
  #plot(gam_CA_N_sm_nochl_S, pages=1)
  
  dat_hist$presxCA_sm_nochl_S <- predict(gam_CA_P_sm_nochl_S, dat_hist, type="response")
  abundxCA_sm_nochl_S <- predict(gam_CA_N_sm_nochl_S, dat_hist, type="response")
  dat_hist$gam_CA_sm_nochl_S <- dat_hist$presxCA_sm_nochl_S * exp(abundxCA_sm_nochl_S)
  
  dat_fcast$presxCA_sm_nochl_S <- predict(gam_CA_P_sm_nochl_S, dat_fcast, type="response")
  abundxCA_sm_nochl_S <- predict(gam_CA_N_sm_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_CA_sm_nochl_S <- dat_fcast$presxCA_sm_nochl_S * exp(abundxCA_sm_nochl_S)
}
#Closed Area - Medium
if("CA_med" %in% sampling) {
  print("Fitting GAM-CAMED_nochl_S")
  gam_CA_P_med_nochl_S <- gam(pres~s(temp) + s(mld) + s(lat,lon), data=dat_hist_CA_med, family=binomial)
  #plot(gam_CA_P_med_nochl_S, pages=1)
  gam_CA_N_med_nochl_S <- gam(log_abundance~s(temp) + s(mld) + s(lat,lon), data=dat_hist_CA_med[dat_hist_CA_med$abundance>0,], family=gaussian)
  #plot(gam_CA_N_med_nochl_S, pages=1)
  
  dat_hist$presxCA_med_nochl_S <- predict(gam_CA_P_med_nochl_S, dat_hist, type="response")
  abundxCA_med_nochl_S <- predict(gam_CA_N_med_nochl_S, dat_hist, type="response")
  dat_hist$gam_CA_med_nochl_S <- dat_hist$presxCA_med_nochl_S * exp(abundxCA_med_nochl_S)
  
  dat_fcast$presxCA_med_nochl_S <- predict(gam_CA_P_med_nochl_S, dat_fcast, type="response")
  abundxCA_med_nochl_S <- predict(gam_CA_N_med_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_CA_med_nochl_S <- dat_fcast$presxCA_med_nochl_S * exp(abundxCA_med_nochl_S)
}
#Closed Area - Large
if("CA_lar" %in% sampling) {
  print("Fitting GAM-CALAR_nochl_S")
  gam_CA_P_lar_nochl_S <- gam(pres~ s(temp) + s(mld) + s(lat,lon), data=dat_hist_CA_lar, family=binomial)
  #plot(gam_CA_P_lar_nochl_S, pages=1)
  gam_CA_N_lar_nochl_S <- gam(log_abundance~ s(temp) + s(mld) +  s(lat,lon), data=dat_hist_CA_lar[dat_hist_CA_lar$abundance>0,], family=gaussian)
  #plot(gam_CA_N_lar_nochl_S, pages=1)
  
  dat_hist$presxCA_lar_nochl_S <- predict(gam_CA_P_lar_nochl_S, dat_hist, type="response")
  abundxCA_lar_nochl_S <- predict(gam_CA_N_lar_nochl_S, dat_hist, type="response")
  dat_hist$gam_CA_lar_nochl_S <- dat_hist$presxCA_lar_nochl_S * exp(abundxCA_lar_nochl_S)
  
  dat_fcast$presxCA_lar_nochl_S <- predict(gam_CA_P_lar_nochl_S, dat_fcast, type="response")
  abundxCA_lar_nochl_S <- predict(gam_CA_N_lar_nochl_S, dat_fcast, type="response")
  dat_fcast$gam_CA_lar_nochl_S<- dat_fcast$presxCA_lar_nochl_S * exp(abundxCA_lar_nochl_S)
}
