### Fishery-dependent SDM (L^3) Estimation Model Code
## Function to run GAMs on data from OM
## Function returns only the fitSd and predicSd values
#modified my M.Karp 1/15/20


###################################################
#   Full Models- with Chla and space Smoother  #
###################################################

##### Random sampling ####
#
if("ran" %in% sampling){
  print("Fitting GAM-Ran-S")
  gam_Ran_PS <- gam(pres~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_random, family=binomial) #noS - chl-surface is log of chl-surface
  #plot(gam_Ran_PS, pages=1)
  gam_Ran_NS <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_random[dat_hist_random$abundance>0,], family=gaussian)
  #plot(gam_Ran_NS, pages=1)
  
  dat_hist$presxRS <- predict(gam_Ran_PS, dat_hist, type="response")
  abundxRS <- predict(gam_Ran_NS, dat_hist, type="response") 
  dat_hist$gam_Ran_S <- dat_hist$presxRS * exp(abundxRS)  #predicted catch
  
  dat_fcast$presxRS <- predict(gam_Ran_PS, dat_fcast, type="response")
  abundxRS <- predict(gam_Ran_NS, dat_fcast, type="response")
  dat_fcast$gam_Ran_S <- dat_fcast$presxRS * exp(abundxRS)
}

#
##### Preferential Sampling #######
#
#Pref sampling - 0.5
if("tar_0.5" %in% sampling) {
  print("Fitting GAM-Pref0.5_S")
  gam_Tar_P_1S<- gam(pres~s(temp)+ s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Tar_1, family=binomial)
  #plot(gam_Tar_P_1S, pages=1)
  gam_Tar_N_1S <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface)+ s(lat,lon), data=dat_hist_Tar_1[dat_hist_Tar_1$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_1S, pages=1)
  
  dat_hist$presxT_1_S <- predict(gam_Tar_P_1S, dat_hist, type="response")
  abundxT_1S <- predict(gam_Tar_N_1S, dat_hist, type="response")
  dat_hist$gam_Tar_0.5_S <- dat_hist$presxT_1_S * exp(abundxT_1S)
  
  dat_fcast$presxT_1_S <- predict(gam_Tar_P_1S, dat_fcast, type="response")
  abundxT_1S <- predict(gam_Tar_N_1S, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.5_S <- dat_fcast$presxT_1_S * exp(abundxT_1S)
}
#Pref sampling -0.6
if("tar_0.6" %in% sampling) {
  print("Fitting GAM-Pref0.6_S")
  gam_Tar_P_2S <- gam(pres~s(temp)+ s(mld) + s(chl_surface)+ s(lat,lon), data=dat_hist_Tar_2, family=binomial)
  #plot(gam_Tar_P_2S, pages=1)
  gam_Tar_N_2S <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface)+ s(lat,lon), data=dat_hist_Tar_2[dat_hist_Tar_2$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_2S, pages=1)
  
  dat_hist$presxT_2S <- predict(gam_Tar_P_2S, dat_hist, type="response")
  abundxT_2S <- predict(gam_Tar_N_2S, dat_hist, type="response")
  dat_hist$gam_Tar_0.6_S <- dat_hist$presxT_2S * exp(abundxT_2S)
  
  dat_fcast$presxT_2S <- predict(gam_Tar_P_2S, dat_fcast, type="response")
  abundxT_2S <- predict(gam_Tar_N_2S, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.6_S <- dat_fcast$presxT_2S * exp(abundxT_2S)
}
#Pref sampling -0.7
if("tar_0.7" %in% sampling) {
  print("Fitting GAM-Pref0.7_S")
  gam_Tar_P_3S <- gam(pres~s(temp)+ s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Tar_3, family=binomial)
  #plot(gam_Tar_P_3S, pages=1)
  gam_Tar_N_3S <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface)+ s(lat,lon), data=dat_hist_Tar_3[dat_hist_Tar_3$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_3S, pages=1)
  
  dat_hist$presxT_3S <- predict(gam_Tar_P_3S, dat_hist, type="response")
  abundxT_3S <- predict(gam_Tar_N_3S, dat_hist, type="response")
  dat_hist$gam_Tar_0.7_S <- dat_hist$presxT_3S * exp(abundxT_3S)
  
  dat_fcast$presxT_3S <- predict(gam_Tar_P_3S, dat_fcast, type="response")
  abundxT_3S <- predict(gam_Tar_N_3S, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.7_S <- dat_fcast$presxT_3S * exp(abundxT_3S)
}
#Pref sampling -0.8
if("tar_0.8" %in% sampling) {
  print("Fitting GAM-Pref0.8_S")
  gam_Tar_P_4S <- gam(pres~s(temp)+ s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Tar_4, family=binomial)
  #plot(gam_Tar_P_4S, pages=1)
  gam_Tar_N_4S <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Tar_4[dat_hist_Tar_4$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_4S, pages=1)
  
  dat_hist$presxT_4S <- predict(gam_Tar_P_4S, dat_hist, type="response")
  abundxT_4S <- predict(gam_Tar_N_4S, dat_hist, type="response")
  dat_hist$gam_Tar_0.8_S <- dat_hist$presxT_4S * exp(abundxT_4S)
  
  dat_fcast$presxT_4S <- predict(gam_Tar_P_4S, dat_fcast, type="response")
  abundxT_4S <- predict(gam_Tar_N_4S, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.8_S <- dat_fcast$presxT_4S * exp(abundxT_4S)
}
#Pref sampling -0.9
if("tar_0.9" %in% sampling) {
  print("Fitting GAM-Pref0.9_S")
  gam_Tar_P_5S <- gam(pres~s(temp)+ s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Tar_5, family=binomial)
  #plot(gam_Tar_P_5S, pages=1)
  gam_Tar_N_5S <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Tar_5[dat_hist_Tar_5$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_5S, pages=1)
  
  dat_hist$presxT_5S <- predict(gam_Tar_P_5S, dat_hist, type="response")
  abundxT_5S <- predict(gam_Tar_N_5S, dat_hist, type="response")
  dat_hist$gam_Tar_0.9_S <- dat_hist$presxT_5S * exp(abundxT_5S)
  
  dat_fcast$presxT_5S <- predict(gam_Tar_P_5S, dat_fcast, type="response")
  abundxT_5S <- predict(gam_Tar_N_5S, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.9_S <- dat_fcast$presxT_5S * exp(abundxT_5S)
}

#
##### Distance from Port Sampling #######
#
#Dist sampling - NPO
if("npo" %in% sampling){
  print("Fitting GAM-NPO_S")
  gam_Dist_P_npoS <- gam(pres~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Dist_npo, family=binomial)
  #plot(gam_Dist_P_npoS, pages=1)
  gam_Dist_N_npoS <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Dist_npo[dat_hist_Dist_npo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_npoS, pages=1)
  
  dat_hist$presxD_npoS <- predict(gam_Dist_P_npoS, dat_hist, type="response")
  abundxD_npoS <- predict(gam_Dist_N_npoS, dat_hist, type="response")
  dat_hist$gam_Dist_npo_S <- dat_hist$presxD_npoS * exp(abundxD_npoS)
  
  dat_fcast$presxD_npoS <- predict(gam_Dist_P_npo, dat_fcast, type="response")
  abundxD_npoS <- predict(gam_Dist_N_npo, dat_fcast, type="response")
  dat_fcast$gam_Dist_npo_S <- dat_fcast$presxD_npoS * exp(abundxD_npoS)
}
#Dist sampling - NPN
if("npn" %in% sampling){
  print("Fitting GAM-NPN_S")
  gam_Dist_P_npnS <- gam(pres~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Dist_npn, family=binomial)
  #plot(gam_Dist_P_npnS, pages=1)
  gam_Dist_N_npnS <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Dist_npn[dat_hist_Dist_npn$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_npnS, pages=1)
  
  dat_hist$presxD_npnS <- predict(gam_Dist_P_npnS, dat_hist, type="response")
  abundxD_npnS <- predict(gam_Dist_N_npnS, dat_hist, type="response")
  dat_hist$gam_Dist_npn_S <- dat_hist$presxD_npnS * exp(abundxD_npnS)
  
  dat_fcast$presxD_npnS <- predict(gam_Dist_P_npnS, dat_fcast, type="response")
  abundxD_npnS <- predict(gam_Dist_N_npnS, dat_fcast, type="response")
  dat_fcast$gam_Dist_npn_S <- dat_fcast$presxD_npnS * exp(abundxD_npnS)
}
#Dist sampling - MPO
if("mpo" %in% sampling){
  print("Fitting GAM-MPO_S")
  gam_Dist_P_mpoS <- gam(pres~s(temp) + s(mld) + s(chl_surface)+ s(lat,lon), data=dat_hist_Dist_mpo, family=binomial)
  #plot(gam_Dist_P_mpoS, pages=1)
  gam_Dist_N_mpoS <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface)+ s(lat,lon), data=dat_hist_Dist_mpo[dat_hist_Dist_mpo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_mpoS, pages=1)
  
  dat_hist$presxD_mpoS <- predict(gam_Dist_P_mpoS, dat_hist, type="response")
  abundxD_mpoS <- predict(gam_Dist_N_mpoS, dat_hist, type="response")
  dat_hist$gam_Dist_mpo_S <- dat_hist$presxD_mpoS * exp(abundxD_mpoS)
  
  dat_fcast$presxD_mpoS <- predict(gam_Dist_P_mpoS, dat_fcast, type="response")
  abundxD_mpoS <- predict(gam_Dist_N_mpoS, dat_fcast, type="response")
  dat_fcast$gam_Dist_mpo_S <- dat_fcast$presxD_mpoS * exp(abundxD_mpoS)
}
#Dist sampling - MPN
if("mpn" %in% sampling){
  print("Fitting GAM-MPN_S")
  gam_Dist_P_mpnS <- gam(pres~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Dist_mpn, family=binomial)
  #plot(gam_Dist_P_mpnS, pages=1)
  gam_Dist_N_mpnS <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Dist_mpn[dat_hist_Dist_mpn$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_mpnS, pages=1)
  
  dat_hist$presxD_mpnS <- predict(gam_Dist_P_mpnS, dat_hist, type="response")
  abundxD_mpnS <- predict(gam_Dist_N_mpnS, dat_hist, type="response")
  dat_hist$gam_Dist_mpn_S <- dat_hist$presxD_mpnS * exp(abundxD_mpnS)
  
  dat_fcast$presxD_mpnS <- predict(gam_Dist_P_mpnS, dat_fcast, type="response")
  abundxD_mpnS <- predict(gam_Dist_N_mpnS, dat_fcast, type="response")
  dat_fcast$gam_Dist_mpn_S <- dat_fcast$presxD_mpnS * exp(abundxD_mpnS)
}
#Dist sampling - SPO
if("spo" %in% sampling){
  print("Fitting GAM-SPO_S")
  gam_Dist_P_spoS <- gam(pres~s(temp) + s(mld) + s(chl_surface)+ s(lat,lon), data=dat_hist_Dist_spo, family=binomial)
  #plot(gam_Dist_P_spoS, pages=1)
  gam_Dist_N_spoS <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface)+ s(lat,lon), data=dat_hist_Dist_spo[dat_hist_Dist_spo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_spoS, pages=1)
  
  dat_hist$presxD_spoS <- predict(gam_Dist_P_spoS, dat_hist, type="response")
  abundxD_spoS <- predict(gam_Dist_N_spoS, dat_hist, type="response")
  dat_hist$gam_Dist_spo_S <- dat_hist$presxD_spoS * exp(abundxD_spoS)
  
  dat_fcast$presxD_spoS <- predict(gam_Dist_P_spoS, dat_fcast, type="response")
  abundxD_spoS <- predict(gam_Dist_N_spoS, dat_fcast, type="response")
  dat_fcast$gam_Dist_spo_S <- dat_fcast$presxD_spoS * exp(abundxD_spoS)
}
#Dist sampling - SPN
if("spn" %in% sampling){
  print("Fitting GAM-SPN_S")
  gam_Dist_P_spnS <- gam(pres~s(temp) + s(mld) + s(chl_surface)+ s(lat,lon), data=dat_hist_Dist_spn, family=binomial)
  #plot(gam_Dist_P_spnS, pages=1)
  gam_Dist_N_spnS <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface)+ s(lat,lon), data=dat_hist_Dist_spn[dat_hist_Dist_spn$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_spnS, pages=1)
  
  dat_hist$presxD_spnS <- predict(gam_Dist_P_spnS, dat_hist, type="response")
  abundxD_spnS <- predict(gam_Dist_N_spnS, dat_hist, type="response")
  dat_hist$gam_Dist_spn_S <- dat_hist$presxD_spnS * exp(abundxD_spnS)
  
  dat_fcast$presxD_spnS <- predict(gam_Dist_P_spnS, dat_fcast, type="response")
  abundxD_spnS <- predict(gam_Dist_N_spnS, dat_fcast, type="response")
  dat_fcast$gam_Dist_spn_S <- dat_fcast$presxD_spnS * exp(abundxD_spnS)
}
#Dist sampling - ALLO
if("allo" %in% sampling){
  print("Fitting GAM-ALLO_S")
  gam_Dist_P_alloS <- gam(pres~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Dist_allo, family=binomial)
  #plot(gam_Dist_P_alloS, pages=1)
  gam_Dist_N_alloS <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Dist_allo[dat_hist_Dist_allo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_alloS, pages=1)
  
  dat_hist$presxD_alloS <- predict(gam_Dist_P_alloS, dat_hist, type="response")
  abundxD_alloS <- predict(gam_Dist_N_alloS, dat_hist, type="response")
  dat_hist$gam_Dist_allo_S <- dat_hist$presxD_alloS * exp(abundxD_alloS)
  
  dat_fcast$presxD_alloS <- predict(gam_Dist_P_alloS, dat_fcast, type="response")
  abundxD_alloS <- predict(gam_Dist_N_alloS, dat_fcast, type="response")
  dat_fcast$gam_Dist_allo_S <- dat_fcast$presxD_alloS * exp(abundxD_alloS)
}
#Dist sampling - Alln
if("alln" %in% sampling){
  print("Fitting GAM-ALLN_S")
  gam_Dist_P_allnS <- gam(pres~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Dist_alln, family=binomial)
  #plot(gam_Dist_P_allnS, pages=1)
  gam_Dist_N_allnS <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_Dist_alln[dat_hist_Dist_alln$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_allnS, pages=1)
  
  dat_hist$presxD_allnS <- predict(gam_Dist_P_allnS, dat_hist, type="response")
  abundxD_allnS <- predict(gam_Dist_N_allnS, dat_hist, type="response")
  dat_hist$gam_Dist_alln_S <- dat_hist$presxD_allnS * exp(abundxD_allnS)
  
  dat_fcast$presxD_allnS <- predict(gam_Dist_P_allnS, dat_fcast, type="response")
  abundxD_allnS <- predict(gam_Dist_N_allnS, dat_fcast, type="response")
  dat_fcast$gam_Dist_alln_S <- dat_fcast$presxD_allnS * exp(abundxD_allnS)
}
#
######BY sampling
#
if("BY" %in% sampling){
  print("Fitting GAM-BY_S")
  gam_BY_PS <- gam(pres~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_BY, family=binomial)
  #plot(gam_BY_PS, pages=1)
  gam_BY_NS <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_BY[dat_hist_BY$abundance>0,], family=gaussian)
  #plot(gam_BY_NS, pages=1)
  
  dat_hist$presxBS <- predict(gam_BY_PS, dat_hist, type="response")
  abundxBS <- predict(gam_BY_NS, dat_hist, type="response")
  dat_hist$gam_BY_S <- dat_hist$presxBS * exp(abundxBS)
  
  dat_fcast$presxBS <- predict(gam_BY_PS, dat_fcast, type="response")
  abundxBS <- predict(gam_BY_NS, dat_fcast, type="response")
  dat_fcast$gam_BY_S <- dat_fcast$presxBS * exp(abundxBS)
}

#
#####Closed Area Sampling
#
#Closed Area - Small
if("CA_sm" %in% sampling) {
  print("Fitting GAM-CASM_S")
  gam_CA_P_smS <- gam(pres~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_CA_sm, family=binomial)
  #plot(gam_CA_P_smS, pages=1)
  gam_CA_N_smS <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_CA_sm[dat_hist_CA_sm$abundance>0,], family=gaussian)
  #plot(gam_CA_N_smS, pages=1)
  
  dat_hist$presxCA_smS <- predict(gam_CA_P_smS, dat_hist, type="response")
  abundxCA_smS <- predict(gam_CA_N_smS, dat_hist, type="response")
  dat_hist$gam_CA_sm_S <- dat_hist$presxCA_smS * exp(abundxCA_smS)
  
  dat_fcast$presxCA_smS <- predict(gam_CA_P_smS, dat_fcast, type="response")
  abundxCA_smS <- predict(gam_CA_N_smS, dat_fcast, type="response")
  dat_fcast$gam_CA_sm_S <- dat_fcast$presxCA_smS * exp(abundxCA_smS)
}
#Closed Area - Medium
if("CA_med" %in% sampling) {
  print("Fitting GAM-CAMED_S")
  gam_CA_P_medS <- gam(pres~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_CA_med, family=binomial)
  #plot(gam_CA_P_medS, pages=1)
  gam_CA_N_medS <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_CA_med[dat_hist_CA_med$abundance>0,], family=gaussian)
  #plot(gam_CA_N_medS, pages=1)
  
  dat_hist$presxCA_medS <- predict(gam_CA_P_medS, dat_hist, type="response")
  abundxCA_medS <- predict(gam_CA_N_medS, dat_hist, type="response")
  dat_hist$gam_CA_med_S <- dat_hist$presxCA_medS * exp(abundxCA_medS)
  
  dat_fcast$presxCA_medS <- predict(gam_CA_P_medS, dat_fcast, type="response")
  abundxCA_medS <- predict(gam_CA_N_medS, dat_fcast, type="response")
  dat_fcast$gam_CA_med_S <- dat_fcast$presxCA_medS * exp(abundxCA_medS)
}
#Closed Area - Large
if("CA_lar" %in% sampling) {
  print("Fitting GAM-CALAR_S")
  gam_CA_P_larS <- gam(pres~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_CA_lar, family=binomial)
  #plot(gam_CA_P_larS, pages=1)
  gam_CA_N_larS <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + s(lat,lon), data=dat_hist_CA_lar[dat_hist_CA_lar$abundance>0,], family=gaussian)
  #plot(gam_CA_N_larS, pages=1)
  
  dat_hist$presxCA_larS <- predict(gam_CA_P_larS, dat_hist, type="response")
  abundxCA_larS <- predict(gam_CA_N_larS, dat_hist, type="response")
  dat_hist$gam_CA_lar_S <- dat_hist$presxCA_larS * exp(abundxCA_larS)
  
  dat_fcast$presxCA_larS <- predict(gam_CA_P_larS, dat_fcast, type="response")
  abundxCA_larS <- predict(gam_CA_N_larS, dat_fcast, type="response")
  dat_fcast$gam_CA_lar_S <- dat_fcast$presxCA_larS * exp(abundxCA_larS)
}
