### Fishery-dependent SDM (L^3) Estimation Model Code
## Function to run GAMs on data from OM -- with missing chla covariate & some w/ spatial...
## Function returns only the fitted and predicted values
#modified my M.Karp 1/04/21

############################################################
#   Partial_models without Chla and with Space-Time tensor #
############################################################

#
##### Random sampling ####
#
if("ran" %in% sampling){
  print("Fitting GAM-Ran_nochl_te")
  gam_Ran_P_nochlte <- gam(pres~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_random, family=binomial) 
  #plot(gam_Ran_P_nochlte, pages=1)
  gam_Ran_N_nochlte <- gam(log_abundance~s(temp)+ s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_random[dat_hist_random$abundance>0,], family=gaussian)
  #plot(gam_Ran_N_nochlte, pages=1)
  
  dat_hist$presxR_nochlte <- predict(gam_Ran_P_nochlte, dat_hist, type="response")
  abundxR_nochlte <- predict(gam_Ran_N_nochlte, dat_hist, type="response") 
  dat_hist$gam_Ran_nochl_te <- dat_hist$presxR_nochlte * exp(abundxR_nochlte)  #predicted catch
  
  dat_fcast$presxR_nochlte <- predict(gam_Ran_P_nochlte, dat_fcast, type="response")
  abundxR_nochlte <- predict(gam_Ran_N_nochlte, dat_fcast, type="response")
  dat_fcast$gam_Ran_nochl_te <- dat_fcast$presxR_nochlte * exp(abundxR_nochlte)
}

#
##### Preferential Sampling #######
#
#Pref sampling - 0.5
if("tar_0.5" %in% sampling) {
  print("Fitting GAM-Pref0.5_nochl_te")
  gam_Tar_P_1_nochlte <- gam(pres~s(temp)+ s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Tar_1, family=binomial)
  #plot(gam_Tar_P_1_nochlte, pages=1)
  gam_Tar_N_1_nochlte <- gam(log_abundance~s(temp)+ s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Tar_1[dat_hist_Tar_1$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_1_nochlte, pages=1)
  
  dat_hist$presxT_1_nochlte <- predict(gam_Tar_P_1_nochlte, dat_hist, type="response")
  abundxT_1_nochlte <- predict(gam_Tar_N_1_nochlte, dat_hist, type="response")
  dat_hist$gam_Tar_0.5_nochl_te <- dat_hist$presxT_1_nochlte * exp(abundxT_1_nochlte)
  
  dat_fcast$presxT_1_nochlte <- predict(gam_Tar_P_1_nochlte, dat_fcast, type="response")
  abundxT_1_nochlte <- predict(gam_Tar_N_1_nochlte, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.5_nochl_te <- dat_fcast$presxT_1_nochlte * exp(abundxT_1_nochlte)
}
#Pref sampling -0.6
if("tar_0.6" %in% sampling) {
  print("Fitting GAM-Pref0.6_nochl_te")
  gam_Tar_P_2_nochlte <- gam(pres~s(temp)+ s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Tar_2, family=binomial)
  #plot(gam_Tar_P_2_nochlte, pages=1)
  gam_Tar_N_2_nochlte <- gam(log_abundance~s(temp)+ s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Tar_2[dat_hist_Tar_2$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_2_nochlte, pages=1)
  
  dat_hist$presxT_2_nochlte <- predict(gam_Tar_P_2_nochlte, dat_hist, type="response")
  abundxT_2_nochlte <- predict(gam_Tar_N_2_nochlte, dat_hist, type="response")
  dat_hist$gam_Tar_0.6_nochl_te <- dat_hist$presxT_2_nochlte * exp(abundxT_2_nochlte)
  
  dat_fcast$presxT_2_nochlte <- predict(gam_Tar_P_2_nochlte, dat_fcast, type="response")
  abundxT_2_nochlte <- predict(gam_Tar_N_2_nochlte, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.6_nochl_te <- dat_fcast$presxT_2_nochlte * exp(abundxT_2_nochlte)
}
#Pref sampling -0.7
if("tar_0.7" %in% sampling) {
  print("Fitting GAM-Pref0.7_nochl_te")
  gam_Tar_P_3_nochlte <- gam(pres~s(temp)+ s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Tar_3, family=binomial)
  #plot(gam_Tar_P_3_nochlte, pages=1)
  gam_Tar_N_3_nochlte <- gam(log_abundance~s(temp)+ s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Tar_3[dat_hist_Tar_3$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_3_nochlte, pages=1)
  
  dat_hist$presxT_3_nochlte <- predict(gam_Tar_P_3_nochlte, dat_hist, type="response")
  abundxT_3_nochlte <- predict(gam_Tar_N_3_nochlte, dat_hist, type="response")
  dat_hist$gam_Tar_0.7_nochl_te <- dat_hist$presxT_3_nochlte * exp(abundxT_3_nochlte)
  
  dat_fcast$presxT_3_nochlte <- predict(gam_Tar_P_3_nochlte, dat_fcast, type="response")
  abundxT_3_nochlte <- predict(gam_Tar_N_3_nochlte, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.7_nochl_te <- dat_fcast$presxT_3_nochlte * exp(abundxT_3_nochlte)
}
#Pref sampling -0.8
if("tar_0.8" %in% sampling) {
  print("Fitting GAM-Pref0.8_nochl_te")
  gam_Tar_P_4_nochlte <- gam(pres~s(temp)+ s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Tar_4, family=binomial)
  #plot(gam_Tar_P_4_nochlte, pages=1)
  gam_Tar_N_4_nochlte <- gam(log_abundance~s(temp)+ s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Tar_4[dat_hist_Tar_4$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_4_nochlte, pages=1)
  
  dat_hist$presxT_4_nochlte <- predict(gam_Tar_P_4_nochlte, dat_hist, type="response")
  abundxT_4_nochlte <- predict(gam_Tar_N_4_nochlte, dat_hist, type="response")
  dat_hist$gam_Tar_0.8_nochl_te <- dat_hist$presxT_4_nochlte * exp(abundxT_4_nochlte)
  
  dat_fcast$presxT_4_nochlte <- predict(gam_Tar_P_4_nochlte, dat_fcast, type="response")
  abundxT_4_nochlte <- predict(gam_Tar_N_4_nochlte, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.8_nochl_te <- dat_fcast$presxT_4_nochlte * exp(abundxT_4_nochlte)
}
#Pref sampling -0.9
if("tar_0.9" %in% sampling) {
  print("Fitting GAM-Pref0.9_nochl_te")
  gam_Tar_P_5_nochlte <- gam(pres~s(temp)+ s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Tar_5, family=binomial)
  #plot(gam_Tar_P_5_nochlte, pages=1)
  gam_Tar_N_5_nochlte <- gam(log_abundance~s(temp)+ s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Tar_5[dat_hist_Tar_5$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_5_nochlte, pages=1)
  
  dat_hist$presxT_5_nochlte <- predict(gam_Tar_P_5_nochlte, dat_hist, type="response")
  abundxT_5_nochlte <- predict(gam_Tar_N_5_nochlte, dat_hist, type="response")
  dat_hist$gam_Tar_0.9_nochl_te <- dat_hist$presxT_5_nochlte* exp(abundxT_5_nochlte)
  
  dat_fcast$presxT_5_nochlte <- predict(gam_Tar_P_5_nochlte, dat_fcast, type="response")
  abundxT_5_nochlte <- predict(gam_Tar_N_5_nochlte, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.9_nochl_te <- dat_fcast$presxT_5_nochlte * exp(abundxT_5_nochlte)
}

#
##### Distance from Port Sampling #######
#
#Dist sampling - NPO
if("npo" %in% sampling){
  print("Fitting GAM-NPO_nochl_te")
  gam_Dist_P_npo_nochlte <- gam(pres~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_npo, family=binomial)
  #plot(gam_Dist_P_npo_nochlte, pages=1)
  gam_Dist_N_npo_nochlte <- gam(log_abundance~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_npo[dat_hist_Dist_npo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_npo_nochlte, pages=1)
  
  dat_hist$presxD_npo_nochlte <- predict(gam_Dist_P_npo_nochlte, dat_hist, type="response")
  abundxD_npo_nochlte <- predict(gam_Dist_N_npo_nochlte, dat_hist, type="response")
  dat_hist$gam_Dist_npo_nochl_te <- dat_hist$presxD_npo_nochlte * exp(abundxD_npo_nochlte)
  
  dat_fcast$presxD_npo_nochlte <- predict(gam_Dist_P_npo_nochlte, dat_fcast, type="response")
  abundxD_npo_nochlte <- predict(gam_Dist_N_npo_nochlte, dat_fcast, type="response")
  dat_fcast$gam_Dist_npo_nochl_te <- dat_fcast$presxD_npo_nochlte * exp(abundxD_npo_nochlte)
}
#Dist sampling - NPN
if("npn" %in% sampling){
  print("Fitting GAM-NPN_nochl_te")
  gam_Dist_P_npn_nochlte <- gam(pres~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_npn, family=binomial)
  #plot(gam_Dist_P_npn_nochlte, pages=1)
  gam_Dist_N_npn_nochlte <- gam(log_abundance~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_npn[dat_hist_Dist_npn$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_npn_nochlte, pages=1)
  
  dat_hist$presxD_npn_nochlte <- predict(gam_Dist_P_npn_nochlte, dat_hist, type="response")
  abundxD_npn_nochlte <- predict(gam_Dist_N_npn_nochlte, dat_hist, type="response")
  dat_hist$gam_Dist_npn_nochl_te <- dat_hist$presxD_npn_nochlte * exp(abundxD_npn_nochlte)
  
  dat_fcast$presxD_npn_nochlte <- predict(gam_Dist_P_npn_nochlte, dat_fcast, type="response")
  abundxD_npn_nochlte <- predict(gam_Dist_N_npn_nochlte, dat_fcast, type="response")
  dat_fcast$gam_Dist_npn_nochl_te <- dat_fcast$presxD_npn_nochlte * exp(abundxD_npn_nochlte)
}

#Dist sampling - MPO
if("mpo" %in% sampling){
  print("Fitting GAM-MPO_nochl_te")
  gam_Dist_P_mpo_nochlte <- gam(pres~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_mpo, family=binomial)
  #plot(gam_Dist_P_mpo_nochlte, pages=1)
  gam_Dist_N_mpo_nochlte <- gam(log_abundance~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_mpo[dat_hist_Dist_mpo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_mpo_nochlte, pages=1)
  
  dat_hist$presxD_mpo_nochlte <- predict(gam_Dist_P_mpo_nochlte, dat_hist, type="response")
  abundxD_mpo_nochlte <- predict(gam_Dist_N_mpo_nochlte, dat_hist, type="response")
  dat_hist$gam_Dist_mpo_nochl_te <- dat_hist$presxD_mpo_nochlte * exp(abundxD_mpo_nochlte)
  
  dat_fcast$presxD_mpo_nochlte <- predict(gam_Dist_P_mpo_nochlte, dat_fcast, type="response")
  abundxD_mpo_nochlte <- predict(gam_Dist_N_mpo_nochlte, dat_fcast, type="response")
  dat_fcast$gam_Dist_mpo_nochl_te <- dat_fcast$presxD_mpo_nochlte * exp(abundxD_mpo_nochlte)
}
#Dist sampling - MPN
if("mpn" %in% sampling){
  print("Fitting GAM-NPN_nochl_te")
  gam_Dist_P_mpn_nochlte <- gam(pres~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_mpn, family=binomial)
  #plot(gam_Dist_P_mpn_nochlte, pages=1)
  gam_Dist_N_mpn_nochlte <- gam(log_abundance~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_mpn[dat_hist_Dist_mpn$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_mpn_nochlte, pages=1)
  
  dat_hist$presxD_mpn_nochlte <- predict(gam_Dist_P_mpn_nochlte, dat_hist, type="response")
  abundxD_mpn_nochlte <- predict(gam_Dist_N_mpn_nochlte, dat_hist, type="response")
  dat_hist$gam_Dist_mpn_nochl_te <- dat_hist$presxD_mpn_nochlte * exp(abundxD_mpn_nochlte)
  
  dat_fcast$presxD_mpn_nochlte <- predict(gam_Dist_P_mpn_nochlte, dat_fcast, type="response")
  abundxD_mpn_nochlte <- predict(gam_Dist_N_mpn_nochlte, dat_fcast, type="response")
  dat_fcast$gam_Dist_mpn_nochl_te <- dat_fcast$presxD_mpn_nochlte * exp(abundxD_mpn_nochlte)
}
#Dist sampling - SPO
if("spo" %in% sampling){
  print("Fitting GAM-SPO_nochl_te")
  gam_Dist_P_spo_nochlte <- gam(pres~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_spo, family=binomial)
  #plot(gam_Dist_P_spo_nochlte, pages=1)
  gam_Dist_N_spo_nochlte <- gam(log_abundance~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_spo[dat_hist_Dist_spo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_spo_nochlte, pages=1)
  
  dat_hist$presxD_spo_nochlte <- predict(gam_Dist_P_spo_nochlte, dat_hist, type="response")
  abundxD_spo_nochlte <- predict(gam_Dist_N_spo_nochlte, dat_hist, type="response")
  dat_hist$gam_Dist_spo_nochl_te <- dat_hist$presxD_spo_nochlte * exp(abundxD_spo_nochlte)
  
  dat_fcast$presxD_spo_nochlte <- predict(gam_Dist_P_spo_nochlte, dat_fcast, type="response")
  abundxD_spo_nochlte <- predict(gam_Dist_N_spo_nochlte, dat_fcast, type="response")
  dat_fcast$gam_Dist_spo_nochl_te <- dat_fcast$presxD_spo_nochlte * exp(abundxD_spo_nochlte)
}
#Dist sampling - SPN
if("spn" %in% sampling){
  print("Fitting GAM-SPN_nochl_te")
  gam_Dist_P_spn_nochlte <- gam(pres~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_spn, family=binomial)
  #plot(gam_Dist_P_spn_nochlte, pages=1)
  gam_Dist_N_spn_nochlte <- gam(log_abundance~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_spn[dat_hist_Dist_spn$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_spn_nochlte, pages=1)
  
  dat_hist$presxD_spn_nochlte <- predict(gam_Dist_P_spn_nochlte, dat_hist, type="response")
  abundxD_spn_nochlte <- predict(gam_Dist_N_spn_nochlte, dat_hist, type="response")
  dat_hist$gam_Dist_spn_nochl_te <- dat_hist$presxD_spn_nochlte * exp(abundxD_spn_nochlte)
  
  dat_fcast$presxD_spn_nochlte <- predict(gam_Dist_P_spn_nochlte, dat_fcast, type="response")
  abundxD_spn_nochlte <- predict(gam_Dist_N_spn_nochlte, dat_fcast, type="response")
  dat_fcast$gam_Dist_spn_nochl_te <- dat_fcast$presxD_spn_nochlte * exp(abundxD_spn_nochlte)
}
#Dist sampling - ALLO
if("allo" %in% sampling){
  print("Fitting GAM-ALLO_nochl_te")
  gam_Dist_P_allo_nochlte <- gam(pres~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_allo, family=binomial)
  #plot(gam_Dist_P_allo_nochlte, pages=1)
  gam_Dist_N_allo_nochlte <- gam(log_abundance~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_allo[dat_hist_Dist_allo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_allo_nochlte, pages=1)
  
  dat_hist$presxD_allo_nochlte <- predict(gam_Dist_P_allo_nochlte, dat_hist, type="response")
  abundxD_allo_nochlte <- predict(gam_Dist_N_allo_nochlte, dat_hist, type="response")
  dat_hist$gam_Dist_allo_nochl_te <- dat_hist$presxD_allo_nochlte * exp(abundxD_allo_nochlte)
  
  dat_fcast$presxD_allo_nochlte <- predict(gam_Dist_P_allo_nochlte, dat_fcast, type="response")
  abundxD_allo_nochlte <- predict(gam_Dist_N_allo_nochlte, dat_fcast, type="response")
  dat_fcast$gam_Dist_allo_nochl_te <- dat_fcast$presxD_allo_nochlte * exp(abundxD_allo_nochlte)
}
#Dist sampling - Alln
if("alln" %in% sampling){
  print("Fitting GAM-ALLN_nochl_te")
  gam_Dist_P_alln_nochlte <- gam(pres~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_alln, family=binomial)
  #plot(gam_Dist_P_alln_nochte, pages=1)
  gam_Dist_N_alln_nochlte <- gam(log_abundance~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_Dist_alln[dat_hist_Dist_alln$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_alln_nochlte, pages=1)
  
  dat_hist$presxD_alln_nochlte <- predict(gam_Dist_P_alln_nochlte, dat_hist, type="response")
  abundxD_alln_nochlte <- predict(gam_Dist_N_alln_nochlte, dat_hist, type="response")
  dat_hist$gam_Dist_alln_nochl_te <- dat_hist$presxD_alln_nochlte * exp(abundxD_alln_nochlte)
  
  dat_fcast$presxD_alln_nochlte <- predict(gam_Dist_P_alln_nochlte, dat_fcast, type="response")
  abundxD_alln_nochlte <- predict(gam_Dist_N_alln_nochlte, dat_fcast, type="response")
  dat_fcast$gam_Dist_alln_nochl_te <- dat_fcast$presxD_alln_nochlte * exp(abundxD_alln_nochlte)
}

#
######BY sampling
#
if("BY" %in% sampling){
  print("Fitting GAM-BY_nochl_te")
  gam_BY_P_nochlte <- gam(pres~s(temp) + s(mld)+ te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_BY, family=binomial)
  #plot(gam_BY_P_nochlte, pages=1)
  gam_BY_N_nochlte<- gam(log_abundance~s(temp) + s(mld)+ te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_BY[dat_hist_BY$abundance>0,], family=gaussian)
  #plot(gam_BY_N_nochlte, pages=1)
  
  dat_hist$presxB_nochlte <- predict(gam_BY_P_nochlte, dat_hist, type="response")
  abundxB_nochlte <- predict(gam_BY_N_nochlte, dat_hist, type="response")
  dat_hist$gam_BY_nochl_te <- dat_hist$presxB_nochlte * exp(abundxB_nochlte)
  
  dat_fcast$presxB_nochlte <- predict(gam_BY_P_nochlte, dat_fcast, type="response")
  abundxB_nochlte <- predict(gam_BY_N_nochlte, dat_fcast, type="response")
  dat_fcast$gam_BY_nochl_te <- dat_fcast$presxB_nochlte * exp(abundxB_nochlte)
}

#
#####Closed Area Sampling
#
#Closed Area - Small
if("CA_sm" %in% sampling) {
  print("Fitting GAM-CASM_nochl_te")
  gam_CA_P_sm_nochlte <- gam(pres~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_CA_sm, family=binomial)
  #plot(gam_CA_P_sm_nochlte, pages=1)
  gam_CA_N_sm_nochlte <- gam(log_abundance~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_CA_sm[dat_hist_CA_sm$abundance>0,], family=gaussian)
  #plot(gam_CA_N_sm_nochlte, pages=1)
  
  dat_hist$presxCA_sm_nochlte <- predict(gam_CA_P_sm_nochlte, dat_hist, type="response")
  abundxCA_sm_nochlte <- predict(gam_CA_N_sm_nochlte, dat_hist, type="response")
  dat_hist$gam_CA_sm_nochl_te <- dat_hist$presxCA_sm_nochlte * exp(abundxCA_sm_nochlte)
  
  dat_fcast$presxCA_sm_nochlte <- predict(gam_CA_P_sm_nochlte, dat_fcast, type="response")
  abundxCA_sm_nochlte <- predict(gam_CA_N_sm_nochlte, dat_fcast, type="response")
  dat_fcast$gam_CA_sm_nochl_te <- dat_fcast$presxCA_sm_nochlte * exp(abundxCA_sm_nochlte)
}
#Closed Area - Medium
if("CA_med" %in% sampling) {
  print("Fitting GAM-CAMED_nochl_te")
  gam_CA_P_med_nochlte <- gam(pres~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_CA_med, family=binomial)
  #plot(gam_CA_P_med_nochlte, pages=1)
  gam_CA_N_med_nochlte <- gam(log_abundance~s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_CA_med[dat_hist_CA_med$abundance>0,], family=gaussian)
  #plot(gam_CA_N_med_nochlte, pages=1)
  
  dat_hist$presxCA_med_nochlte <- predict(gam_CA_P_med_nochlte, dat_hist, type="response")
  abundxCA_med_nochlte <- predict(gam_CA_N_med_nochlte, dat_hist, type="response")
  dat_hist$gam_CA_med_nochl_te <- dat_hist$presxCA_med_nochlte * exp(abundxCA_med_nochlte)
  
  dat_fcast$presxCA_med_nochlte <- predict(gam_CA_P_med_nochlte, dat_fcast, type="response")
  abundxCA_med_nochlte <- predict(gam_CA_N_med_nochlte, dat_fcast, type="response")
  dat_fcast$gam_CA_med_nochl_te <- dat_fcast$presxCA_med_nochlte * exp(abundxCA_med_nochlte)
}
#Closed Area - Large
if("CA_lar" %in% sampling) {
  print("Fitting GAM-CALAR_nochl_te")
  gam_CA_P_lar_nochlte <- gam(pres~ s(temp) + s(mld) + te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_CA_lar, family=binomial)
  #plot(gam_CA_P_lar_nochlte, pages=1)
  gam_CA_N_lar_nochlte <- gam(log_abundance~ s(temp) + s(mld) +  te(lat, lon, year, d=c(2,1), bs=c('ds','cr'), m=c(1,.5), k=c(30,10)), data=dat_hist_CA_lar[dat_hist_CA_lar$abundance>0,], family=gaussian)
  #plot(gam_CA_N_lar_nochlte, pages=1)
  
  dat_hist$presxCA_lar_nochlte <- predict(gam_CA_P_lar_nochlte, dat_hist, type="response")
  abundxCA_lar_nochlte <- predict(gam_CA_N_lar_nochlte, dat_hist, type="response")
  dat_hist$gam_CA_lar_nochl_te <- dat_hist$presxCA_lar_nochlte * exp(abundxCA_lar_nochlte)
  
  dat_fcast$presxCA_lar_nochlte <- predict(gam_CA_P_lar_nochlte, dat_fcast, type="response")
  abundxCA_lar_nochlte <- predict(gam_CA_N_lar_nochlte, dat_fcast, type="response")
  dat_fcast$gam_CA_lar_nochl_te<- dat_fcast$presxCA_lar_nochlte * exp(abundxCA_lar_nochlte)
}