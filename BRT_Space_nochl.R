### Fishery-dependent SDM (L^3) Estimation Model Code
## Function to run BRTs on data from OM, with all Covariates (temp, chla, mld)
## Function returns only the fitted and predicted values
## modified my M.Karp 1/15/21

#############################
#   Full Models- with Chla and Space component #
#############################
#
###Random sampling
#
if("ran" %in% sampling) {
  print("Fitting BRT-Ran_nochl_S")
  brt_R_P_nochl_S <- gbm.step(data=dat_hist_random,
                        gbm.x = c(1:2, 25, 27),
                        gbm.y = 'pres',
                        family = "bernoulli",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  brt_R_N_nochl_S <- gbm.step(data=dat_hist_random[dat_hist_random$abundance>0,], 
                        gbm.x = c(1:2, 25, 27),
                        gbm.y = 'log_abundance',
                        family = "gaussian",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  presRx_nochl_nochl_S <- predict(brt_R_P_nochl_nochl_S, dat_hist, n.trees=brt_R_P_nochl_nochl_S$gbm.call$best.trees, type="response")
  abundRx_nochl_nochl_S <- exp(predict(brt_R_N_nochl_nochl_S, dat_hist, n.trees=brt_R_N_nochl_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Ran_nochl__S <- presRx_nochl_nochl_S * abundRx_nochl_nochl_S
  
  presRx_nochl_S <- predict(brt_R_P_nochl_S, dat_fcast, n.trees=brt_R_P_nochl_S$gbm.call$best.trees, type="response")
  abundRx_nochl_S <- exp(predict(brt_R_N_nochl_S, dat_fcast, n.trees=brt_R_N_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Ran_nochl_S <- presRx_nochl_S * abundRx_nochl_S
}

#
### Preferential (target) sampling
#
#Target sampling-0.5
if("tar_0.5" %in% sampling){
  print("Fitting BRT-Pref0.5_nochl_S")
  brt_T_P_1_nochl_S <- gbm.step(data=dat_hist_Tar_1,
                          gbm.x = c(1:2, 25, 27),
                          gbm.y = 'pres',
                          family = "bernoulli",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_1_nochl_S <- gbm.step(data=dat_hist_Tar_1[dat_hist_Tar_1$abundance>0,],
                          gbm.x = c(1:2, 25, 27),
                          gbm.y = 'log_abundance',
                          family = "gaussian",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  presTx_1_nochl_S <- predict(brt_T_P_1_nochl_S, dat_hist, n.trees=brt_T_P_1_nochl_S$gbm.call$best.trees, type="response")
  abundTx_1_nochl_S <- exp(predict(brt_T_N_1_nochl_S, dat_hist, n.trees=brt_T_N_1_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.5_nochl_S <- presTx_1_nochl_S * abundTx_1_nochl_S
  
  presTx_1_nochl_S <- predict(brt_T_P_1_nochl_S, dat_fcast, n.trees=brt_T_P_1_nochl_S$gbm.call$best.trees, type="response")
  abundTx_1_nochl_S <- exp(predict(brt_T_N_1_nochl_S, dat_fcast, n.trees=brt_T_N_1_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.5_nochl_S <- presTx_1_nochl_S * abundTx_1_nochl_S
}

#Target sampling-0.6
if("tar_0.6" %in% sampling){
  print("Fitting BRT-Pref0.6_nochl_S")
  brt_T_P_2_nochl_S <- gbm.step(data=dat_hist_Tar_2,
                          gbm.x = c(1:2, 25, 27),
                          gbm.y = 'pres',
                          family = "bernoulli",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_2_nochl_S <- gbm.step(data=dat_hist_Tar_2[dat_hist_Tar_2$abundance>0,],
                          gbm.x = c(1:2, 25, 27),
                          gbm.y = 'log_abundance',
                          family = "gaussian",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  presTx_2_nochl_S <- predict(brt_T_P_2_nochl_S, dat_hist, n.trees=brt_T_P_2_nochl_S$gbm.call$best.trees, type="response")
  abundTx_2_nochl_S <- exp(predict(brt_T_N_2_nochl_S, dat_hist, n.trees=brt_T_N_2_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.6_nochl_S <- presTx_2_nochl_S* abundTx_2_nochl_S
  
  presTx_2_nochl_S <- predict(brt_T_P_2_nochl_S, dat_fcast, n.trees=brt_T_P_2_nochl_S$gbm.call$best.trees, type="response")
  abundTx_2_nochl_S <- exp(predict(brt_T_N_2_nochl_S, dat_fcast, n.trees=brt_T_N_2_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.6_nochl_S <- presTx_2_nochl_S * abundTx_2_nochl_S
}

#Target sampling-0.7
if("tar_0.7" %in% sampling){
  print("Fitting BRT-Pref0.7_nochl_S")
  brt_T_P_3_nochl_S <- gbm.step(data=dat_hist_Tar_3,
                          gbm.x = c(1:2, 25, 27),
                          gbm.y = 'pres',
                          family = "bernoulli",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_3_nochl_S <- gbm.step(data=dat_hist_Tar_3[dat_hist_Tar_3$abundance>0,],
                          gbm.x = c(1:2, 25, 27),
                          gbm.y = 'log_abundance',
                          family = "gaussian",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  presTx_3_nochl_S <- predict(brt_T_P_3_nochl_S, dat_hist, n.trees=brt_T_P_3_nochl_S$gbm.call$best.trees, type="response")
  abundTx_3_nochl_S <- exp(predict(brt_T_N_3_nochl_S, dat_hist, n.trees=brt_T_N_3_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.7_nochl_S <- presTx_3_nochl_S * abundTx_3_nochl_S
  
  presTx_3_nochl_S <- predict(brt_T_P_3_nochl_S, dat_fcast, n.trees=brt_T_P_3_nochl_S$gbm.call$best.trees, type="response")
  abundTx_3_nochl_S <- exp(predict(brt_T_N_3_nochl_S, dat_fcast, n.trees=brt_T_N_3_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.7_nochl_S <- presTx_3_nochl_S * abundTx_3_nochl_S
}

#Target sampling-0.8
if("tar_0.8" %in% sampling){
  print("Fitting BRT-Pref0.8_nochl_S")
  brt_T_P_4_nochl_S <- gbm.step(data=dat_hist_Tar_4,
                          gbm.x = c(1:2, 25, 27),
                          gbm.y = 'pres',
                          family = "bernoulli",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_4_nochl_S <- gbm.step(data=dat_hist_Tar_4[dat_hist_Tar_4$abundance>0,],
                          gbm.x = c(1:2, 25, 27),
                          gbm.y = 'log_abundance',
                          family = "gaussian",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  presTx_4_nochl_S <- predict(brt_T_P_4_nochl_S, dat_hist, n.trees=brt_T_P_4_nochl_S$gbm.call$best.trees, type="response")
  abundTx_4_nochl_S <- exp(predict(brt_T_N_4_nochl_S, dat_hist, n.trees=brt_T_N_4_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.8_nochl_S <- presTx_4_nochl_S * abundTx_4_nochl_S
  
  presTx_4_nochl_S <- predict(brt_T_P_4_nochl_S, dat_fcast, n.trees=brt_T_P_4_nochl_S$gbm.call$best.trees, type="response")
  abundTx_4_nochl_S <- exp(predict(brt_T_N_4_nochl_S, dat_fcast, n.trees=brt_T_N_4_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.8_nochl_S <- presTx_4_nochl_S * abundTx_4_nochl_S
}

#Target sampling-0.9
if("tar_0.9" %in% sampling){
  print("Fitting BRT-Pref0.9_nochl_S")
  brt_T_P_5_nochl_S <- gbm.step(data=dat_hist_Tar_5,
                          gbm.x = c(1:2, 25, 27),
                          gbm.y = 'pres',
                          family = "bernoulli",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_5_nochl_S <- gbm.step(data=dat_hist_Tar_5[dat_hist_Tar_5$abundance>0,],
                          gbm.x = c(1:2, 25, 27),
                          gbm.y = 'log_abundance',
                          family = "gaussian",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  presTx_5_nochl_S <- predict(brt_T_P_5_nochl_S, dat_hist, n.trees=brt_T_P_5_nochl_S$gbm.call$best.trees, type="response")
  abundTx_5_nochl_S <- exp(predict(brt_T_N_5_nochl_S, dat_hist, n.trees=brt_T_N_5_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.9_nochl_S <- presTx_5_nochl_S * abundTx_5_nochl_S
  
  presTx_5_nochl_S <- predict(brt_T_P_5_nochl_S, dat_fcast, n.trees=brt_T_P_5_nochl_S$gbm.call$best.trees, type="response")
  abundTx_5_nochl_S <- exp(predict(brt_T_N_5_nochl_S, dat_fcast, n.trees=brt_T_N_5_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.9_nochl_S <- presTx_5_nochl_S * abundTx_5_nochl_S
}

#
### Distance from Port Sampling
#
#Northern Offshore
if("npo" %in% sampling){
  print("Fitting BRT-NPO_nochl_S")
  brt_dist_P_npo_nochl_S <- gbm.step(data=dat_hist_Dist_npo,
                               gbm.x = c(1:2, 25, 27),
                               gbm.y = 'pres',
                               family = "bernoulli",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_npo_nochl_S <- gbm.step(data=dat_hist_Dist_npo[dat_hist_Dist_npo$abundance>0,],
                               gbm.x = c(1:2, 25, 27),
                               gbm.y = 'log_abundance',
                               family = "gaussian",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  presDx_npo_nochl_S <- predict(brt_dist_P_npo_nochl_S, dat_hist, n.trees=brt_dist_P_npo_nochl_S$gbm.call$best.trees, type="response")
  abundDx_npo_nochl_S <- exp(predict(brt_dist_N_npo_nochl_S, dat_hist, n.trees=brt_dist_N_npo_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_npo_nochl_S <- presDx_npo_nochl_S * abundDx_npo_nochl_S
  
  presDx_npo_nochl_S <- predict(brt_dist_P_npo_nochl_S, dat_fcast, n.trees=brt_dist_P_npo_nochl_S$gbm.call$best.trees, type="response")
  abundDx_npo_nochl_S <- exp(predict(brt_dist_N_npo_nochl_S, dat_fcast, n.trees=brt_dist_N_npo_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_npo_nochl_S <- presDx_npo_nochl_S * abundDx_npo_nochl_S
}

#Northern Nearshore
if("npn" %in% sampling){
  print("Fitting BRT-NPN_nochl_S")
  brt_dist_P_npn_nochl_S <- gbm.step(data=dat_hist_Dist_npn,
                               gbm.x = c(1:2, 25, 27),
                               gbm.y = 'pres',
                               family = "bernoulli",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_npn_nochl_S <- gbm.step(data=dat_hist_Dist_npn[dat_hist_Dist_npn$abundance>0,],
                               gbm.x = c(1:2, 25, 27),
                               gbm.y = 'log_abundance',
                               family = "gaussian",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  presDx_npn_nochl_S <- predict(brt_dist_P_npn_nochl_S, dat_hist, n.trees=brt_dist_P_npn_nochl_S$gbm.call$best.trees, type="response")
  abundDx_npn_nochl_S <- exp(predict(brt_dist_N_npn_nochl_S, dat_hist, n.trees=brt_dist_N_npn_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_npn_nochl_S <- presDx_npn_nochl_S * abundDx_npn_nochl_S
  
  presDx_npn_nochl_S <- predict(brt_dist_P_npn_nochl_S, dat_fcast, n.trees=brt_dist_P_npn_nochl_S$gbm.call$best.trees, type="response")
  abundDx_npn_nochl_S <- exp(predict(brt_dist_N_npn_nochl_S, dat_fcast, n.trees=brt_dist_N_npn_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_npn_nochl_S <- presDx_npn_nochl_S * abundDx_npn_nochl_S
}  

#Middle Offshore
if("mpo" %in% sampling){
  print("Fitting BRT-MPO_nochl_S")
  brt_dist_P_mpo_nochl_S <- gbm.step(data=dat_hist_Dist_mpo,
                               gbm.x = c(1:2, 25, 27),
                               gbm.y = 'pres',
                               family = "bernoulli",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_mpo_nochl_S <- gbm.step(data=dat_hist_Dist_mpo[dat_hist_Dist_mpo$abundance>0,],
                               gbm.x = c(1:2, 25, 27),
                               gbm.y = 'log_abundance',
                               family = "gaussian",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  presDx_mpo_nochl_S <- predict(brt_dist_P_mpo_nochl_S, dat_hist, n.trees=brt_dist_P_mpo_nochl_S$gbm.call$best.trees, type="response")
  abundDx_mpo_nochl_S <- exp(predict(brt_dist_N_mpo_nochl_S, dat_hist, n.trees=brt_dist_N_mpo_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_mpo_nochl_S <- presDx_mpo_nochl_S * abundDx_mpo_nochl_S
  
  presDx_mpo_nochl_S <- predict(brt_dist_P_mpo_nochl_S, dat_fcast, n.trees=brt_dist_P_mpo_nochl_S$gbm.call$best.trees, type="response")
  abundDx_mpo_nochl_S <- exp(predict(brt_dist_N_mpo_nochl_S, dat_fcast, n.trees=brt_dist_N_mpo_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_mpo_nochl_S <- presDx_mpo_nochl_S * abundDx_mpo_nochl_S
}

#Middle nearshore
if("mpn" %in% sampling){
  print("Fitting BRT-MPN_nochl_S")
  brt_dist_P_mpn_nochl_S <- gbm.step(data=dat_hist_Dist_mpn,
                               gbm.x = c(1:2, 25, 27),
                               gbm.y = 'pres',
                               family = "bernoulli",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_mpn_nochl_S <- gbm.step(data=dat_hist_Dist_mpn[dat_hist_Dist_mpn$abundance>0,],
                               gbm.x = c(1:2, 25, 27),
                               gbm.y = 'log_abundance',
                               family = "gaussian",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  presDx_mpn_nochl_S <- predict(brt_dist_P_mpn_nochl_S, dat_hist, n.trees=brt_dist_P_mpn_nochl_S$gbm.call$best.trees, type="response")
  abundDx_mpn_nochl_S <- exp(predict(brt_dist_N_mpn_nochl_S, dat_hist, n.trees=brt_dist_N_mpn_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_mpn_nochl_S <- presDx_mpn_nochl_S * abundDx_mpn_nochl_S
  
  presDx_mpn_nochl_S <- predict(brt_dist_P_mpn_nochl_S, dat_fcast, n.trees=brt_dist_P_mpn_nochl_S$gbm.call$best.trees, type="response")
  abundDx_mpn_nochl_S <- exp(predict(brt_dist_N_mpn_nochl_S, dat_fcast, n.trees=brt_dist_N_mpn_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_mpn_nochl_S <- presDx_mpn_nochl_S * abundDx_mpn_nochl_S
}

#Southern Offshore
if("spo" %in% sampling){
  print("Fitting BRT-SPO_nochl_S")
  brt_dist_P_nochl_Spo_nochl_S <- gbm.step(data=dat_hist_Dist_nochl_Spo,
                               gbm.x = c(1:2, 25, 27),
                               gbm.y = 'pres',
                               family = "bernoulli",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_nochl_Spo_nochl_S <- gbm.step(data=dat_hist_Dist_nochl_Spo[dat_hist_Dist_nochl_Spo$abundance>0,],
                               gbm.x = c(1:2, 25, 27),
                               gbm.y = 'log_abundance',
                               family = "gaussian",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  presDx_nochl_Spo_nochl_S <- predict(brt_dist_P_nochl_Spo_nochl_S, dat_hist, n.trees=brt_dist_P_nochl_Spo_nochl_S$gbm.call$best.trees, type="response")
  abundDx_nochl_Spo_nochl_S <- exp(predict(brt_dist_N_nochl_Spo_nochl_S, dat_hist, n.trees=brt_dist_N_nochl_Spo_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_nochl_Spo_nochl_S <- presDx_nochl_Spo_nochl_S * abundDx_nochl_Spo_nochl_S
  
  presDx_nochl_Spo_nochl_S <- predict(brt_dist_P_nochl_Spo_nochl_S, dat_fcast, n.trees=brt_dist_P_nochl_Spo_nochl_S$gbm.call$best.trees, type="response")
  abundDx_nochl_Spo_nochl_S <- exp(predict(brt_dist_N_nochl_Spo_nochl_S, dat_fcast, n.trees=brt_dist_N_nochl_Spo_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_nochl_Spo_nochl_S <- presDx_nochl_Spo_nochl_S * abundDx_nochl_Spo_nochl_S
}

#Southern Nearshore
if("spn" %in% sampling){
  print("Fitting BRT-SPN_nochl_S")
  brt_dist_P_nochl_Spn_nochl_S <- gbm.step(data=dat_hist_Dist_nochl_Spn,
                               gbm.x = c(1:2, 25, 27),
                               gbm.y = 'pres',
                               family = "bernoulli",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_nochl_Spn_nochl_S <- gbm.step(data=dat_hist_Dist_nochl_Spn[dat_hist_Dist_nochl_Spn$abundance>0,],
                               gbm.x = c(1:2, 25, 27),
                               gbm.y = 'log_abundance',
                               family = "gaussian",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  presDx_nochl_Spn_nochl_S <- predict(brt_dist_P_nochl_Spn_nochl_S, dat_hist, n.trees=brt_dist_P_nochl_Spn_nochl_S$gbm.call$best.trees, type="response")
  abundDx_nochl_Spn_nochl_S <- exp(predict(brt_dist_N_nochl_Spn_nochl_S, dat_hist, n.trees=brt_dist_N_nochl_Spn_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_nochl_Spn_nochl_S <- presDx_nochl_Spn_nochl_S * abundDx_nochl_Spn_nochl_S
  
  presDx_nochl_Spn_nochl_S <- predict(brt_dist_P_nochl_Spn_nochl_S, dat_fcast, n.trees=brt_dist_P_nochl_Spn_nochl_S$gbm.call$best.trees, type="response")
  abundDx_nochl_Spn_nochl_S <- exp(predict(brt_dist_N_nochl_Spn_nochl_S, dat_fcast, n.trees=brt_dist_N_nochl_Spn_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_nochl_Spn_nochl_S <- presDx_nochl_Spn_nochl_S * abundDx_nochl_Spn_nochl_S
}

#All Offshore
if("allo" %in% sampling){
  print("Fitting BRT-ALLO_nochl_S")
  brt_dist_P_allo_nochl_S <- gbm.step(data=dat_hist_Dist_allo,
                                gbm.x = c(1:2, 25, 27),
                                gbm.y = 'pres',
                                family = "bernoulli",
                                tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_allo_nochl_S <- gbm.step(data=dat_hist_Dist_allo[dat_hist_Dist_allo$abundance>0,],
                                gbm.x = c(1:2, 25, 27),
                                gbm.y = 'log_abundance',
                                family = "gaussian",
                                tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                plot.main=FALSE, verbose = FALSE)
  
  presDx_allo_nochl_S <- predict(brt_dist_P_allo_nochl_S, dat_hist, n.trees=brt_dist_P_allo_nochl_S$gbm.call$best.trees, type="response")
  abundDx_allo_nochl_S <- exp(predict(brt_dist_N_allo_nochl_S, dat_hist, n.trees=brt_dist_N_allo_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_allo_nochl_S <- presDx_allo_nochl_S * abundDx_allo_nochl_S
  
  presDx_allo_nochl_S <- predict(brt_dist_P_allo_nochl_S, dat_fcast, n.trees=brt_dist_P_allo_nochl_S$gbm.call$best.trees, type="response")
  abundDx_allo_nochl_S <- exp(predict(brt_dist_N_allo_nochl_S, dat_fcast, n.trees=brt_dist_N_allo_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_allo_nochl_S <- presDx_allo_nochl_S * abundDx_allo_nochl_S
}

#All Nearshore
if("alln" %in% sampling){
  print("Fitting BRT-ALLN_nochl_S")
  brt_dist_P_alln_nochl_S <- gbm.step(data=dat_hist_Dist_alln,
                                gbm.x = c(1:2, 25, 27),
                                gbm.y = 'pres',
                                family = "bernoulli",
                                tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_alln_nochl_S <- gbm.step(data=dat_hist_Dist_alln[dat_hist_Dist_alln$abundance>0,],
                                gbm.x = c(1:2, 25, 27),
                                gbm.y = 'log_abundance',
                                family = "gaussian",
                                tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                plot.main=FALSE, verbose = FALSE)
  
  presDx_alln_nochl_S <- predict(brt_dist_P_alln_nochl_S, dat_hist, n.trees=brt_dist_P_alln_nochl_S$gbm.call$best.trees, type="response")
  abundDx_alln_nochl_S <- exp(predict(brt_dist_N_alln_nochl_S, dat_hist, n.trees=brt_dist_N_alln_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_alln_nochl_S <- presDx_alln_nochl_S * abundDx_alln_nochl_S
  
  presDx_alln_nochl_S <- predict(brt_dist_P_alln_nochl_S, dat_fcast, n.trees=brt_dist_P_alln_nochl_S$gbm.call$best.trees, type="response")
  abundDx_alln_nochl_S <- exp(predict(brt_dist_N_alln_nochl_S, dat_fcast, n.trees=brt_dist_N_alln_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_alln_nochl_S <- presDx_alln_nochl_S * abundDx_alln_nochl_S
}
#
## Bycatch + Opt Target Sampling
#
if("BY" %in% sampling){
  print("Fitting BRT-BY_nochl_S")
  brt_B_P_nochl_S <- gbm.step(data=dat_hist_BY,
                        gbm.x = c(1:2, 25, 27),
                        gbm.y = 'pres',
                        family = "bernoulli",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  brt_B_N_nochl_S <- gbm.step(data=dat_hist_BY[dat_hist_BY$abundance>0,], 
                        gbm.x = c(1:2, 25, 27),
                        gbm.y = 'log_abundance',
                        family = "gaussian",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  presBx_nochl_S <- predict(brt_B_P_nochl_S, dat_hist, n.trees=brt_B_P_nochl_S$gbm.call$best.trees, type="response")
  abundBx_nochl_S <- exp(predict(brt_B_N_nochl_S, dat_hist, n.trees=brt_B_N_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_BY_nochl_S <- presBx_nochl_S * abundBx_nochl_S
  
  presBx_nochl_S <- predict(brt_B_P_nochl_S, dat_fcast, n.trees=brt_B_P_nochl_S$gbm.call$best.trees, type="response")
  abundBx_nochl_S <- exp(predict(brt_B_N_nochl_S, dat_fcast, n.trees=brt_B_N_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_BY_nochl_S<- presBx_nochl_S * abundBx_nochl_S
}

#
###Closed Areas + Opt Target Species
#
# Small Closed Area
if("CA_nochl_Sm" %in% sampling){
  print("Fitting BRT-CASM_nochl_S")
  brt_CA_P_nochl_Sm_nochl_S <- gbm.step(data=dat_hist_CA_nochl_Sm,
                            gbm.x = c(1:2, 25, 27),
                            gbm.y = 'pres',
                            family = "bernoulli",
                            tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                            plot.main=FALSE, verbose = FALSE)
  
  brt_CA_N_nochl_Sm_nochl_S <- gbm.step(data=dat_hist_CA_nochl_Sm[dat_hist_CA_nochl_Sm$abundance>0,], 
                            gbm.x = c(1:2, 25, 27),
                            gbm.y = 'log_abundance',
                            family = "gaussian",
                            tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                            plot.main=FALSE, verbose = FALSE)
  
  presCAx_nochl_Sm_nochl_S <- predict(brt_CA_P_nochl_Sm_nochl_S, dat_hist, n.trees=brt_CA_P_nochl_Sm_nochl_S$gbm.call$best.trees, type="response")
  abundCAx_nochl_Sm_nochl_S <- exp(predict(brt_CA_N_nochl_Sm_nochl_S, dat_hist, n.trees=brt_CA_N_nochl_Sm_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_CA_nochl_Sm_nochl_S <- presCAx_nochl_Sm_nochl_S * abundCAx_nochl_Sm_nochl_S
  
  presCAx_nochl_Sm_nochl_S <- predict(brt_CA_P_nochl_Sm_nochl_S, dat_fcast, n.trees=brt_CA_P_nochl_Sm_nochl_S$gbm.call$best.trees, type="response")
  abundCAx_nochl_Sm_nochl_S <- exp(predict(brt_CA_N_nochl_Sm_nochl_S, dat_fcast, n.trees=brt_CA_N_nochl_Sm_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_CA_nochl_Sm_nochl_S<- presCAx_nochl_Sm_nochl_S * abundCAx_nochl_Sm_nochl_S
}

# Medium Closed Area
if("CA_med" %in% sampling){
  print("Fitting BRT-CAMED_nochl_S")
  brt_CA_P_med_nochl_S <- gbm.step(data=dat_hist_CA_med,
                             gbm.x = c(1:2, 25, 27),
                             gbm.y = 'pres',
                             family = "bernoulli",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  brt_CA_N_med_nochl_S <- gbm.step(data=dat_hist_CA_med[dat_hist_CA_med$abundance>0,], 
                             gbm.x = c(1:2, 25, 27),
                             gbm.y = 'log_abundance',
                             family = "gaussian",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  presCAx_med_nochl_S <- predict(brt_CA_P_med_nochl_S, dat_hist, n.trees=brt_CA_P_med_nochl_S$gbm.call$best.trees, type="response")
  abundCAx_med_nochl_S <- exp(predict(brt_CA_N_med, dat_hist, n.trees=brt_CA_N_med_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_CA_med_nochl_S <- presCAx_med_nochl_S * abundCAx_med_nochl_S
  
  presCAx_med_nochl_S <- predict(brt_CA_P_med_nochl_S, dat_fcast, n.trees=brt_CA_P_med_nochl_S$gbm.call$best.trees, type="response")
  abundCAx_med_nochl_S <- exp(predict(brt_CA_N_med_nochl_S, dat_fcast, n.trees=brt_CA_N_med_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_CA_med_nochl_S<- presCAx_med_nochl_S * abundCAx_med_nochl_S
}

# Large Closed Area
if("CA_lar" %in% sampling){
  print("Fitting BRT-CALAR_nochl_S")
  brt_CA_P_lar_nochl_S <- gbm.step(data=dat_hist_CA_lar,
                             gbm.x = c(1:2, 25, 27),
                             gbm.y = 'pres',
                             family = "bernoulli",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  brt_CA_N_lar_nochl_S <- gbm.step(data=dat_hist_CA_lar[dat_hist_CA_lar$abundance>0,], 
                             gbm.x = c(1:2, 25, 27),
                             gbm.y = 'log_abundance',
                             family = "gaussian",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  presCAx_lar_nochl_S <- predict(brt_CA_P_lar_nochl_S, dat_hist, n.trees=brt_CA_P_lar_nochl_S$gbm.call$best.trees, type="response")
  abundCAx_lar_nochl_S <- exp(predict(brt_CA_N_lar_nochl_S, dat_hist, n.trees=brt_CA_N_lar_nochl_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_CA_lar_nochl_S <- presCAx_lar_nochl_S * abundCAx_lar_nochl_S
  
  presCAx_lar_nochl_S <- predict(brt_CA_P_lar_nochl_S, dat_fcast, n.trees=brt_CA_P_lar_nochl_S$gbm.call$best.trees, type="response")
  abundCAx_lar_nochl_S <- exp(predict(brt_CA_N_lar_nochl_S, dat_fcast, n.trees=brt_CA_N_lar_nochl_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_CA_lar_nochl_S<- presCAx_lar_nochl_S * abundCAx_lar_nochl_S
}
