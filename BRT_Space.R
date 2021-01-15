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
  print("Fitting BRT-Ran_S")
  brt_R_P_S <- gbm.step(data=dat_hist_random,
                      gbm.x = c(1:2, 25, 27, 28),
                      gbm.y = 'pres',
                      family = "bernoulli",
                      tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                      plot.main=FALSE, verbose = FALSE)
  
  brt_R_N_S <- gbm.step(data=dat_hist_random[dat_hist_random$abundance>0,], 
                      gbm.x = c(1:2, 25, 27, 28),
                      gbm.y = 'log_abundance',
                      family = "gaussian",
                      tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                      plot.main=FALSE, verbose = FALSE)
  
  presRx_S <- predict(brt_R_P_S, dat_hist, n.trees=brt_R_P_S$gbm.call$best.trees, type="response")
  abundRx_S <- exp(predict(brt_R_N_S, dat_hist, n.trees=brt_R_N_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Ran_S <- presRx_S * abundRx_S
  
  presRx_S <- predict(brt_R_P_S, dat_fcast, n.trees=brt_R_P_S$gbm.call$best.trees, type="response")
  abundRx_S <- exp(predict(brt_R_N_S, dat_fcast, n.trees=brt_R_N_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Ran_S <- presRx_S * abundRx_S
}

#
### Preferential (target) sampling
#
#Target sampling-0.5
if("tar_0.5" %in% sampling){
  print("Fitting BRT-Pref0.5_S")
  brt_T_P_1_S <- gbm.step(data=dat_hist_Tar_1,
                        gbm.x = c(1:2, 25, 27, 28),
                        gbm.y = 'pres',
                        family = "bernoulli",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_1_S <- gbm.step(data=dat_hist_Tar_1[dat_hist_Tar_1$abundance>0,],
                        gbm.x = c(1:2, 25, 27, 28),
                        gbm.y = 'log_abundance',
                        family = "gaussian",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  presTx_1_S <- predict(brt_T_P_1_S, dat_hist, n.trees=brt_T_P_1_S$gbm.call$best.trees, type="response")
  abundTx_1_S <- exp(predict(brt_T_N_1_S, dat_hist, n.trees=brt_T_N_1_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.5_S <- presTx_1_S * abundTx_1_S
  
  presTx_1_S <- predict(brt_T_P_1_S, dat_fcast, n.trees=brt_T_P_1_S$gbm.call$best.trees, type="response")
  abundTx_1_S <- exp(predict(brt_T_N_1_S, dat_fcast, n.trees=brt_T_N_1_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.5_S <- presTx_1_S * abundTx_1_S
}

#Target sampling-0.6
if("tar_0.6" %in% sampling){
  print("Fitting BRT-Pref0.6_S")
  brt_T_P_2_S <- gbm.step(data=dat_hist_Tar_2,
                        gbm.x = c(1:2, 25, 27, 28),
                        gbm.y = 'pres',
                        family = "bernoulli",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_2_S <- gbm.step(data=dat_hist_Tar_2[dat_hist_Tar_2$abundance>0,],
                        gbm.x = c(1:2, 25, 27, 28),
                        gbm.y = 'log_abundance',
                        family = "gaussian",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  presTx_2_S <- predict(brt_T_P_2_S, dat_hist, n.trees=brt_T_P_2_S$gbm.call$best.trees, type="response")
  abundTx_2_S <- exp(predict(brt_T_N_2_S, dat_hist, n.trees=brt_T_N_2_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.6_S <- presTx_2_S* abundTx_2_S
  
  presTx_2_S <- predict(brt_T_P_2_S, dat_fcast, n.trees=brt_T_P_2_S$gbm.call$best.trees, type="response")
  abundTx_2_S <- exp(predict(brt_T_N_2_S, dat_fcast, n.trees=brt_T_N_2_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.6_S <- presTx_2_S * abundTx_2_S
}

#Target sampling-0.7
if("tar_0.7" %in% sampling){
  print("Fitting BRT-Pref0.7_S")
  brt_T_P_3_S <- gbm.step(data=dat_hist_Tar_3,
                        gbm.x = c(1:2, 25, 27, 28),
                        gbm.y = 'pres',
                        family = "bernoulli",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_3_S <- gbm.step(data=dat_hist_Tar_3[dat_hist_Tar_3$abundance>0,],
                        gbm.x = c(1:2, 25, 27, 28),
                        gbm.y = 'log_abundance',
                        family = "gaussian",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  presTx_3_S <- predict(brt_T_P_3_S, dat_hist, n.trees=brt_T_P_3_S$gbm.call$best.trees, type="response")
  abundTx_3_S <- exp(predict(brt_T_N_3_S, dat_hist, n.trees=brt_T_N_3_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.7_S <- presTx_3_S * abundTx_3_S
  
  presTx_3_S <- predict(brt_T_P_3_S, dat_fcast, n.trees=brt_T_P_3_S$gbm.call$best.trees, type="response")
  abundTx_3_S <- exp(predict(brt_T_N_3_S, dat_fcast, n.trees=brt_T_N_3_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.7_S <- presTx_3_S * abundTx_3_S
}

#Target sampling-0.8
if("tar_0.8" %in% sampling){
  print("Fitting BRT-Pref0.8_S")
  brt_T_P_4_S <- gbm.step(data=dat_hist_Tar_4,
                        gbm.x = c(1:2, 25, 27, 28),
                        gbm.y = 'pres',
                        family = "bernoulli",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_4_S <- gbm.step(data=dat_hist_Tar_4[dat_hist_Tar_4$abundance>0,],
                        gbm.x = c(1:2, 25, 27, 28),
                        gbm.y = 'log_abundance',
                        family = "gaussian",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  presTx_4_S <- predict(brt_T_P_4_S, dat_hist, n.trees=brt_T_P_4_S$gbm.call$best.trees, type="response")
  abundTx_4_S <- exp(predict(brt_T_N_4_S, dat_hist, n.trees=brt_T_N_4_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.8_S <- presTx_4_S * abundTx_4_S
  
  presTx_4_S <- predict(brt_T_P_4_S, dat_fcast, n.trees=brt_T_P_4_S$gbm.call$best.trees, type="response")
  abundTx_4_S <- exp(predict(brt_T_N_4_S, dat_fcast, n.trees=brt_T_N_4_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.8_S <- presTx_4_S * abundTx_4_S
}

#Target sampling-0.9
if("tar_0.9" %in% sampling){
  print("Fitting BRT-Pref0.9_S")
  brt_T_P_5_S <- gbm.step(data=dat_hist_Tar_5,
                        gbm.x = c(1:2, 25, 27, 28),
                        gbm.y = 'pres',
                        family = "bernoulli",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_5_S <- gbm.step(data=dat_hist_Tar_5[dat_hist_Tar_5$abundance>0,],
                        gbm.x = c(1:2, 25, 27, 28),
                        gbm.y = 'log_abundance',
                        family = "gaussian",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  presTx_5_S <- predict(brt_T_P_5_S, dat_hist, n.trees=brt_T_P_5_S$gbm.call$best.trees, type="response")
  abundTx_5_S <- exp(predict(brt_T_N_5_S, dat_hist, n.trees=brt_T_N_5_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.9_S <- presTx_5_S * abundTx_5_S
  
  presTx_5_S <- predict(brt_T_P_5_S, dat_fcast, n.trees=brt_T_P_5_S$gbm.call$best.trees, type="response")
  abundTx_5_S <- exp(predict(brt_T_N_5_S, dat_fcast, n.trees=brt_T_N_5_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.9_S <- presTx_5_S * abundTx_5_S
}

#
### Distance from Port Sampling
#
#Northern Offshore
if("npo" %in% sampling){
  print("Fitting BRT-NPO_S")
  brt_dist_P_npo_S <- gbm.step(data=dat_hist_Dist_npo,
                             gbm.x = c(1:2, 25, 27, 28),
                             gbm.y = 'pres',
                             family = "bernoulli",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_npo_S <- gbm.step(data=dat_hist_Dist_npo[dat_hist_Dist_npo$abundance>0,],
                             gbm.x = c(1:2, 25, 27, 28),
                             gbm.y = 'log_abundance',
                             family = "gaussian",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  presDx_npo_S <- predict(brt_dist_P_npo_S, dat_hist, n.trees=brt_dist_P_npo_S$gbm.call$best.trees, type="response")
  abundDx_npo_S <- exp(predict(brt_dist_N_npo_S, dat_hist, n.trees=brt_dist_N_npo_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_npo_S <- presDx_npo_S * abundDx_npo_S
  
  presDx_npo_S <- predict(brt_dist_P_npo_S, dat_fcast, n.trees=brt_dist_P_npo_S$gbm.call$best.trees, type="response")
  abundDx_npo_S <- exp(predict(brt_dist_N_npo_S, dat_fcast, n.trees=brt_dist_N_npo_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_npo_S <- presDx_npo_S * abundDx_npo_S
}

#Northern Nearshore
if("npn" %in% sampling){
  print("Fitting BRT-NPN_S")
  brt_dist_P_npn_S <- gbm.step(data=dat_hist_Dist_npn,
                             gbm.x = c(1:2, 25, 27, 28),
                             gbm.y = 'pres',
                             family = "bernoulli",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_npn_S <- gbm.step(data=dat_hist_Dist_npn[dat_hist_Dist_npn$abundance>0,],
                             gbm.x = c(1:2, 25, 27, 28),
                             gbm.y = 'log_abundance',
                             family = "gaussian",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  presDx_npn_S <- predict(brt_dist_P_npn_S, dat_hist, n.trees=brt_dist_P_npn_S$gbm.call$best.trees, type="response")
  abundDx_npn_S <- exp(predict(brt_dist_N_npn_S, dat_hist, n.trees=brt_dist_N_npn_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_npn_S <- presDx_npn_S * abundDx_npn_S
  
  presDx_npn_S <- predict(brt_dist_P_npn_S, dat_fcast, n.trees=brt_dist_P_npn_S$gbm.call$best.trees, type="response")
  abundDx_npn_S <- exp(predict(brt_dist_N_npn_S, dat_fcast, n.trees=brt_dist_N_npn_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_npn_S <- presDx_npn_S * abundDx_npn_S
}  

#Middle Offshore
if("mpo" %in% sampling){
  print("Fitting BRT-MPO_S")
  brt_dist_P_mpo_S <- gbm.step(data=dat_hist_Dist_mpo,
                             gbm.x = c(1:2, 25, 27, 28),
                             gbm.y = 'pres',
                             family = "bernoulli",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_mpo_S <- gbm.step(data=dat_hist_Dist_mpo[dat_hist_Dist_mpo$abundance>0,],
                             gbm.x = c(1:2, 25, 27, 28),
                             gbm.y = 'log_abundance',
                             family = "gaussian",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  presDx_mpo_S <- predict(brt_dist_P_mpo_S, dat_hist, n.trees=brt_dist_P_mpo_S$gbm.call$best.trees, type="response")
  abundDx_mpo_S <- exp(predict(brt_dist_N_mpo_S, dat_hist, n.trees=brt_dist_N_mpo_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_mpo_S <- presDx_mpo_S * abundDx_mpo_S
  
  presDx_mpo_S <- predict(brt_dist_P_mpo_S, dat_fcast, n.trees=brt_dist_P_mpo_S$gbm.call$best.trees, type="response")
  abundDx_mpo_S <- exp(predict(brt_dist_N_mpo_S, dat_fcast, n.trees=brt_dist_N_mpo_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_mpo_S <- presDx_mpo_S * abundDx_mpo_S
}

#Middle nearshore
if("mpn" %in% sampling){
  print("Fitting BRT-MPN_S")
  brt_dist_P_mpn_S <- gbm.step(data=dat_hist_Dist_mpn,
                             gbm.x = c(1:2, 25, 27, 28),
                             gbm.y = 'pres',
                             family = "bernoulli",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_mpn_S <- gbm.step(data=dat_hist_Dist_mpn[dat_hist_Dist_mpn$abundance>0,],
                             gbm.x = c(1:2, 25, 27, 28),
                             gbm.y = 'log_abundance',
                             family = "gaussian",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  presDx_mpn_S <- predict(brt_dist_P_mpn_S, dat_hist, n.trees=brt_dist_P_mpn_S$gbm.call$best.trees, type="response")
  abundDx_mpn_S <- exp(predict(brt_dist_N_mpn_S, dat_hist, n.trees=brt_dist_N_mpn_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_mpn_S <- presDx_mpn_S * abundDx_mpn_S
  
  presDx_mpn_S <- predict(brt_dist_P_mpn_S, dat_fcast, n.trees=brt_dist_P_mpn_S$gbm.call$best.trees, type="response")
  abundDx_mpn_S <- exp(predict(brt_dist_N_mpn_S, dat_fcast, n.trees=brt_dist_N_mpn_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_mpn_S <- presDx_mpn_S * abundDx_mpn_S
}

#Southern Offshore
if("spo" %in% sampling){
  print("Fitting BRT-SPO_S")
  brt_dist_P_spo_S <- gbm.step(data=dat_hist_Dist_spo,
                             gbm.x = c(1:2, 25, 27, 28),
                             gbm.y = 'pres',
                             family = "bernoulli",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_spo_S <- gbm.step(data=dat_hist_Dist_spo[dat_hist_Dist_spo$abundance>0,],
                             gbm.x = c(1:2, 25, 27, 28),
                             gbm.y = 'log_abundance',
                             family = "gaussian",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  presDx_spo_S <- predict(brt_dist_P_spo_S, dat_hist, n.trees=brt_dist_P_spo_S$gbm.call$best.trees, type="response")
  abundDx_spo_S <- exp(predict(brt_dist_N_spo_S, dat_hist, n.trees=brt_dist_N_spo_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_spo_S <- presDx_spo_S * abundDx_spo_S
  
  presDx_spo_S <- predict(brt_dist_P_spo_S, dat_fcast, n.trees=brt_dist_P_spo_S$gbm.call$best.trees, type="response")
  abundDx_spo_S <- exp(predict(brt_dist_N_spo_S, dat_fcast, n.trees=brt_dist_N_spo_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_spo_S <- presDx_spo_S * abundDx_spo_S
}

#Southern Nearshore
if("spn" %in% sampling){
  print("Fitting BRT-SPN_S")
  brt_dist_P_spn_S <- gbm.step(data=dat_hist_Dist_spn,
                             gbm.x = c(1:2, 25, 27, 28),
                             gbm.y = 'pres',
                             family = "bernoulli",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_spn_S <- gbm.step(data=dat_hist_Dist_spn[dat_hist_Dist_spn$abundance>0,],
                             gbm.x = c(1:2, 25, 27, 28),
                             gbm.y = 'log_abundance',
                             family = "gaussian",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  presDx_spn_S <- predict(brt_dist_P_spn_S, dat_hist, n.trees=brt_dist_P_spn_S$gbm.call$best.trees, type="response")
  abundDx_spn_S <- exp(predict(brt_dist_N_spn_S, dat_hist, n.trees=brt_dist_N_spn_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_spn_S <- presDx_spn_S * abundDx_spn_S
  
  presDx_spn_S <- predict(brt_dist_P_spn_S, dat_fcast, n.trees=brt_dist_P_spn_S$gbm.call$best.trees, type="response")
  abundDx_spn_S <- exp(predict(brt_dist_N_spn_S, dat_fcast, n.trees=brt_dist_N_spn_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_spn_S <- presDx_spn_S * abundDx_spn_S
}

#All Offshore
if("allo" %in% sampling){
  print("Fitting BRT-ALLO_S")
  brt_dist_P_allo_S <- gbm.step(data=dat_hist_Dist_allo,
                              gbm.x = c(1:2, 25, 27, 28),
                              gbm.y = 'pres',
                              family = "bernoulli",
                              tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                              plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_allo_S <- gbm.step(data=dat_hist_Dist_allo[dat_hist_Dist_allo$abundance>0,],
                              gbm.x = c(1:2, 25, 27, 28),
                              gbm.y = 'log_abundance',
                              family = "gaussian",
                              tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                              plot.main=FALSE, verbose = FALSE)
  
  presDx_allo_S <- predict(brt_dist_P_allo_S, dat_hist, n.trees=brt_dist_P_allo_S$gbm.call$best.trees, type="response")
  abundDx_allo_S <- exp(predict(brt_dist_N_allo_S, dat_hist, n.trees=brt_dist_N_allo_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_allo_S <- presDx_allo_S * abundDx_allo_S
  
  presDx_allo_S <- predict(brt_dist_P_allo_S, dat_fcast, n.trees=brt_dist_P_allo_S$gbm.call$best.trees, type="response")
  abundDx_allo_S <- exp(predict(brt_dist_N_allo_S, dat_fcast, n.trees=brt_dist_N_allo_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_allo_S <- presDx_allo_S * abundDx_allo_S
}

#All Nearshore
if("alln" %in% sampling){
  print("Fitting BRT-ALLN_S")
  brt_dist_P_alln_S <- gbm.step(data=dat_hist_Dist_alln,
                              gbm.x = c(1:2, 25, 27, 28),
                              gbm.y = 'pres',
                              family = "bernoulli",
                              tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                              plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_alln_S <- gbm.step(data=dat_hist_Dist_alln[dat_hist_Dist_alln$abundance>0,],
                              gbm.x = c(1:2, 25, 27, 28),
                              gbm.y = 'log_abundance',
                              family = "gaussian",
                              tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                              plot.main=FALSE, verbose = FALSE)
  
  presDx_alln_S <- predict(brt_dist_P_alln_S, dat_hist, n.trees=brt_dist_P_alln_S$gbm.call$best.trees, type="response")
  abundDx_alln_S <- exp(predict(brt_dist_N_alln_S, dat_hist, n.trees=brt_dist_N_alln_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_alln_S <- presDx_alln_S * abundDx_alln_S
  
  presDx_alln_S <- predict(brt_dist_P_alln_S, dat_fcast, n.trees=brt_dist_P_alln_S$gbm.call$best.trees, type="response")
  abundDx_alln_S <- exp(predict(brt_dist_N_alln_S, dat_fcast, n.trees=brt_dist_N_alln_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_alln_S <- presDx_alln_S * abundDx_alln_S
}
#
## Bycatch + Opt Target Sampling
#
if("BY" %in% sampling){
  print("Fitting BRT-BY_S")
  brt_B_P_S <- gbm.step(data=dat_hist_BY,
                      gbm.x = c(1:2, 25, 27, 28),
                      gbm.y = 'pres',
                      family = "bernoulli",
                      tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                      plot.main=FALSE, verbose = FALSE)
  
  brt_B_N_S <- gbm.step(data=dat_hist_BY[dat_hist_BY$abundance>0,], 
                      gbm.x = c(1:2, 25, 27, 28),
                      gbm.y = 'log_abundance',
                      family = "gaussian",
                      tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                      plot.main=FALSE, verbose = FALSE)
  
  presBx_S <- predict(brt_B_P_S, dat_hist, n.trees=brt_B_P_S$gbm.call$best.trees, type="response")
  abundBx_S <- exp(predict(brt_B_N_S, dat_hist, n.trees=brt_B_N_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_BY_S <- presBx_S * abundBx_S
  
  presBx_S <- predict(brt_B_P_S, dat_fcast, n.trees=brt_B_P_S$gbm.call$best.trees, type="response")
  abundBx_S <- exp(predict(brt_B_N_S, dat_fcast, n.trees=brt_B_N_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_BY_S<- presBx_S * abundBx_S
}

#
###Closed Areas + Opt Target Species
#
# Small Closed Area
if("CA_sm" %in% sampling){
  print("Fitting BRT-CASM_S")
  brt_CA_P_sm_S <- gbm.step(data=dat_hist_CA_sm,
                          gbm.x = c(1:2, 25, 27, 28),
                          gbm.y = 'pres',
                          family = "bernoulli",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  brt_CA_N_sm_S <- gbm.step(data=dat_hist_CA_sm[dat_hist_CA_sm$abundance>0,], 
                          gbm.x = c(1:2, 25, 27, 28),
                          gbm.y = 'log_abundance',
                          family = "gaussian",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  presCAx_sm_S <- predict(brt_CA_P_sm_S, dat_hist, n.trees=brt_CA_P_sm_S$gbm.call$best.trees, type="response")
  abundCAx_sm_S <- exp(predict(brt_CA_N_sm_S, dat_hist, n.trees=brt_CA_N_sm_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_CA_sm_S <- presCAx_sm_S * abundCAx_sm_S
  
  presCAx_sm_S <- predict(brt_CA_P_sm_S, dat_fcast, n.trees=brt_CA_P_sm_S$gbm.call$best.trees, type="response")
  abundCAx_sm_S <- exp(predict(brt_CA_N_sm_S, dat_fcast, n.trees=brt_CA_N_sm_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_CA_sm_S<- presCAx_sm_S * abundCAx_sm_S
}

# Medium Closed Area
if("CA_med" %in% sampling){
  print("Fitting BRT-CAMED_S")
  brt_CA_P_med_S <- gbm.step(data=dat_hist_CA_med,
                           gbm.x = c(1:2, 25, 27, 28),
                           gbm.y = 'pres',
                           family = "bernoulli",
                           tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                           plot.main=FALSE, verbose = FALSE)
  
  brt_CA_N_med_S <- gbm.step(data=dat_hist_CA_med[dat_hist_CA_med$abundance>0,], 
                           gbm.x = c(1:2, 25, 27, 28),
                           gbm.y = 'log_abundance',
                           family = "gaussian",
                           tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                           plot.main=FALSE, verbose = FALSE)
  
  presCAx_med_S <- predict(brt_CA_P_med_S, dat_hist, n.trees=brt_CA_P_med_S$gbm.call$best.trees, type="response")
  abundCAx_med_S <- exp(predict(brt_CA_N_med, dat_hist, n.trees=brt_CA_N_med_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_CA_med_S <- presCAx_med_S * abundCAx_med_S
  
  presCAx_med_S <- predict(brt_CA_P_med_S, dat_fcast, n.trees=brt_CA_P_med_S$gbm.call$best.trees, type="response")
  abundCAx_med_S <- exp(predict(brt_CA_N_med_S, dat_fcast, n.trees=brt_CA_N_med_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_CA_med_S<- presCAx_med_S * abundCAx_med_S
}

# Large Closed Area
if("CA_lar" %in% sampling){
  print("Fitting BRT-CALAR_S")
  brt_CA_P_lar_S <- gbm.step(data=dat_hist_CA_lar,
                           gbm.x = c(1:2, 25, 27, 28),
                           gbm.y = 'pres',
                           family = "bernoulli",
                           tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                           plot.main=FALSE, verbose = FALSE)
  
  brt_CA_N_lar_S <- gbm.step(data=dat_hist_CA_lar[dat_hist_CA_lar$abundance>0,], 
                           gbm.x = c(1:2, 25, 27, 28),
                           gbm.y = 'log_abundance',
                           family = "gaussian",
                           tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                           plot.main=FALSE, verbose = FALSE)
  
  presCAx_lar_S <- predict(brt_CA_P_lar_S, dat_hist, n.trees=brt_CA_P_lar_S$gbm.call$best.trees, type="response")
  abundCAx_lar_S <- exp(predict(brt_CA_N_lar_S, dat_hist, n.trees=brt_CA_N_lar_S$gbm.call$best.trees, type="response"))
  dat_hist$brt_CA_lar_S <- presCAx_lar_S * abundCAx_lar_S
  
  presCAx_lar_S <- predict(brt_CA_P_lar_S, dat_fcast, n.trees=brt_CA_P_lar_S$gbm.call$best.trees, type="response")
  abundCAx_lar_S <- exp(predict(brt_CA_N_lar_S, dat_fcast, n.trees=brt_CA_N_lar_S$gbm.call$best.trees, type="response"))
  dat_fcast$brt_CA_lar_S<- presCAx_lar_S * abundCAx_lar_S
}
