### Fishery-dependent SDM (L^3) Estimation Model Code
## Function to run BRTs on data from OM, with all Covariates (temp, chla, mld)
## Function returns only the fitted and predicted values
## modified my M.Karp 1/04/21

#############################
#   Full Models- with Chla  #
#############################
#
###Random sampling
#
    if("ran" %in% sampling) {
    print("Fitting BRT-Ran")
    brt_R_P <- gbm.step(data=dat_hist_random,
                        gbm.x = c(25, 27, 28),
                        gbm.y = 'pres',
                        family = "bernoulli",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)

    brt_R_N <- gbm.step(data=dat_hist_random[dat_hist_random$abundance>0,], 
                        gbm.x = c(25, 27, 28),
                        gbm.y = 'log_abundance',
                        family = "gaussian",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)

    presRx <- predict(brt_R_P, dat_hist, n.trees=brt_R_P$gbm.call$best.trees, type="response")
    abundRx <- exp(predict(brt_R_N, dat_hist, n.trees=brt_R_N$gbm.call$best.trees, type="response"))
    dat_hist$brt_Ran <- presRx * abundRx
    
    presRx <- predict(brt_R_P, dat_fcast, n.trees=brt_R_P$gbm.call$best.trees, type="response")
    abundRx <- exp(predict(brt_R_N, dat_fcast, n.trees=brt_R_N$gbm.call$best.trees, type="response"))
    dat_fcast$brt_Ran <- presRx * abundRx
    }

#
### Preferential (target) sampling
#
    #Target sampling-0.5
    if("tar_0.5" %in% sampling){
      print("Fitting BRT-Pref0.5")
      brt_T_P_1 <- gbm.step(data=dat_hist_Tar_1,
                          gbm.x = c(25,27,28),
                          gbm.y = 'pres',
                          family = "bernoulli",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
      
      brt_T_N_1 <- gbm.step(data=dat_hist_Tar_1[dat_hist_Tar_1$abundance>0,],
                          gbm.x = c(25,27,28),
                          gbm.y = 'log_abundance',
                          family = "gaussian",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
      
      presTx_1 <- predict(brt_T_P_1, dat_hist, n.trees=brt_T_P_1$gbm.call$best.trees, type="response")
      abundTx_1 <- exp(predict(brt_T_N_1, dat_hist, n.trees=brt_T_N_1$gbm.call$best.trees, type="response"))
      dat_hist$brt_Tar_0.5 <- presTx_1 * abundTx_1
      
      presTx_1 <- predict(brt_T_P_1, dat_fcast, n.trees=brt_T_P_1$gbm.call$best.trees, type="response")
      abundTx_1 <- exp(predict(brt_T_N_1, dat_fcast, n.trees=brt_T_N_1$gbm.call$best.trees, type="response"))
      dat_fcast$brt_Tar_0.5 <- presTx_1 * abundTx_1
    }

    #Target sampling-0.6
    if("tar_0.6" %in% sampling){
      print("Fitting BRT-Pref0.6")
      brt_T_P_2 <- gbm.step(data=dat_hist_Tar_2,
                            gbm.x = c(25,27,28),
                            gbm.y = 'pres',
                            family = "bernoulli",
                            tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                            plot.main=FALSE, verbose = FALSE)
      
      brt_T_N_2 <- gbm.step(data=dat_hist_Tar_2[dat_hist_Tar_2$abundance>0,],
                            gbm.x = c(25,27,28),
                            gbm.y = 'log_abundance',
                            family = "gaussian",
                            tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                            plot.main=FALSE, verbose = FALSE)
      
      presTx_2 <- predict(brt_T_P_2, dat_hist, n.trees=brt_T_P_2$gbm.call$best.trees, type="response")
      abundTx_2 <- exp(predict(brt_T_N_2, dat_hist, n.trees=brt_T_N_2$gbm.call$best.trees, type="response"))
      dat_hist$brt_Tar_0.6 <- presTx_2 * abundTx_2
      
      presTx_2 <- predict(brt_T_P_2, dat_fcast, n.trees=brt_T_P_2$gbm.call$best.trees, type="response")
      abundTx_2 <- exp(predict(brt_T_N_2, dat_fcast, n.trees=brt_T_N_2$gbm.call$best.trees, type="response"))
      dat_fcast$brt_Tar_0.6 <- presTx_2 * abundTx_2
    }

    #Target sampling-0.7
    if("tar_0.7" %in% sampling){
      print("Fitting BRT-Pref0.7")
      brt_T_P_3 <- gbm.step(data=dat_hist_Tar_3,
                            gbm.x = c(25,27,28),
                            gbm.y = 'pres',
                            family = "bernoulli",
                            tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                            plot.main=FALSE, verbose = FALSE)
      
      brt_T_N_3 <- gbm.step(data=dat_hist_Tar_3[dat_hist_Tar_3$abundance>0,],
                            gbm.x = c(25,27,28),
                            gbm.y = 'log_abundance',
                            family = "gaussian",
                            tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                            plot.main=FALSE, verbose = FALSE)
      
      presTx_3 <- predict(brt_T_P_3, dat_hist, n.trees=brt_T_P_3$gbm.call$best.trees, type="response")
      abundTx_3 <- exp(predict(brt_T_N_3, dat_hist, n.trees=brt_T_N_3$gbm.call$best.trees, type="response"))
      dat_hist$brt_Tar_0.7 <- presTx_3 * abundTx_3
      
      presTx_3 <- predict(brt_T_P_3, dat_fcast, n.trees=brt_T_P_3$gbm.call$best.trees, type="response")
      abundTx_3 <- exp(predict(brt_T_N_3, dat_fcast, n.trees=brt_T_N_3$gbm.call$best.trees, type="response"))
      dat_fcast$brt_Tar_0.7 <- presTx_3 * abundTx_3
    }

    #Target sampling-0.8
    if("tar_0.8" %in% sampling){
      print("Fitting BRT-Pref0.8")
      brt_T_P_4 <- gbm.step(data=dat_hist_Tar_4,
                            gbm.x = c(25,27,28),
                            gbm.y = 'pres',
                            family = "bernoulli",
                            tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                            plot.main=FALSE, verbose = FALSE)
      
      brt_T_N_4 <- gbm.step(data=dat_hist_Tar_4[dat_hist_Tar_4$abundance>0,],
                            gbm.x = c(25,27,28),
                            gbm.y = 'log_abundance',
                            family = "gaussian",
                            tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                            plot.main=FALSE, verbose = FALSE)
      
      presTx_4 <- predict(brt_T_P_4, dat_hist, n.trees=brt_T_P_4$gbm.call$best.trees, type="response")
      abundTx_4 <- exp(predict(brt_T_N_4, dat_hist, n.trees=brt_T_N_4$gbm.call$best.trees, type="response"))
      dat_hist$brt_Tar_0.8 <- presTx_4 * abundTx_4
      
      presTx_4 <- predict(brt_T_P_4, dat_fcast, n.trees=brt_T_P_4$gbm.call$best.trees, type="response")
      abundTx_4 <- exp(predict(brt_T_N_4, dat_fcast, n.trees=brt_T_N_4$gbm.call$best.trees, type="response"))
      dat_fcast$brt_Tar_0.8 <- presTx_4 * abundTx_4
    }

    #Target sampling-0.9
    if("tar_0.9" %in% sampling){
      print("Fitting BRT-Pref0.9")
      brt_T_P_5 <- gbm.step(data=dat_hist_Tar_5,
                            gbm.x = c(25,27,28),
                            gbm.y = 'pres',
                            family = "bernoulli",
                            tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                            plot.main=FALSE, verbose = FALSE)
      
      brt_T_N_5 <- gbm.step(data=dat_hist_Tar_5[dat_hist_Tar_5$abundance>0,],
                            gbm.x = c(25,27,28),
                            gbm.y = 'log_abundance',
                            family = "gaussian",
                            tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                            plot.main=FALSE, verbose = FALSE)
      
      presTx_5 <- predict(brt_T_P_5, dat_hist, n.trees=brt_T_P_5$gbm.call$best.trees, type="response")
      abundTx_5 <- exp(predict(brt_T_N_5, dat_hist, n.trees=brt_T_N_5$gbm.call$best.trees, type="response"))
      dat_hist$brt_Tar_0.9 <- presTx_5 * abundTx_5
      
      presTx_5 <- predict(brt_T_P_5, dat_fcast, n.trees=brt_T_P_5$gbm.call$best.trees, type="response")
      abundTx_5 <- exp(predict(brt_T_N_5, dat_fcast, n.trees=brt_T_N_5$gbm.call$best.trees, type="response"))
      dat_fcast$brt_Tar_0.9 <- presTx_5 * abundTx_5
    }

#
### Distance from Port Sampling
#
    #Northern Offshore
    if("npo" %in% sampling){
      print("Fitting BRT-NPO")
      brt_dist_P_npo <- gbm.step(data=dat_hist_Dist_npo,
                            gbm.x = c(25,27,28),
                            gbm.y = 'pres',
                            family = "bernoulli",
                            tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                            plot.main=FALSE, verbose = FALSE)
      
      brt_dist_N_npo <- gbm.step(data=dat_hist_Dist_npo[dat_hist_Dist_npo$abundance>0,],
                            gbm.x = c(25,27,28),
                            gbm.y = 'log_abundance',
                            family = "gaussian",
                            tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                            plot.main=FALSE, verbose = FALSE)
      
      presDx_npo <- predict(brt_dist_P_npo, dat_hist, n.trees=brt_dist_P_npo$gbm.call$best.trees, type="response")
      abundDx_npo <- exp(predict(brt_dist_N_npo, dat_hist, n.trees=brt_dist_N_npo$gbm.call$best.trees, type="response"))
      dat_hist$brt_Dist_npo <- presDx_npo * abundDx_npo
      
      presDx_npo <- predict(brt_dist_P_npo, dat_fcast, n.trees=brt_dist_P_npo$gbm.call$best.trees, type="response")
      abundDx_npo <- exp(predict(brt_dist_N_npo, dat_fcast, n.trees=brt_dist_N_npo$gbm.call$best.trees, type="response"))
      dat_fcast$brt_Dist_npo <- presDx_npo * abundDx_npo
    }
  
    #Northern Nearshore
    if("npn" %in% sampling){
      print("Fitting BRT-NPN")
      brt_dist_P_npn <- gbm.step(data=dat_hist_Dist_npn,
                                 gbm.x = c(25,27,28),
                                 gbm.y = 'pres',
                                 family = "bernoulli",
                                 tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                 plot.main=FALSE, verbose = FALSE)
      
      brt_dist_N_npn <- gbm.step(data=dat_hist_Dist_npn[dat_hist_Dist_npn$abundance>0,],
                                 gbm.x = c(25,27,28),
                                 gbm.y = 'log_abundance',
                                 family = "gaussian",
                                 tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                 plot.main=FALSE, verbose = FALSE)
      
      presDx_npn <- predict(brt_dist_P_npn, dat_hist, n.trees=brt_dist_P_npn$gbm.call$best.trees, type="response")
      abundDx_npn <- exp(predict(brt_dist_N_npn, dat_hist, n.trees=brt_dist_N_npn$gbm.call$best.trees, type="response"))
      dat_hist$brt_Dist_npn <- presDx_npn * abundDx_npn
      
      presDx_npn <- predict(brt_dist_P_npn, dat_fcast, n.trees=brt_dist_P_npn$gbm.call$best.trees, type="response")
      abundDx_npn <- exp(predict(brt_dist_N_npn, dat_fcast, n.trees=brt_dist_N_npn$gbm.call$best.trees, type="response"))
      dat_fcast$brt_Dist_npn <- presDx_npn * abundDx_npn
    }  

    #Middle Offshore
    if("mpo" %in% sampling){
      print("Fitting BRT-MPO")
      brt_dist_P_mpo <- gbm.step(data=dat_hist_Dist_mpo,
                                 gbm.x = c(25,27,28),
                                 gbm.y = 'pres',
                                 family = "bernoulli",
                                 tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                 plot.main=FALSE, verbose = FALSE)
      
      brt_dist_N_mpo <- gbm.step(data=dat_hist_Dist_mpo[dat_hist_Dist_mpo$abundance>0,],
                                 gbm.x = c(25,27,28),
                                 gbm.y = 'log_abundance',
                                 family = "gaussian",
                                 tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                 plot.main=FALSE, verbose = FALSE)
      
      presDx_mpo <- predict(brt_dist_P_mpo, dat_hist, n.trees=brt_dist_P_mpo$gbm.call$best.trees, type="response")
      abundDx_mpo <- exp(predict(brt_dist_N_mpo, dat_hist, n.trees=brt_dist_N_mpo$gbm.call$best.trees, type="response"))
      dat_hist$brt_Dist_mpo <- presDx_mpo * abundDx_mpo
      
      presDx_mpo <- predict(brt_dist_P_mpo, dat_fcast, n.trees=brt_dist_P_mpo$gbm.call$best.trees, type="response")
      abundDx_mpo <- exp(predict(brt_dist_N_mpo, dat_fcast, n.trees=brt_dist_N_mpo$gbm.call$best.trees, type="response"))
      dat_fcast$brt_Dist_mpo <- presDx_mpo * abundDx_mpo
    }

    #Middle nearshore
    if("mpn" %in% sampling){
      print("Fitting BRT-MPN")
      brt_dist_P_mpn <- gbm.step(data=dat_hist_Dist_mpn,
                                 gbm.x = c(25,27,28),
                                 gbm.y = 'pres',
                                 family = "bernoulli",
                                 tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                 plot.main=FALSE, verbose = FALSE)
      
      brt_dist_N_mpn <- gbm.step(data=dat_hist_Dist_mpn[dat_hist_Dist_mpn$abundance>0,],
                                 gbm.x = c(25,27,28),
                                 gbm.y = 'log_abundance',
                                 family = "gaussian",
                                 tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                 plot.main=FALSE, verbose = FALSE)
      
      presDx_mpn <- predict(brt_dist_P_mpn, dat_hist, n.trees=brt_dist_P_mpn$gbm.call$best.trees, type="response")
      abundDx_mpn <- exp(predict(brt_dist_N_mpn, dat_hist, n.trees=brt_dist_N_mpn$gbm.call$best.trees, type="response"))
      dat_hist$brt_Dist_mpn <- presDx_mpn * abundDx_mpn
      
      presDx_mpn <- predict(brt_dist_P_mpn, dat_fcast, n.trees=brt_dist_P_mpn$gbm.call$best.trees, type="response")
      abundDx_mpn <- exp(predict(brt_dist_N_mpn, dat_fcast, n.trees=brt_dist_N_mpn$gbm.call$best.trees, type="response"))
      dat_fcast$brt_Dist_mpn <- presDx_mpn * abundDx_mpn
    }

    #Southern Offshore
    if("spo" %in% sampling){
      print("Fitting BRT-SPO")
      brt_dist_P_spo <- gbm.step(data=dat_hist_Dist_spo,
                                 gbm.x = c(25,27,28),
                                 gbm.y = 'pres',
                                 family = "bernoulli",
                                 tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                 plot.main=FALSE, verbose = FALSE)
      
      brt_dist_N_spo <- gbm.step(data=dat_hist_Dist_spo[dat_hist_Dist_spo$abundance>0,],
                                 gbm.x = c(25,27,28),
                                 gbm.y = 'log_abundance',
                                 family = "gaussian",
                                 tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                 plot.main=FALSE, verbose = FALSE)
      
      presDx_spo <- predict(brt_dist_P_spo, dat_hist, n.trees=brt_dist_P_spo$gbm.call$best.trees, type="response")
      abundDx_spo <- exp(predict(brt_dist_N_spo, dat_hist, n.trees=brt_dist_N_spo$gbm.call$best.trees, type="response"))
      dat_hist$brt_Dist_spo <- presDx_spo * abundDx_spo
      
      presDx_spo <- predict(brt_dist_P_spo, dat_fcast, n.trees=brt_dist_P_spo$gbm.call$best.trees, type="response")
      abundDx_spo <- exp(predict(brt_dist_N_spo, dat_fcast, n.trees=brt_dist_N_spo$gbm.call$best.trees, type="response"))
      dat_fcast$brt_Dist_spo <- presDx_spo * abundDx_spo
    }

    #Southern Nearshore
    if("spn" %in% sampling){
      print("Fitting BRT-SPN")
      brt_dist_P_spn <- gbm.step(data=dat_hist_Dist_spn,
                                 gbm.x = c(25,27,28),
                                 gbm.y = 'pres',
                                 family = "bernoulli",
                                 tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                 plot.main=FALSE, verbose = FALSE)
      
      brt_dist_N_spn <- gbm.step(data=dat_hist_Dist_spn[dat_hist_Dist_spn$abundance>0,],
                                 gbm.x = c(25,27,28),
                                 gbm.y = 'log_abundance',
                                 family = "gaussian",
                                 tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                 plot.main=FALSE, verbose = FALSE)
      
      presDx_spn <- predict(brt_dist_P_spn, dat_hist, n.trees=brt_dist_P_spn$gbm.call$best.trees, type="response")
      abundDx_spn <- exp(predict(brt_dist_N_spn, dat_hist, n.trees=brt_dist_N_spn$gbm.call$best.trees, type="response"))
      dat_hist$brt_Dist_spn <- presDx_spn * abundDx_spn
      
      presDx_spn <- predict(brt_dist_P_spn, dat_fcast, n.trees=brt_dist_P_spn$gbm.call$best.trees, type="response")
      abundDx_spn <- exp(predict(brt_dist_N_spn, dat_fcast, n.trees=brt_dist_N_spn$gbm.call$best.trees, type="response"))
      dat_fcast$brt_Dist_spn <- presDx_spn * abundDx_spn
    }

    #All Offshore
    if("allo" %in% sampling){
      print("Fitting BRT-ALLO")
      brt_dist_P_allo <- gbm.step(data=dat_hist_Dist_allo,
                                 gbm.x = c(25,27,28),
                                 gbm.y = 'pres',
                                 family = "bernoulli",
                                 tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                 plot.main=FALSE, verbose = FALSE)
      
      brt_dist_N_allo <- gbm.step(data=dat_hist_Dist_allo[dat_hist_Dist_allo$abundance>0,],
                                 gbm.x = c(25,27,28),
                                 gbm.y = 'log_abundance',
                                 family = "gaussian",
                                 tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                 plot.main=FALSE, verbose = FALSE)
      
      presDx_allo <- predict(brt_dist_P_allo, dat_hist, n.trees=brt_dist_P_allo$gbm.call$best.trees, type="response")
      abundDx_allo <- exp(predict(brt_dist_N_allo, dat_hist, n.trees=brt_dist_N_allo$gbm.call$best.trees, type="response"))
      dat_hist$brt_Dist_allo <- presDx_allo * abundDx_allo
      
      presDx_allo <- predict(brt_dist_P_allo, dat_fcast, n.trees=brt_dist_P_allo$gbm.call$best.trees, type="response")
      abundDx_allo <- exp(predict(brt_dist_N_allo, dat_fcast, n.trees=brt_dist_N_allo$gbm.call$best.trees, type="response"))
      dat_fcast$brt_Dist_allo <- presDx_allo * abundDx_allo
    }

    #All Nearshore
    if("alln" %in% sampling){
      print("Fitting BRT-ALLN")
      brt_dist_P_alln <- gbm.step(data=dat_hist_Dist_alln,
                                 gbm.x = c(25,27,28),
                                 gbm.y = 'pres',
                                 family = "bernoulli",
                                 tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                 plot.main=FALSE, verbose = FALSE)
      
      brt_dist_N_alln <- gbm.step(data=dat_hist_Dist_alln[dat_hist_Dist_alln$abundance>0,],
                                 gbm.x = c(25,27,28),
                                 gbm.y = 'log_abundance',
                                 family = "gaussian",
                                 tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                 plot.main=FALSE, verbose = FALSE)
      
      presDx_alln <- predict(brt_dist_P_alln, dat_hist, n.trees=brt_dist_P_alln$gbm.call$best.trees, type="response")
      abundDx_alln <- exp(predict(brt_dist_N_alln, dat_hist, n.trees=brt_dist_N_alln$gbm.call$best.trees, type="response"))
      dat_hist$brt_Dist_alln <- presDx_alln * abundDx_alln
      
      presDx_alln <- predict(brt_dist_P_alln, dat_fcast, n.trees=brt_dist_P_alln$gbm.call$best.trees, type="response")
      abundDx_alln <- exp(predict(brt_dist_N_alln, dat_fcast, n.trees=brt_dist_N_alln$gbm.call$best.trees, type="response"))
      dat_fcast$brt_Dist_alln <- presDx_alln * abundDx_alln
    }
#
## Bycatch + Opt Target Sampling
#
    if("BY" %in% sampling){
      print("Fitting BRT-BY")
      brt_B_P <- gbm.step(data=dat_hist_BY,
                          gbm.x = c(25,27,28),
                          gbm.y = 'pres',
                          family = "bernoulli",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)

      brt_B_N <- gbm.step(data=dat_hist_BY[dat_hist_BY$abundance>0,], 
                          gbm.x = c(25,27,28),
                          gbm.y = 'log_abundance',
                          family = "gaussian",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)

      presBx <- predict(brt_B_P, dat_hist, n.trees=brt_B_P$gbm.call$best.trees, type="response")
      abundBx <- exp(predict(brt_B_N, dat_hist, n.trees=brt_B_N$gbm.call$best.trees, type="response"))
      dat_hist$brt_BY <- presBx * abundBx
      
      presBx <- predict(brt_B_P, dat_fcast, n.trees=brt_B_P$gbm.call$best.trees, type="response")
      abundBx <- exp(predict(brt_B_N, dat_fcast, n.trees=brt_B_N$gbm.call$best.trees, type="response"))
      dat_fcast$brt_BY<- presBx * abundBx
    }

#
###Closed Areas + Opt Target Species
#
    # Small Closed Area
    if("CA_sm" %in% sampling){
      print("Fitting BRT-CASM")
      brt_CA_P_sm <- gbm.step(data=dat_hist_CA_sm,
                           gbm.x = c(25,27,28),
                           gbm.y = 'pres',
                           family = "bernoulli",
                           tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                           plot.main=FALSE, verbose = FALSE)

      brt_CA_N_sm <- gbm.step(data=dat_hist_CA_sm[dat_hist_CA_sm$abundance>0,], 
                           gbm.x = c(25,27,28),
                           gbm.y = 'log_abundance',
                           family = "gaussian",
                           tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                           plot.main=FALSE, verbose = FALSE)

      presCAx_sm <- predict(brt_CA_P_sm, dat_hist, n.trees=brt_CA_P_sm$gbm.call$best.trees, type="response")
      abundCAx_sm <- exp(predict(brt_CA_N_sm, dat_hist, n.trees=brt_CA_N_sm$gbm.call$best.trees, type="response"))
      dat_hist$brt_CA_sm <- presCAx_sm * abundCAx_sm
      
      presCAx_sm <- predict(brt_CA_P_sm, dat_fcast, n.trees=brt_CA_P_sm$gbm.call$best.trees, type="response")
      abundCAx_sm <- exp(predict(brt_CA_N_sm, dat_fcast, n.trees=brt_CA_N_sm$gbm.call$best.trees, type="response"))
      dat_fcast$brt_CA_sm<- presCAx_sm * abundCAx_sm
    }

    # Medium Closed Area
    if("CA_med" %in% sampling){
      print("Fitting BRT-CAMED")
      brt_CA_P_med <- gbm.step(data=dat_hist_CA_med,
                              gbm.x = c(25,27,28),
                              gbm.y = 'pres',
                              family = "bernoulli",
                              tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                              plot.main=FALSE, verbose = FALSE)

      brt_CA_N_med <- gbm.step(data=dat_hist_CA_med[dat_hist_CA_med$abundance>0,], 
                              gbm.x = c(25,27,28),
                              gbm.y = 'log_abundance',
                              family = "gaussian",
                              tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                              plot.main=FALSE, verbose = FALSE)

      presCAx_med <- predict(brt_CA_P_med, dat_hist, n.trees=brt_CA_P_med$gbm.call$best.trees, type="response")
      abundCAx_med <- exp(predict(brt_CA_N_med, dat_hist, n.trees=brt_CA_N_med$gbm.call$best.trees, type="response"))
      dat_hist$brt_CA_med <- presCAx_med * abundCAx_med
      
      presCAx_med <- predict(brt_CA_P_med, dat_fcast, n.trees=brt_CA_P_med$gbm.call$best.trees, type="response")
      abundCAx_med <- exp(predict(brt_CA_N_med, dat_fcast, n.trees=brt_CA_N_med$gbm.call$best.trees, type="response"))
      dat_fcast$brt_CA_med<- presCAx_med * abundCAx_med
    }

    # Large Closed Area
    if("CA_lar" %in% sampling){
      print("Fitting BRT-CALAR")
      brt_CA_P_lar <- gbm.step(data=dat_hist_CA_lar,
                              gbm.x = c(25,27,28),
                              gbm.y = 'pres',
                              family = "bernoulli",
                              tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                              plot.main=FALSE, verbose = FALSE)

      brt_CA_N_lar <- gbm.step(data=dat_hist_CA_lar[dat_hist_CA_lar$abundance>0,], 
                              gbm.x = c(25,27,28),
                              gbm.y = 'log_abundance',
                              family = "gaussian",
                              tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                              plot.main=FALSE, verbose = FALSE)

      presCAx_lar <- predict(brt_CA_P_lar, dat_hist, n.trees=brt_CA_P_lar$gbm.call$best.trees, type="response")
      abundCAx_lar <- exp(predict(brt_CA_N_lar, dat_hist, n.trees=brt_CA_N_lar$gbm.call$best.trees, type="response"))
      dat_hist$brt_CA_lar <- presCAx_lar * abundCAx_lar
      
      presCAx_lar <- predict(brt_CA_P_lar, dat_fcast, n.trees=brt_CA_P_lar$gbm.call$best.trees, type="response")
      abundCAx_lar <- exp(predict(brt_CA_N_lar, dat_fcast, n.trees=brt_CA_N_lar$gbm.call$best.trees, type="response"))
      dat_fcast$brt_CA_lar<- presCAx_lar * abundCAx_lar
    }

##PLOTS
gbm.plot(brt_R_P, write.title=F, main="brt_R_P", plot.layout = c(2,2))
gbm.plot(brt_R_N, write.title=F, main="brt_R_N", plot.layout = c(2,2))
gbm.plot(brt_T_P_1, write.title=F, main="brt_T_P_1", plot.layout = c(2,2))
gbm.plot(brt_T_N_1, write.title=F, main="brt_T_N_1", plot.layout = c(2,2))
gbm.plot(brt_T_P_2, write.title=F, main="brt_T_P_2", plot.layout = c(2,2))
gbm.plot(brt_T_N_2, write.title=F, main="brt_T_N_2", plot.layout = c(2,2))
gbm.plot(brt_T_P_3, write.title=F, main="brt_T_P_3", plot.layout = c(2,2))
gbm.plot(brt_T_N_3, write.title=F, main="brt_T_N_3", plot.layout = c(2,2))
gbm.plot(brt_T_P_4, write.title=F, main="brt_T_P_4", plot.layout = c(2,2))
gbm.plot(brt_T_N_4, write.title=F, main="brt_T_N_4", plot.layout = c(2,2))
gbm.plot(brt_T_P_5, write.title=F, main="brt_T_P_5", plot.layout = c(2,2))
gbm.plot(brt_T_N_5, write.title=F, main="brt_T_N_5", plot.layout = c(2,2))
gbm.plot(brt_dist_P_npo, write.title=F, main="brt_dist_P_npo", plot.layout = c(2,2))
gbm.plot(brt_dist_N_npo, write.title=F, main="brt_dist_N_npo", plot.layout = c(2,2))
gbm.plot(brt_dist_P_mpo, write.title=F, main="brt_dist_P_mpo", plot.layout = c(2,2))
gbm.plot(brt_dist_N_mpo, write.title=F, main="brt_dist_N_mpo", plot.layout = c(2,2))
gbm.plot(brt_dist_P_spo, write.title=F, main="brt_dist_P_spo", plot.layout = c(2,2))
gbm.plot(brt_dist_N_spo, write.title=F, main="brt_dist_N_spo", plot.layout = c(2,2))
gbm.plot(brt_dist_P_allo, write.title=F, main="brt_dist_P_allo", plot.layout = c(2,2))
gbm.plot(brt_dist_N_allo, write.title=F, main="brt_dist_N_allo", plot.layout = c(2,2))
gbm.plot(brt_dist_P_npn, write.title=F, main="brt_dist_P_npn", plot.layout = c(2,2))
gbm.plot(brt_dist_N_npn, write.title=F, main="brt_dist_N_npn", plot.layout = c(2,2))
gbm.plot(brt_dist_P_mpn, write.title=F, main="brt_dist_P_mpn", plot.layout = c(2,2))
gbm.plot(brt_dist_N_mpn, write.title=F, main="brt_dist_N_mpn", plot.layout = c(2,2))
gbm.plot(brt_dist_P_spn, write.title=F, main="brt_dist_P_spn", plot.layout = c(2,2))
gbm.plot(brt_dist_N_spn, write.title=F, main="brt_dist_N_spn", plot.layout = c(2,2))
gbm.plot(brt_dist_P_alln, write.title=F, main="brt_dist_P_alln", plot.layout = c(2,2))
gbm.plot(brt_dist_N_alln, write.title=F, main="brt_dist_N_alln", plot.layout = c(2,2))
gbm.plot(brt_CA_P_sm, write.title=F, main="brt_CA_P_sm", plot.layout = c(2,2))
gbm.plot(brt_CA_N_sm, write.title=F, main="brt_CA_N_sm", plot.layout = c(2,2))
gbm.plot(brt_CA_P_med, write.title=F, main="brt_CA_P_med", plot.layout = c(2,2))
gbm.plot(brt_CA_N_med, write.title=F, main="brt_CA_N_med", plot.layout = c(2,2))
gbm.plot(brt_CA_P_lar, write.title=F, main="brt_CA_P_lar", plot.layout = c(2,2))
gbm.plot(brt_CA_N_lar, write.title=F, main="brt_CA_N_lar", plot.layout = c(2,2))

