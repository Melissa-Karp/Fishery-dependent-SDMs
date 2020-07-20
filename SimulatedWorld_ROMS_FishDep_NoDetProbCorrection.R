#This code is the operating model for a HMS archetype representing a Pelagic mobile predator-like species
#It uses average spring conditions from downscales ROMS projections using Hadley 
#IMPORTANT: Download average spring ROMS data from here: https://www.dropbox.com/sh/pj1vjz4h1n27hnb/AAD9sySDudMKmq3l4yLPvjQra?dl=0

#We include a trophic interaction between predator (PMP) and prey, but note the estmation model will use chl-a as a proxy for prey information:
#prey: distribution and suitability driven by SST and zooplankton integrated across top 200m
#predator: distribution and abundance driven by SST, MLD, and Species A
#Simulates fishery dependent sampling bias by building a fishing location suitability raster using virtual species
#Biases Include:
# Random Sampling (control)
# Preferential Sampling - fishing suitability a function of target species habitat suitability in Time Y
# Distance from Port - fishing location suitability function of distance from closest port
# Bycatch avoidance - adds behavior of fishermen actively trying to avoid areas of high bycatch risk 
# Spatial closure - fishing locations limited by spatial closure 

dir <- "~/DisMAP project/Location, Location, Location/Location Workshop/ROMS/hadley"

SimulateWorld_ROMS_FishDepFun_NODetProb <- function(dir, nsamples){
  #dir is the local directory that points to where ROMS data is stored 
  #nsamples let's you choose the number of samples each year
  
  #required libraries
  library(virtualspecies)
  library(scales)
  library(raster)
  library(glmmfields)
  library(dplyr)
  
  #----Create output file----
  #This will be the information passed to the estimation model
  output <- as.data.frame(matrix(NA, nrow=21912*121,ncol=14)) #21912 non-NA grid cells in ROMS
  colnames(output) <- c("lon","lat","year","pres_t","suitability_t","suitability_t1","random_sampled","pref_sampled", "Dist_sampled", "BY_sampled", "ByDist_sampled", "sst","zoo_200", "mld")
  
  #----Load in rasters----
  #These are the average spring conditions from the downscaled 'had' earth system model
  files_sst <- list.files(paste0(dir,'/sst_spring_avg'), full.names = TRUE, pattern=".grd") #should be 1452 files
  files_mld <- list.files(paste0(dir,'/ild_0.5C'), full.names = TRUE, pattern=".grd")
  files_zoo <- list.files(paste0(dir,'/zoo_200m'), full.names = TRUE, pattern=".grd")
  files_chl_surface <- list.files(paste0(dir,'/chl_surface'), full.names = TRUE, pattern=".grd")
  years <- seq(1980,2100,1)
  
  #---Load in Distance to port file -----
  dist_to_ports<-read.csv("~/DisMAP project/Location, Location, Location/Location Workshop/Dist_to_Ports.csv")
  
  #----Loop through each year----
  for (y in 1:121){
    print(paste0("Creating environmental simulation for Year ",years[y]))
    
    #Load in environmental rasters for a specific year
    sst <- raster(files_sst[y])
    mld <- raster(files_mld[y])
    zoo <- raster(files_zoo[y])
    chl_surface<-raster(files_chl_surface[y])
    chla_surface <- log(chla_surface)
    
    #Optional: plot environmental layers
    # par(mfrow=c(1,4))
    # plot(sst, main= 'SST')
    # plot(chl_surface, main = 'Chl-a Surface')
    # plot(mld, main = 'MLD')
    # plot(zoo, main = 'Zooplankton')
    
    #Environmental at time t-1 (previous year).
    #If year 1, then just use at time t.
    if (y>1){
      sst_t1 <- raster(files_sst[y-1])
      mld_t1<-raster(files_mld[y-1])
      zoo_t1<-raster(files_zoo[y-1])
      
    } else {
      sst_t1 <- sst
      mld_t1 <- mld
      zoo_t1 <- zoo
    }
    
    #---- ASSIGN RESPONSE CURVES ---- ####
    
    #----SPECIES A (prey): assign response curves----
    #Species A: likes high zooplankton and medium temps
    
    #Stack rasters
    spA_stack <- stack(sst, zoo)
    names(spA_stack) <- c('sst', 'zoo')
    
    #Assign preferences
    spA_parameters <- formatFunctions(sst = c(fun="dnorm",mean=15,sd=5),
                                      zoo = c(fun="logisticFun",alpha=-6,beta=50))
    spA_suitability <- generateSpFromFun(spA_stack,parameters=spA_parameters, rescale = FALSE,rescale.each.response = FALSE) #Important: make sure rescaling is false. Doesn't work well in the 'for' loop.
    # plot(spA_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spA_suitability) #plot response curves
    
    #manually rescale
    ref_max_sst <- dnorm(spA_parameters$sst$args[1], mean=spA_parameters$sst$args[1], sd=spA_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_zoo <- 1 / (1 + exp(((zoo@data@max) - spA_parameters$zoo$args[2])/spA_parameters$zoo$args[1]))
    ref_max <- ref_max_sst * ref_max_zoo #simple multiplication of layers.
    spA_suitability$suitab.raster <- (1/ref_max)*spA_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
    # spA_suitability$suitab.raster
    # plot(spA_suitability$suitab.raster) #plot habitat suitability
    
    
    #----SPECIES B (Pelagic mobile predator): assign response curves----
    #Species B: likes to eat Species A, and warmer temperatues & shallow MLD
    
    #Stack rasters
    spB_stack <- stack(sst, mld, spA_suitability$suitab.raster)
    names(spB_stack) <- c('sst',"mld", "spA")
    
    #Assign preferences
    spB_parameters <- formatFunctions(sst = c(fun="dnorm",mean=17,sd=5),
                                      mld = c(fun="dnorm",mean=50,sd=25),
                                      spA = c(fun="logisticFun",alpha=-0.05,beta=0.5))
    spB_suitability <- generateSpFromFun(spB_stack,parameters=spB_parameters, rescale = FALSE,rescale.each.response = FALSE)
    #plot(spB_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spB_suitability) #plot response curves
    
    #manually rescale
    ref_max_sst <- dnorm(spB_parameters$sst$args[1], mean=spB_parameters$sst$args[1], sd=spB_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_mld <- dnorm(spB_parameters$mld$args[1], mean=spB_parameters$mld$args[1], sd=spB_parameters$mld$args[2])
    ref_max_spA <- 1 / (1 + exp(((spA_suitability$suitab.raster@data@max) - spB_parameters$spA$args[2])/spB_parameters$spA$args[1]))
    ref_max <- ref_max_sst * ref_max_mld * ref_max_spA
    spB_suitability$suitab.raster <- (1/ref_max)*spB_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum temp is encountered
    # print(spB_suitability$suitab.raster)
    plot(spB_suitability$suitab.raster) #plot habitat suitability
    
    #-------Species B: Suitability for PREVIOUS YEAR (y-1)---------#
    #Sp A suitability y-1
    
    #Stack rasters
    spA_stack_t1 <- stack(sst_t1, zoo_t1)
    names(spA_stack_t1) <- c('sst_t1', 'zoo_t1')
    
    #Assign preferences
    spA_parameters_t1 <- formatFunctions(sst_t1 = c(fun="dnorm",mean=15,sd=5),
                                         zoo_t1 = c(fun="logisticFun",alpha=-6,beta=50))
    spA_suitability_t1 <- generateSpFromFun(spA_stack_t1,parameters=spA_parameters_t1, rescale = FALSE,rescale.each.response = FALSE) #Important: make sure rescaling is false. Doesn't work well in the 'for' loop.
    # plot(spA_suitability_t1$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spA_suitability_t1) #plot response curves
    
    #manually rescale
    ref_max_sst_t1 <- dnorm(spA_parameters_t1$sst_t1$args[1], mean=spA_parameters_t1$sst_t1$args[1], sd=spA_parameters_t1$sst_t1$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_zoo_t1 <- 1 / (1 + exp(((zoo_t1@data@max) - spA_parameters_t1$zoo_t1$args[2])/spA_parameters_t1$zoo_t1$args[1]))
    ref_max_t1 <- ref_max_sst_t1 * ref_max_zoo_t1 #simple multiplication of layers.
    spA_suitability_t1$suitab.raster <- (1/ref_max_t1)*spA_suitability_t1$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
    # spA_suitability_t1$suitab.raster
    #plot(spA_suitability_t1$suitab.raster) #plot habitat suitability
    
    #Sp B suitability in year y-1
    spB_stack_t1 <- stack(sst_t1, mld_t1, spA_suitability_t1$suitab.raster)
    names(spB_stack_t1) <- c('sst_t1',"mld_t1", "spA_t1")
    
    #Assign preferences
    spB_parameters_t1 <- formatFunctions(sst_t1 = c(fun="dnorm",mean=17,sd=5),
                                         mld_t1 = c(fun="dnorm",mean=50,sd=25),
                                         spA_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.5))
    spB_suitability_t1 <- generateSpFromFun(spB_stack_t1,parameters=spB_parameters_t1, rescale = FALSE,rescale.each.response = FALSE)
    #plot(spB_suitability_t1$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spB_suitability) #plot response curves
    
    #manually rescale
    ref_max_sst_t1 <- dnorm(spB_parameters_t1$sst_t1$args[1], mean=spB_parameters_t1$sst_t1$args[1], sd=spB_parameters_t1$sst_t1$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_mld_t1 <- dnorm(spB_parameters_t1$mld_t1$args[1], mean=spB_parameters_t1$mld_t1$args[1], sd=spB_parameters_t1$mld_t1$args[2])
    ref_max_spA_t1 <- 1 / (1 + exp(((spA_suitability_t1$suitab.raster@data@max) - spB_parameters_t1$spA_t1$args[2])/spB_parameters_t1$spA_t1$args[1]))
    ref_max_t1 <- ref_max_sst_t1 * ref_max_mld_t1 * ref_max_spA_t1
    spB_suitability_t1$suitab.raster <- (1/ref_max_t1)*spB_suitability_t1$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum temp is encountered
    # print(spB_suitability_t1$suitab.raster)
    plot(spB_suitability_t1$suitab.raster) #plot habitat suitability
    
    #-----Species C: Turtle-like (BYCATCH SPECIES)------
    #Species C: likes high zooplankton and warm temps. 
    
    #Stack rasters
    spC_stack <- stack(sst, zoo)
    names(spC_stack) <- c('sst', 'zoo')
    
    #Assign preferences
    spC_parameters <- formatFunctions(sst = c(fun="dnorm",mean=25,sd=10),
                                      zoo = c(fun="logisticFun",alpha=-6,beta=50))
    spC_suitability <- generateSpFromFun(spC_stack,parameters=spC_parameters, rescale = FALSE,rescale.each.response = FALSE) #Important: make sure rescaling is false. Doesn't work well in the 'for' loop.
    #plot(spC_suitability$suitab.raster) #plot habitat suitability
    #virtualspecies::plotResponse(spC_suitability) #plot response curves
    
    #manually rescale
    ref_max_sst <- dnorm(spC_parameters$sst$args[1], mean=spC_parameters$sst$args[1], sd=spC_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_zoo <- 1 / (1 + exp(((zoo@data@max) - spC_parameters$zoo$args[2])/spC_parameters$zoo$args[1]))
    ref_max <- ref_max_sst * ref_max_zoo #simple multiplication of layers.
    spC_suitability$suitab.raster <- (1/ref_max)*spC_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
    plot(spC_suitability$suitab.raster) #plot habitat suitability
    
    #----CONVERT SP B SUITABILITY TO PRESENCE-ABSENCE----####
    
    #Use a specific function to convert suitability (0-1) to presence or absence (1 or 0)
    set.seed(y) #*********set seed to change every year, based on value of y. This will keep within year constant, but allow for noise between years
    #this will insure that the PA raster created for year Y and year Y-1 will be the same as we are using the pa.method = "probability"
    suitability_PA <- virtualspecies::convertToPA(spB_suitability, PA.method = "probability", beta = 0.5,
                                                  alpha = -0.05, species.prevalence = NULL, plot = FALSE)
    # plotSuitabilityToProba(suitability_PA) #Let's you plot the shape of conversion function
    # plot(suitability_PA$pa.raster)
    
    #-------BUILD Fishing Suitability Rasters for BIASING SAMPLING--------------------######
    
    #get "minimum distance from a port" for each cell
    dist_to_ports$dist_min <- apply(dist_to_ports[,c("dp1", "dp2",
                                                     "dp3", "dp3",
                                                     "dp5")], 1, FUN=min)
    dist_min_raster<- rasterFromXYZ(dist_to_ports[,c("lon","lat","dist_min")])
    dist_min_km<-dist_min_raster/1000
    # plot(dist_min_rescale)
    # plot(dist_min_raster)
    
    ##### 1. Just Target Species Habitat Suitability considered 
    Tar_stack <- stack(spB_suitability_t1$suitab.raster)
    names(Tar_stack) <- c('spB_t1')
    
    #Assign preferences
    Tar_params <- formatFunctions(spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.5))
    Tar_suitability <- generateSpFromFun(Tar_stack,parameters=Tar_params, rescale = FALSE,rescale.each.response = FALSE)
    #plot(Tar_suitability$suitab.raster) 
    #virtualspecies::plotResponse(Tar_suitability) #plot response curves
    #manually rescale
    Tar_rescale<- reclassify(Tar_suitability$suitab.raster, c(0.01, 0.25, 0.01))
    Tar_reclassify<- Tar_rescale/maxValue(Tar_rescale)
    #plot(Tar_rescale)
    plot(Tar_reclassify)
    
    ##### 2. Distance to Port and Target Species Suitability considered
    Dist_stack <- stack(spB_suitability_t1$suitab.raster, dist_min_km)
    names(Dist_stack) <- c('spB_t1', "min_dist")
    
    #Assign preferences
    linear.function <- function(x, a, b)
    {
      a*x + b
    }
    Dist_params <- formatFunctions(min_dist = c(fun="linear.function",a=-0.0007, b=1),
                                   spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.5))
    Dist_suitability <- generateSpFromFun(Dist_stack,parameters=Dist_params, rescale = FALSE,rescale.each.response = FALSE)
    #plot(Dist_suitability$suitab.raster) #plot habitat suitability
    Dist_rescale<- reclassify(Dist_suitability$suitab.raster, c(0.01, 0.25, 0.01))
    #plot(Dist_rescale)
    Dist_reclassify<-Dist_rescale/maxValue(Dist_rescale)
    plot(Dist_reclassify)
    
    ##### 3. BYCATCH RISK and target spp suitability considered
    Bycatch_stack <- stack(spB_suitability_t1$suitab.raster, spC_suitability$suitab.raster)
    names(Bycatch_stack) <- c('spB_t1', "spC")
    
    #Assign preferences
    Bycatch_params <- formatFunctions(spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.5),
                                      spC = c(fun="logisticFun", alpha= 0.05, beta=0.5))
    Bycatch_suitability <- generateSpFromFun(Bycatch_stack,parameters=Bycatch_params, rescale = FALSE,rescale.each.response = FALSE)
    #plot(Bycatch_suitability$suitab.raster) 
    #virtualspecies::plotResponse(Bycatch_suitability) #plot response curves
    
    BY_rescale <- reclassify(Bycatch_suitability$suitab.raster, c(0.01, 0.25, 0.01))
    #plot(BY_rescale)
    By_reclassify<-BY_rescale/maxValue(BY_rescale)
    plot(By_reclassify)
    
    ##### 4. Bycatch risk AND Distance to Port AND Target Species Suitability considered 
    ByDist_stack <- stack(spB_suitability_t1$suitab.raster, spC_suitability$suitab.raster, dist_min_km)
    names(ByDist_stack) <- c('spB_t1', "spC", "min_dist")
    
    #Assign preferences
    ByDist_params <- formatFunctions(min_dist = c(fun="linear.function",a=-0.0007, b=1),
                                     spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.5),
                                     spC = c(fun="logisticFun", alpha= 0.05, beta=0.5))
    ByDist_suitability <- generateSpFromFun(ByDist_stack,parameters=ByDist_params, rescale = FALSE,rescale.each.response = FALSE)
    #plot(ByDist_suitability$suitab.raster) #RASTER USED TO WEIGHT SAMPLE OCCURRENCES
    #virtualspecies::plotResponse(ByDist_suitability) #plot response curves
    BD_rescale<-reclassify(ByDist_suitability$suitab.raster, c(0.01, 0.25, 0.01))
    #plot(BD_rescale)
    ByDist_reclassify<-BD_rescale/maxValue(BD_rescale)
    plot(ByDist_reclassify)
    
    ###-----SAMPLE PRESENCES AND ABSENCES-----#####
    
    #******Random Sampling of nsamples*******
    presence.points.random <- sampleOccurrences(suitability_PA,n = nsamples,type = "presence-absence",
                                                detection.probability = 1, error.probability=0, plot = TRUE,
                                                sample.prevalence = NULL)
    
    #convert to dataframe
    pres_df_random <- cbind(as.data.frame(presence.points.random$sample.points$x),as.data.frame(presence.points.random$sample.points$y), as.data.frame(presence.points.random$sample.points$Observed))
    #pres_df_random <- cbind(as.data.frame(Ran_data$lons[Ran_data$cell_num %in% cells_ran]),as.data.frame(Ran_data$lats[Ran_data$cell_num %in% cells_ran]))
    colnames(pres_df_random) <- c("x","y", "Obs_random")
    pres_df_random$random_sampled <- 1
    #expand dataframe to include all possible locations
    df_full <- as.data.frame(rasterToPoints(sst)[,1:2]) #picking an example raster to extract lat and lon from
    df_full_2 <- left_join(df_full, pres_df_random, by=c('x','y'))
    df_full_2$random_sampled <- ifelse(is.na(df_full_2$random_sampled),0,df_full_2$random_sampled)
    
    #****** 1. Preferential Sampling of nsamples - based on HABITAT SUITABILITY of target species**********
    presence.points.pref<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                            detection.probability = 1, bias = "manual", weights = Tar_reclassify, plot = TRUE)
    #convert to dataframe
    pres_df_pref <- cbind(as.data.frame(presence.points.pref$sample.points$x),as.data.frame(presence.points.pref$sample.points$y), as.data.frame(presence.points.pref$sample.points$Observed))
    #pres_df_pref <- cbind(as.data.frame(Tar_data$lons[Tar_data$cell_num %in% cells_tar]),as.data.frame(Tar_data$lats[Tar_data$cell_num %in% cells_tar]))
    colnames(pres_df_pref) <- c("x","y", "Obs_pref")
    pres_df_pref$pref_sampled <- 1
    #add to dataframe
    df_full_3 <- left_join(df_full_2, pres_df_pref, by=c('x','y'))
    df_full_3$pref_sampled <- ifelse(is.na(df_full_3$pref_sampled),0,df_full_3$pref_sampled)
    
    
    #******* 2. Sampling of nsamples based on Habitat Suitability & DISTANCE TO PORT********** 
    presence.points.dist<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                            detection.probability = 1, bias = "manual", weights= Dist_reclassify, plot = TRUE)
    #convert to dataframe
    pres_df_dist <- cbind(as.data.frame(presence.points.dist$sample.points$x),as.data.frame(presence.points.dist$sample.points$y), as.data.frame(presence.points.dist$sample.points$Observed))
    #pres_df_dist <- cbind(as.data.frame(Dist_data$lons[Dist_data$cell_num %in% cells_dist]),as.data.frame(Dist_data$lats[Dist_data$cell_num %in% cells_dist]))
    colnames(pres_df_dist) <- c("x","y", "Obs_dist")
    pres_df_dist$dist_sampled <- 1
    #add to dataframe
    df_full_4 <- left_join(df_full_3, pres_df_dist, by=c('x','y'))
    df_full_4$dist_sampled <- ifelse(is.na(df_full_4$dist_sampled),0,df_full_4$dist_sampled)
    
    
    #******* 3. Sampling of nsamples based on Habitat Suitability & BYCATCH RISK********** 
    presence.points.Bycatch<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                               detection.probability = 1, bias = "manual", weights= By_reclassify, plot = TRUE)
    #convert to dataframe
    pres_df_BY <- cbind(as.data.frame(presence.points.Bycatch$sample.points$x),as.data.frame(presence.points.Bycatch$sample.points$y), as.data.frame(presence.points.Bycatch$sample.points$Observed))
    #pres_df_BY <- cbind(as.data.frame(BY_data$lons[BY_data$cell_num %in% cells_BY]),as.data.frame(BY_data$lats[BY_data$cell_num %in% cells_BY]))
    colnames(pres_df_BY) <- c("x","y", "Obs_BY")
    pres_df_BY$BY_sampled <- 1
    #add to dataframe
    df_full_5 <- left_join(df_full_4, pres_df_BY, by=c('x','y'))
    df_full_5$BY_sampled <- ifelse(is.na(df_full_5$BY_sampled),0,df_full_5$BY_sampled)
    
    
    #******* 4. Sampling of nsamples based on Habitat Suitability & DISTANCE to Port & BYCATCH RISK********** 
    presence.points.ByDist<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                              detection.probability = 1, bias = "manual", weights= ByDist_reclassify, plot = TRUE)
    #convert to dataframe
    pres_df_ByDist <- cbind(as.data.frame(presence.points.ByDist$sample.points$x),as.data.frame(presence.points.ByDist$sample.points$y), as.data.frame(presence.points.ByDist$sample.points$Observed))
    #pres_df_ByDist <- cbind(as.data.frame(ByDist_data$lons[ByDist_data$cell_num %in% cells_BD]),as.data.frame(ByDist_data$lats[ByDist_data$cell_num %in% cells_BD]))
    colnames(pres_df_ByDist) <- c("x","y", "Obs_ByDist")
    pres_df_ByDist$ByDist_sampled <- 1
    #add to dataframe
    df_full_6 <- left_join(df_full_5, pres_df_ByDist, by=c('x','y'))
    df_full_6$ByDist_sampled <- ifelse(is.na(df_full_6$ByDist_sampled),0,df_full_6$ByDist_sampled)
    
    #----EXTRACT DATA for each year----
    
    print("Extracting suitability")
    ei <- 21912*y #end location in output grid to index to
    se <- ei - (21912-1) #start location in output grid to index to
    output$lat[se:ei] <- df_full_6$y
    output$lon[se:ei] <- df_full_6$x
    output$year[se:ei] <- rep(years[y],21912)
    output$pres_t[se:ei] <- rasterToPoints(suitability_PA$pa.raster)[,3]
    output$suitability_t[se:ei] <- rasterToPoints(spB_suitability$suitab.raster)[,3]  #extract points from suitability file
    output$suitability_t1[se:ei] <- rasterToPoints(spB_suitability_t1$suitab.raster)[,3]
    output$Obs_random <- df_full_6$Obs_random
    output$Obs_pref <- df_full_6$Obs_pref
    output$Obs_dist<-df_full_6$Obs_dist
    output$Obs_BY<-df_full_6$Obs_BY
    output$Obs_ByDist<-df_full_6$Obs_ByDist
    output$random_sampled <- df_full_6$random_sampled
    output$pref_sampled <- df_full_6$pref_sampled
    output$Dist_sampled<-df_full_6$dist_sampled
    output$BY_sampled<-df_full_6$BY_sampled
    output$ByDist_sampled<-df_full_6$ByDist_sampled
    output$sst[se:ei] <-  rasterToPoints(sst)[,3]   #extract points from suitability file
    output$zoo_200[se:ei] <-  rasterToPoints(zoo)[,3]
    output$mld[se:ei] <-  rasterToPoints(mld)[,3]
    output$chl_surface[se:ei]<-rasterToPoints(chl_surface)[,3]
  }
  
  #Average monthly biomass available to CCS is: 1.18x10^5 ± (0.13x10^5 se) mt (from Desiree Tommasi)
  mean_spatial <- round(118000/140, 1)
  se_spatial <- round((13000/140) ,2)
  output$abundance_t <- ifelse(output$pres_t==1,rnorm(nrow(output),mean_spatial, se_spatial)*output$suitability_t,0)
  output$abund_Rand <- ifelse(output$Obs_random==1, output$abundance_t,0)
  output$abund_pref <- ifelse(output$Obs_pref==1, output$abundance_t,0)
  output$abund_dist <- ifelse(output$Obs_dist==1, output$abundance_t,0)
  output$abund_BY <- ifelse(output$Obs_BY==1, output$abundance_t,0)
  output$abund_ByDist <- ifelse(output$Obs_ByDist==1, output$abundance_t,0)
  
  ##CPUE=aN^(B+1) => hyperstability 
  # catch_1<-'^'(output$abundance_t, 0.4)
  # output$catch_hyp<-2.5*catch_1
  # #plot(output$abundance_t, output$catch)
  # output$catch_Rand <- ifelse(output$Obs_random==1, output$catch_hyp,0)
  # output$catch_pref <- ifelse(output$Obs_pref==1, output$catch_hyp,0)
  # output$catch_dist <- ifelse(output$Obs_dist==1, output$catch_hyp,0)
  # output$catch_BY <- ifelse(output$Obs_BY==1, output$catch_hyp,0)
  # output$catch_ByDist <- ifelse(output$Obs_ByDist==1, output$catch_hyp,0)
  
  
  # #--------PRINT and SAVE CSV------------------#####
  
  print('Saving csv to working directory')
  write.csv(output, 'FisheryDependent_OM_Simulation.csv',row.names = FALSE)
  return(output)
}

#test<-SimulateWorld_ROMS_FishDepFun_NODetProb(dir=dir, nsamples = 100)
