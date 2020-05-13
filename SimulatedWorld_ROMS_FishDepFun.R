#This code is the operating model for a HMS archetype representing a Pelagic mobile predator-like species
#It uses average spring conditions from downscales ROMS projections using Hadley 
#IMPORTANT: Download average spring ROMS data from here: https://www.dropbox.com/sh/pj1vjz4h1n27hnb/AAD9sySDudMKmq3l4yLPvjQra?dl=0

#We include a trophic interaction between predator (PMP) and prey, but note the estmation model will use chl-a as a proxy for prey information:
#prey: distribution and suitability driven by SST and zooplankton integrated across top 200m
#predator: distribution and abundance driven by SST, MLD, and Species A
#Simulates fishery dependent sampling bias where fishing locations (sample locations) is skewed to be more common
#in areas where fish abundance was high in the year before (y-1 suitability was high)

dir <- "~/DisMAP project/Location, Location, Location/Location Workshop/ROMS"

SimulateWorld_ROMS_PMP <- function(dir, nsamples){
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
  colnames(output) <- c("lon","lat","year","pres_t","pres_t1","suitability_t","suitability_t1","sst","zoo_200","chla_surface", "mld","random_sampled","pref_sampled", "RUM_sampled")
  
  #----Load in rasters and datasets----
  #These are the average spring conditions from the downscaled gfdl earth system model
  files_sst <- list.files(paste0(dir,'/sst_monthly'), full.names = TRUE, pattern=".grd") #should be 1452 files
  files_chl_surface <- list.files(paste0(dir,'/chl_surface'), full.names = TRUE, pattern=".grd") #should be 1452 files
  files_mld <- list.files(paste0(dir,'/ild_0.5C'), full.names = TRUE, pattern=".grd")
  files_zoo <- list.files(paste0(dir,'/zoo_200m'), full.names = TRUE, pattern=".grd")
  years <- seq(1980,2100,1)
  
  #load distance to ports dataset
  dist_to_ports<-read.csv("~/DisMAP project/Location, Location, Location/Location Workshop/Dist_to_Ports.csv")

  #----Loop through each year----
  for (y in 1:121){
    print(paste0("Creating environmental simulation for Year ",years[y]))
    
    #Load in environmental rasters for a specific year
    sst <- raster(files_sst[y])
    mld <- raster(files_mld[y])
    zoo <- raster(files_zoo[y])
    chla_surface <- raster(files_chl_surface[y])
    chla_surface <- log(chla_surface)
    
    #Optional: plot environmental layers
    # par(mfrow=c(1,4))
    # plot(sst, main= 'SST')
    # plot(chla_surface, main = 'Chl-a Surface')
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
    # plot(spB_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spB_suitability) #plot response curves
    
    #manually rescale
    ref_max_sst <- dnorm(spB_parameters$sst$args[1], mean=spB_parameters$sst$args[1], sd=spB_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_mld <- dnorm(spB_parameters$mld$args[1], mean=spB_parameters$mld$args[1], sd=spB_parameters$mld$args[2])
    ref_max_spA <- 1 / (1 + exp(((spA_suitability$suitab.raster@data@max) - spB_parameters$spA$args[2])/spB_parameters$spA$args[1]))
    ref_max <- ref_max_sst * ref_max_mld * ref_max_spA
    spB_suitability$suitab.raster <- (1/ref_max)*spB_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum temp is encountered
    # print(spB_suitability$suitab.raster)
    # plot(spB_suitability$suitab.raster) #plot habitat suitability
    
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
    # plot(spB_suitability_t1$suitab.raster) #plot habitat suitability
    
    #-----Species C: Turtle-like------
    #Species C: likes high zooplankton and warm temps. 
    
    #Stack rasters
    spC_stack <- stack(sst, zoo)
    names(spC_stack) <- c('sst', 'zoo')
    
    #Assign preferences
    spC_parameters <- formatFunctions(sst = c(fun="dnorm",mean=25,sd=10),
                                      zoo = c(fun="logisticFun",alpha=-6,beta=50))
    spC_suitability <- generateSpFromFun(spC_stack,parameters=spC_parameters, rescale = FALSE,rescale.each.response = FALSE) #Important: make sure rescaling is false. Doesn't work well in the 'for' loop.
    # plot(spC_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spC_suitability) #plot response curves
    
    #manually rescale
    ref_max_sst <- dnorm(spC_parameters$sst$args[1], mean=spC_parameters$sst$args[1], sd=spC_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_zoo <- 1 / (1 + exp(((zoo@data@max) - spC_parameters$zoo$args[2])/spC_parameters$zoo$args[1]))
    ref_max <- ref_max_sst * ref_max_zoo #simple multiplication of layers.
    spC_suitability$suitab.raster <- (1/ref_max)*spC_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
    # spC_suitability$suitab.raster
    # plot(spC_suitability$suitab.raster) #plot habitat suitability
    
    #----CONVERT SP B SUITABILITY TO PRESENCE-ABSENCE----####
    
    #Use a specific function to convert suitability (0-1) to presence or absence (1 or 0)
    set.seed(y) #*********set seed to change every year, based on value of y. This will keep within year constant, but allow for noise between years
    #this will insure that the PA raster created for year Y and year Y-1 will be the same as we are using the pa.method = "probability"
    suitability_PA <- virtualspecies::convertToPA(spB_suitability, PA.method = "probability", beta = 0.5,
                                                  alpha = -0.05, species.prevalence = NULL, plot = FALSE)
    # plotSuitabilityToProba(suitability_PA) #Let's you plot the shape of conversion function
    # plot(suitability_PA$pa.raster)
    
    #y-1 presence absence 
    set.seed(y)
    suitability_PA_t1 <- virtualspecies::convertToPA(spB_suitability_t1, PA.method = "probability", beta = 0.5,
                                                     alpha = -0.05, species.prevalence = NULL, plot = FALSE)
    # plot(suitability_PA_t1$pa.raster)
    
   
       ##-------BUILD Random Utility Function--------------------######
 
   #add presence, suitability, and abundance values from [y-1] to the dist_to_ports dataframe
    dist_to_ports$pres_t1 <- extract(suitability_PA_t1$pa.raster, dist_to_ports[,c("lon","lat")])  #*** JS
    dist_to_ports$suit_t1 <- extract(spB_suitability_t1$suitab.raster, dist_to_ports[,c("lon","lat")])  #*** JS

    mean_spatial <- round(118000/140, 1)
    se_spatial <- round((13000/140) ,2)
    dist_to_ports$abund_t1 <- ifelse(dist_to_ports$pres_t1==1, rnorm(nrow(dist_to_ports), mean_spatial, se_spatial)*dist_to_ports$suit_t1, 0)

    #calculate utility [note: we assume vessels return to port of departure, hence '*2']
    price<-4000 #just picked a random number for now
    cost_per_km<-1500 #random number for now  *** JS
    # port 1
    dist_to_ports$utility_p1 <- price*dist_to_ports$abund_t1 - #revenue
      dist_to_ports$dp1/1000*cost_per_km*2             #cost
    plot(rasterFromXYZ(dist_to_ports[,c("lon","lat","utility_p1")]))  #*** JS: plot P1 utility
    #points(-117.1441, 32.6717, pch=0, cex=2)
    # port 2
    dist_to_ports$utility_p2 <- price*dist_to_ports$abund_t1 - #revenue  #*** JS
      (dist_to_ports$dp2/1000)*cost_per_km*2             #cost
    plot(rasterFromXYZ(dist_to_ports[,c("lon","lat","utility_p2")]))
    # port 3
    dist_to_ports$utility_p3 <- price*dist_to_ports$abund_t1 - #revenue
      (dist_to_ports$dp3/1000)*cost_per_km*2             #cost
    plot(rasterFromXYZ(dist_to_ports[,c("lon","lat","utility_p3")]))  #*** JS: plot P1 utility
    # port 4
    dist_to_ports$utility_p4 <- price*dist_to_ports$abund_t1 - #revenue
      (dist_to_ports$dp4/1000)*cost_per_km*2             #cost
    plot(rasterFromXYZ(dist_to_ports[,c("lon","lat","utility_p4")]))  #*** JS: plot P1 utility
    # port 5
    dist_to_ports$utility_p5 <- price*dist_to_ports$abund_t1 - #revenue
      (dist_to_ports$dp5/1000)*cost_per_km*2
    plot(rasterFromXYZ(dist_to_ports[,c("lon","lat","utility_p5")]))  #*** JS: plot P1 utility

    dist_to_ports$utility_max <- apply(dist_to_ports[,c("utility_p1", "utility_p2",
                                                         "utility_p3", "utility_p3",
                                                         "utility_p5")], 1, FUN=max)

    dist_to_ports$utility_SUM <- apply(dist_to_ports[,c("utility_p1", "utility_p2",
                                                        "utility_p3", "utility_p3",
                                                         "utility_p5")], 1, FUN=sum)

    # head(dist_to_ports)
    # # take make of all ports for coast-wide utility
    df_util_raster <- rasterFromXYZ(dist_to_ports[,c("lon","lat","utility_max")])  #*** JS: save coast-wide utility
    plot(df_util_raster, asp=1, main="Coast-wide Utility")

    df_sum_raster<-rasterFromXYZ(dist_to_ports[,c("lon","lat","utility_SUM")])  #*** JS: save coast-wide utility
    plot(df_sum_raster, asp=1, main="Coast-wide Utility")

    # # *** Rescale utility_max between 0 and 1
    df_util_raster2 <- df_util_raster + abs(minValue(df_util_raster))
    df_util_raster2 <- df_util_raster2/maxValue(df_util_raster2)  ## df_util_raster2 is the raster we use in the next section to sample presence-absences
    plot(df_util_raster2, asp=1, main="Coast-wide Utility")
    points(-117.1441, 32.6717, pch=0, cex=2)
    points(-122.001620, 36.965719, pch=0, cex=2)
    points(-123.050618, 38.334302, pch=0, cex=2)
    points(-124.292000, 43.383975, pch=0, cex=2)
    points(-124.114934, 46.911534, pch=0, cex=2)

    # # *** Rescale utility_SUM between 0 and 1
    df_sum_raster2 <- df_sum_raster + abs(minValue(df_sum_raster))
    df_sum_raster2 <- df_sum_raster2/maxValue(df_sum_raster2)  ## df_util_raster2 is the raster we use in the next section to sample presence-absences
    plot(df_sum_raster2, asp=1, main="Coast-wide Utility")
    points(-117.1441, 32.6717, pch=0, cex=2)
    points(-122.001620, 36.965719, pch=0, cex=2)
    points(-123.050618, 38.334302, pch=0, cex=2)
    points(-124.292000, 43.383975, pch=0, cex=2)
    points(-124.114934, 46.911534, pch=0, cex=2)

   ## Add cost of catching BYCATCH species into RUM **********
    dist_to_ports$suitSpC <- extract(spC_suitability$suitab.raster, dist_to_ports[,c("lon","lat")])  #*** JS
    bycost=1000
    # port 1
    dist_to_ports$utility_p1_BY <- price*dist_to_ports$abund_t1 - #revenue
      dist_to_ports$dp1/1000*cost_per_km*2 + dist_to_ports$suitSpC*bycost           #cost
    plot(rasterFromXYZ(dist_to_ports[,c("lon","lat","utility_p1_BY")]))  #*** JS: plot P1 utility
    #points(-117.1441, 32.6717, pch=0, cex=2)
    # port 2
    dist_to_ports$utility_p2_BY <- price*dist_to_ports$abund_t1 - #revenue  #*** JS
      (dist_to_ports$dp2/1000)*cost_per_km*2 + dist_to_ports$suitSpC*bycost        #cost
    plot(rasterFromXYZ(dist_to_ports[,c("lon","lat","utility_p2_BY")]))
    # port 3
    dist_to_ports$utility_p3_BY <- price*dist_to_ports$abund_t1 - #revenue
      (dist_to_ports$dp3/1000)*cost_per_km*2 + dist_to_ports$suitSpC*bycost             #cost
    plot(rasterFromXYZ(dist_to_ports[,c("lon","lat","utility_p3_BY")]))  #*** JS: plot P1 utility
    # port 4
    dist_to_ports$utility_p4_BY <- price*dist_to_ports$abund_t1 - #revenue
      (dist_to_ports$dp4/1000)*cost_per_km*2 + dist_to_ports$suitSpC*bycost             #cost
    plot(rasterFromXYZ(dist_to_ports[,c("lon","lat","utility_p4_BY")]))  #*** JS: plot P1 utility
    # port 5
    dist_to_ports$utility_p5_BY <- price*dist_to_ports$abund_t1 - #revenue
      (dist_to_ports$dp5/1000)*cost_per_km*2 + dist_to_ports$suitSpC*bycost
    plot(rasterFromXYZ(dist_to_ports[,c("lon","lat","utility_p5_BY")]))  #*** JS: plot P1 utility

    dist_to_ports$utility_max_BY <- apply(dist_to_ports[,c("utility_p1_BY", "utility_p2_BY",
                                                         "utility_p3_BY", "utility_p3_BY",
                                                         "utility_p5_BY")], 1, FUN=max)

    dist_to_ports$utility_SUM_BY <- apply(dist_to_ports[,c("utility_p1_BY", "utility_p2_BY",
                                                         "utility_p3_BY", "utility_p3_BY",
                                                         "utility_p5_BY")], 1, FUN=sum)

    # head(dist_to_ports)
    # # take make of all ports for coast-wide utility
    df_util_raster_BY <- rasterFromXYZ(dist_to_ports[,c("lon","lat","utility_max_BY")])  #*** JS: save coast-wide utility
    plot(df_util_raster_BY, asp=1, main="Coast-wide Utility")

    df_sum_raster_BY<-rasterFromXYZ(dist_to_ports[,c("lon","lat","utility_SUM_BY")])  #*** JS: save coast-wide utility
    plot(df_sum_raster_BY, asp=1, main="Coast-wide Utility")

    # # *** Rescale utility_max between 0 and 1
    df_util_raster2_BY <- df_util_raster_BY + abs(minValue(df_util_raster_BY))
    df_util_raster2_BY <- df_util_raster2_BY/maxValue(df_util_raster2_BY)  ## df_util_raster2 is the raster we use in the next section to sample presence-absences
    plot(df_util_raster2_BY, asp=1, main="Coast-wide Utility")
    points(-117.1441, 32.6717, pch=0, cex=2)
    points(-122.001620, 36.965719, pch=0, cex=2)
    points(-123.050618, 38.334302, pch=0, cex=2)
    points(-124.292000, 43.383975, pch=0, cex=2)
    points(-124.114934, 46.911534, pch=0, cex=2)

    # # *** Rescale utility_SUM between 0 and 1
    df_sum_raster2_BY <- df_sum_raster_BY + abs(minValue(df_sum_raster_BY))
    df_sum_raster2_BY <- df_sum_raster2_BY/maxValue(df_sum_raster2_BY)  ## df_util_raster2_BY is the raster we use in the next section to sample presence-absences
    plot(df_sum_raster2_BY, asp=1, main="Coast-wide Utility")
    points(-117.1441, 32.6717, pch=0, cex=2)
    points(-122.001620, 36.965719, pch=0, cex=2)
    points(-123.050618, 38.334302, pch=0, cex=2)
    points(-124.292000, 43.383975, pch=0, cex=2)
    points(-124.114934, 46.911534, pch=0, cex=2)


   #-----SAMPLE PRESENCES AND ABSENCES-----#####
    
      #******Random Sampling of nsamples*******
    set.seed(y)
    presence.points.random <- sampleOccurrences(suitability_PA,n = nsamples,type = "presence-absence",
                                                detection.probability = 1,error.probability=0, plot = TRUE,
                                                sample.prevalence = 0.5)  
    #convert to dataframe
    pres_df_random <- cbind(as.data.frame(presence.points.random$sample.points$x),as.data.frame(presence.points.random$sample.points$y))
    colnames(pres_df_random) <- c("x","y")
    pres_df_random$random_sampled <- 1
    #expand dataframe to include all possible locations
    df_full <- as.data.frame(rasterToPoints(sst)[,1:2]) #picking an example raster to extract lat and lon from
    df_full_2 <- left_join(df_full, pres_df_random, by=c('x','y'))
    df_full_2$random_sampled <- ifelse(is.na(df_full_2$random_sampled),0,df_full_2$random_sampled)
    df_full_2$lon <- df_full_2$x
    df_full_2$lat <- df_full_2$y
    
    #******Preferential Sampling of nsamples - based on Habitat suitability of target species**********
    set.seed(y)
    presence.points.preferential <-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                                     detection.probability = 1, bias = "manual", weights = spB_suitability_t1$suitab.raster, plot = TRUE)
    #convert to dataframe
    pres_df_pref <- cbind(as.data.frame(presence.points.preferential$sample.points$x),as.data.frame(presence.points.preferential$sample.points$y))
    colnames(pres_df_pref) <- c("x","y")
    pres_df_pref$pref_sampled <- 1
    #add to dataframe
    df_full_3 <- left_join(df_full_2, pres_df_pref, by=c('x','y'))
    df_full_3$pref_sampled <- ifelse(is.na(df_full_3$pref_sampled),0,df_full_3$pref_sampled)

    
    #*******Sampling of nsamples based on RUM (expected catch & Distance to Port)********** 
    set.seed(y)
    presence.points.RUM<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                                    detection.probability = 1, bias = "manual", weights= df_sum_raster2, plot=TRUE)
    #convert to dataframe
    pres_df_RUM <- cbind(as.data.frame(presence.points.RUM$sample.points$x),as.data.frame(presence.points.RUM$sample.points$y))
    colnames(pres_df_RUM) <- c("x","y")
    pres_df_RUM$RUM_sampled <- 1
    #add to dataframe
    df_full_4 <- left_join(df_full_3, pres_df_RUM, by=c('x','y'))
    df_full_4$RUM_sampled <- ifelse(is.na(df_full_4$RUM_sampled),0,df_full_4$RUM_sampled)

   
     #*******Sampling of nsamples based on Habitat Suitability & Distance to Port & Bycatch********** 
     set.seed(y)
     presence.points.BY<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                            detection.probability = 1, bias = "manual", weights= df_sum_raster2_BY , plot=TRUE)
     #convert to dataframe
     pres_df_BY <- cbind(as.data.frame(presence.points.BY$sample.points$x),as.data.frame(presence.points.BY$sample.points$y))
     colnames(pres_df_BY) <- c("x","y")
     pres_df_BY$BY_sampled <- 1
     #add to dataframe
     df_full_5 <- left_join(df_full_4, pres_df_BY, by=c('x','y'))
     df_full_5$BY_sampled <- ifelse(is.na(df_full_5$BY_sampled),0,df_full_5$BY_sampled)

    #----EXTRACT DATA for each year----
    
    print("Extracting suitability")
    ei <- 21912*y #end location in output grid to index to
    se <- ei - (21912-1) #start location in output grid to index to
    output$lat[se:ei] <- df_full_4$y
    output$lon[se:ei] <- df_full_4$x
    output$year[se:ei] <- rep(years[y],21912)
    output$pres_t[se:ei] <- rasterToPoints(suitability_PA$pa.raster)[,3] 
    output$pres_t1[se:ei] <- rasterToPoints(suitability_PA_t1$pa.raster)[,3] 
    output$suitability_t[se:ei] <- rasterToPoints(spB_suitability$suitab.raster)[,3]  #extract points from suitability file
    output$suitability_t1[se:ei] <- rasterToPoints(spB_suitability_t1$suitab.raster)[,3] 
    output$sst[se:ei] <-  rasterToPoints(sst)[,3]   #extract points from suitability file
    output$zoo_200[se:ei] <-  rasterToPoints(zoo)[,3] 
    output$mld[se:ei] <-  rasterToPoints(mld)[,3] 
    #temporary storage: make sure sample coordinates match output coordinates
    temp_random <-   left_join(output[se:ei,1:11], df_full_4[,c(3,4,5)], by=c('lon','lat'))
    temp_pref <-   left_join(output[se:ei,1:11], df_full_4[,c(4,5,6)], by=c('lon','lat'))
    temp_RUM <-   left_join(output[se:ei,1:11], df_full_4[,c(4,5,7)], by=c('lon','lat'))
    #assign to output
    output$random_sampled[se:ei] <- temp_random$random_sampled
    output$pref_sampled[se:ei] <- temp_random$pref_sampled
    output$RUM_sampled[se:ei]<- temp_random$RUM_sampled
  }
  
  #Average monthly biomass available to CCS is: 1.18x10^5 ± (0.13x10^5 se) mt (from Desiree Tommasi)
  mean_spatial <- round(118000/140, 1) 
  se_spatial <- round((13000/140) ,2) 
  output$abundance_t <- ifelse(output$pres_t==1,rnorm(nrow(output),mean_spatial, se_spatial)*output$suitability_t,0)
  
  print('Saving csv to working directory')
  write.csv(output, 'FisheryDependent_OM_Simulation.csv',row.names = FALSE)
  return(output)
}

# test<-SimulateWorld_ROMS_PMP(dir=dir, nsamples = 50)


