#This code is the operating model for a HMS archetype representing a Pelagic mobile predator-like species
#It uses average spring conditions from downscales ROMS projections using Hadley 
#IMPORTANT: Download average spring ROMS data from here: https://www.dropbox.com/sh/pj1vjz4h1n27hnb/AAD9sySDudMKmq3l4yLPvjQra?dl=0

#We include a trophic interaction between predator (PMP) and prey, but note the estmation model will use chl-a as a proxy for prey information:
#prey: distribution and suitability driven by SST and zooplankton integrated across top 200m
#predator: distribution and abundance driven by SST, MLD, and Species A
#Simulates fishery dependent sampling bias where fishing locations (sample locations) is skewed to be more common
#in areas where fish abundance was high in the year before (y-1 suitability was high)

dir <- "~/DisMAP project/Location, Location, Location/Location Workshop/ROMS"

SimulateWorld_ROMS_PMP_PrefSamp_400 <- function(dir){
  #dir is the local directory that points to where ROMS data is stored 
  
  #required libraries
  library(virtualspecies)
  library(scales)
  library(raster)
  library(glmmfields)
  
  #----Create output file----
  #Assuming 400 'samples' are taken each year, from 1980-2100
  #This will be the information passed to the estimation model
  output <- as.data.frame(matrix(NA, nrow=48400,ncol=9))
  colnames(output) <- c("lon","lat","year","pres","suitability","sst","zoo_200","chla_surface", "mld")
  
  #----Load in rasters----
  #These are the average spring conditions from the downscaled gfdl earth system model
  files_sst <- list.files(paste0(dir,'/sst_spring_avg'), full.names = TRUE, pattern=".grd") #should be 1452 files
  files_chl_surface <- list.files(paste0(dir,'/chl_surface'), full.names = TRUE, pattern=".grd") #should be 1452 files
  files_mld <- list.files(paste0(dir,'/ild_0.5C'), full.names = TRUE, pattern=".grd")
  files_zoo <- list.files(paste0(dir,'/zoo_200m'), full.names = TRUE, pattern=".grd")
  years <- seq(1980,2100,1)
  
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
    #plot(spB_suitability$suitab.raster) #plot habitat suitability
    
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
    #plot(spB_suitability_t1$suitab.raster) #plot habitat suitability
    
    #----CONVERT SP B SUITABILITY TO PRESENCE-ABSENCE----####
    
    #Use a specific function to convert suitability (0-1) to presence or absence (1 or 0)
    
    suitability_PA <- virtualspecies::convertToPA(spB_suitability, PA.method = "probability", beta = 0.5,
                                                  alpha = -0.05, species.prevalence = NULL, plot = FALSE)
    # plotSuitabilityToProba(suitability_PA) #Let's you plot the shape of conversion function
    #plot(suitability_PA$pa.raster)
    
    #y-1 presence absence 
    suitability_PA_t1 <- virtualspecies::convertToPA(spB_suitability_t1, PA.method = "probability", beta = 0.5,
                                                     alpha = -0.05, species.prevalence = NULL, plot = FALSE)
    #plot(suitability_PA_t1$pa.raster)
    
    #-----SAMPLE PRESENCES AND ABSENCES-----#####
    
    #Convert Suitability_PA to a polygon of presence area, and bias sampling intensity to polygon of high suitability for Sp B in y-1

    # pres_poly_t1<- rasterToPolygons(suitability_PA_t1$pa.raster, function(x){x == 1}, dissolve = TRUE) 
    # plot(pres_poly_t1)
    
    presence.points<-sampleOccurrences(suitability_PA, n=100, type="presence-absence",
                                       detection.probability = 1, bias = "manual", weights = suitability_PA_t1$pa.raster, plot = TRUE)
   
    # presence.points<-sampleOccurrences(suitability_PA, n=400, type="presence-absence",
    #                                   detection.probability = 1, bias = "polygon", bias.strength = 50, bias.area= pres_poly_t1, plot = TRUE)
    df <- cbind(as.data.frame(round(presence.points$sample.points$x,1)),as.data.frame(round(presence.points$sample.points$y,1)))
    colnames(df) <- c("x","y")
    
    
    #----EXTRACT DATA for each year----
    
    print("Extracting suitability")
    ei <- 400*y #end location in output grid to index to
    se <- ei - 399 #start location in output grid to index to
    output$lat[se:ei] <- df$y
    output$lon[se:ei] <- df$x
    output$year[se:ei] <- rep(years[y],400)
    output$pres[se:ei] <- presence.points$sample.points$Real
    output$suitability[se:ei] <- raster::extract(spB_suitability$suitab.raster, y= df)  #extract points from suitability file
    output$sst[se:ei] <-  raster::extract(sst, y= df)  #extract points from suitability file
    output$zoo_200[se:ei] <-  raster::extract(zoo, y= df)
    output$chla_surface[se:ei] <-  raster::extract(chla_surface, y= df)
    output$mld[se:ei] <-  raster::extract(mld, y= df)
  }
  
  #Average monthly biomass available to CCS is: 1.18x10^5 ± (0.13x10^5 se) mt (from Desiree Tommasi)
  mean_spatial <- round(118000/140, 1) 
  se_spatial <- round((13000/140) ,2) 
  output$abundance <- ifelse(output$pres==1,rnorm(nrow(output),mean_spatial, se_spatial)*output$suitability,0)
  
  return(output)
}

#test<-SimulateWorld_ROMS_PMP_PrefSamp_400(dir=dir)
