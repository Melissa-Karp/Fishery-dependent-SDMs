## Calculating distance from port for each cell in ROMS raster layer 
library(raster)
library(rgdal)
setwd("~/DisMAP project/Location, Location, Location/Location Workshop/ROMS/sst_spring_avg")

#create output data frame for all the distance to port data
output <- as.data.frame(matrix(NA, nrow=33666,ncol=7))
colnames(output) <- c("lat", "lon", "dp1", "dp2", "dp3", "dp4", "dp5")

p1<- c(-117.1441, 32.6717)#San Diego bay, CA
p2<- c(-122.001620, 36.965719)#Santa Cruz, CA
p3<- c(-123.050618, 38.334302)#Bodega Bay, CA
p4<- c(-124.292000, 43.383975)#Coos Bay, OR
p5<- c(-124.114934, 46.911534)#Westport, WA

#read in ROMS sst data and convert to raster layer, and plot Port locations on Map
ROMS<-raster("sst_monthly_gfdl_SpringMean_1982.grd")
plot(ROMS)
points(-117.1441, 32.6717, pch=0, cex=2)
points(-122.001620, 36.965719, pch=0, cex=2)
points(-123.050618, 38.334302, pch=0, cex=2)
points(-124.292000, 43.383975, pch=0, cex=2)
points(-124.114934, 46.911534, pch=0, cex=2)

## Calculate distance of every cell to each Port 
dp1<- distanceFromPoints(ROMS, p1) ## distances are in meters I think?
dp2<- distanceFromPoints(ROMS, p2)
dp3<- distanceFromPoints(ROMS, p3)
dp4<- distanceFromPoints(ROMS, p4)
dp5<- distanceFromPoints(ROMS, p5)

#Extract distance values for each cell and put in output file
output$lat <- rasterToPoints(dp1)[,2]
output$lon <- rasterToPoints(dp1)[,1]
output$dp1 <- rasterToPoints(dp1)[,3]
output$dp2 <- rasterToPoints(dp2)[,3]
output$dp3 <- rasterToPoints(dp3)[,3]
output$dp4 <- rasterToPoints(dp4)[,3]
output$dp5 <- rasterToPoints (dp5)[,3]
head(output)

write.csv(output, "~/DisMAP project/Location, Location, Location/Location Workshop/Dist_to_Ports.csv")


