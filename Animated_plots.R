## This will make animated gif maps of sampling locations for each bias situation across years
setwd("~/DisMAP project/Location, Location, Location/Location Workshop")
dat<-read.csv('FisheryDependent_OM_Simulation_Final.csv', header=T, sep=",")
dat_historic_1 <- dat[dat$year<=2010,]
dat_historic<-dat_historic_1[dat_historic_1$year>1984,]

library(ggplot2)
library(gganimate)
library(maps)
library(ggthemes)
library(gifski)
library(gapminder)
library(dplyr)
require(viridis)
library(mapdata)
library(cmocean)

### USE THIS CODE ###

#Random Sampling GIF
hms_animate_map_RS <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$random_sampled==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Random Sample {frame_time}") #takes some time 

hms_animated_RS <- gganimate::animate(hms_animate_map_RS,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('RandomSample_2.gif', hms_animated_RS)

#Target Species Optimum Sampling GIF - 0.5
hms_animate_map_PS <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$pref_sampled_1==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Preferential Sampling 0.5 {frame_time}") #takes some time 

hms_animated_PS <- gganimate::animate(hms_animate_map_PS,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('Pref_sample_0.5.gif', hms_animated_PS)

#Target Species Optimum Sampling GIF - 0.6
hms_animate_map_PS <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$pref_sampled_2==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Preferential Sampling 0.6 {frame_time}") #takes some time 

hms_animated_PS <- gganimate::animate(hms_animate_map_PS,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('Pref_sample_0.6.gif', hms_animated_PS)

#Target Species Optimum Sampling GIF - 0.7
hms_animate_map_PS <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$pref_sampled_3==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Preferential Sampling 0.7 {frame_time}") #takes some time 

hms_animated_PS <- gganimate::animate(hms_animate_map_PS,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('Pref_sample_0.7.gif', hms_animated_PS)

#Target Species Optimum Sampling GIF - 0.8
hms_animate_map_PS <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$pref_sampled_4==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Preferential Sampling 0.8 {frame_time}") #takes some time 

hms_animated_PS <- gganimate::animate(hms_animate_map_PS,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('Pref_sample_0.8.gif', hms_animated_PS)

#Target Species Optimum Sampling GIF - 0.9
hms_animate_map_PS <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$pref_sampled_5==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Preferential Sampling 0.9 {frame_time}") #takes some time 

hms_animated_PS <- gganimate::animate(hms_animate_map_PS,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('Pref_sample_0.9.gif', hms_animated_PS)

#Distance From Port - NPO
hms_animate_map_dist <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$dist_sampled_npo==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Distance Sample NPO {frame_time}") #takes some time 

hms_animated_DS <- gganimate::animate(hms_animate_map_dist,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('DistSample_npo.gif', hms_animated_DS)

#Distance NPNearshore
hms_animate_map_dist <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$dist_sampled_npn==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Distance Sample NPN {frame_time}") #takes some time 

hms_animated_DS <- gganimate::animate(hms_animate_map_dist,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('DistSample_npn.gif', hms_animated_DS)

#Distance - MPO
hms_animate_map_dist <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$dist_sampled_mpo==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Distance Sample MPO {frame_time}") #takes some time 

hms_animated_DS <- gganimate::animate(hms_animate_map_dist,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('DistSample_mpo.gif', hms_animated_DS)

#Distance - MPN
hms_animate_map_dist <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$dist_sampled_mpn==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Distance Sample MPN {frame_time}") #takes some time 

hms_animated_DS <- gganimate::animate(hms_animate_map_dist,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('DistSample_mpn.gif', hms_animated_DS)

#Distance - SPO
hms_animate_map_dist <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$dist_sampled_spo==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Distance Sample SPO {frame_time}") #takes some time 

hms_animated_DS <- gganimate::animate(hms_animate_map_dist,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('DistSample_spo.gif', hms_animated_DS)

#Distance - SPN
hms_animate_map_dist <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$dist_sampled_spn==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Distance Sample SPN {frame_time}") #takes some time 

hms_animated_DS <- gganimate::animate(hms_animate_map_dist,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('DistSample_spn.gif', hms_animated_DS)

#Distance - AllO
hms_animate_map_dist <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$dist_sampled_allo==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Distance Sample ALLO {frame_time}") #takes some time 

hms_animated_DS <- gganimate::animate(hms_animate_map_dist,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('DistSample_allo.gif', hms_animated_DS)

#Distance - AllN
hms_animate_map_dist <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$dist_sampled_alln==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Distance Sample AllN {frame_time}") #takes some time 

hms_animated_DS <- gganimate::animate(hms_animate_map_dist,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('DistSample_alln.gif', hms_animated_DS)


#Bycatch Sampling
hms_animate_map_BY <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$BY_sampled==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Bycatch Avoidance Sample {frame_time}") #takes some time 

hms_animated_BS <- gganimate::animate(hms_animate_map_BY,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('BycatchSample_12_7.gif', hms_animated_BS)

#closed Area Sampling -small
hms_animate_map_closed <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$Closed_sampled_1==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Closed Area Small {frame_time}") #takes some time 

hms_animated_Closed<- gganimate::animate(hms_animate_map_closed,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('ClosedSample_small.gif', hms_animated_Closed)

#closed Area Sampling -medium
hms_animate_map_closed <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$Closed_sampled_2==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Closed Area Medium {frame_time}") #takes some time 

hms_animated_Closed<- gganimate::animate(hms_animate_map_closed,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('ClosedSample_medium.gif', hms_animated_Closed)

#closed Area Sampling -Large
hms_animate_map_closed <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_point(data = dat_historic[dat_historic$pres==1,], color="green")+
  geom_point(data = dat_historic[dat_historic$Closed_sampled_3==1,], color="black")+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  #scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Closed Area Large {frame_time}") #takes some time 

hms_animated_Closed<- gganimate::animate(hms_animate_map_closed,nframes = 25, fps=2, renderer = gifski_renderer())#renders in 
anim_save('ClosedSample_large.gif', hms_animated_Closed)

..........############now have it animated sample locations over temperature...
hms_animate_Temp_Ran <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_tile(aes(fill=sst))+
  geom_point(data = dat_historic[dat_historic$random_sampled==1,], color="black")+
  theme_classic() +  labs(y="", x="") +
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = cmocean("haline")(256)) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Random Sampling {frame_time} with SST") #takes some time 

hms_animated_TR <- gganimate::animate(hms_animate_Temp_Ran,nframes = 25, fps=2, renderer = gifski_renderer())#renders in
anim_save('Temp_Ran_sampled.gif', hms_animated_TR)

hms_animate_Temp_Opt <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_tile(aes(fill=sst))+
  geom_point(data = dat_historic[dat_historic$opt_sampled==1,], color="black")+
  theme_classic() +  labs(y="", x="") +
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = cmocean("haline")(256)) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Opt Sampling {frame_time} with SST") #takes some time 

hms_animated_TO <- gganimate::animate(hms_animate_Temp_Opt,nframes = 25, fps=2, renderer = gifski_renderer())#renders in
anim_save('Temp_Opt_sampled.gif', hms_animated_TO)

hms_animate_Temp_Dist <- ggplot(data = dat_historic, aes(x=lon,y=lat))+
  geom_tile(aes(fill=sst))+
  geom_point(data = dat_historic[dat_historic$Dist_sampled==1,], color="black")+
  theme_classic() +  labs(y="", x="") +
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = cmocean("haline")(256)) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HMS Presence: Dist Sampling {frame_time} with SST") #takes some time 

hms_animated_TD <- gganimate::animate(hms_animate_Temp_Dist,nframes = 25, fps=2, renderer = gifski_renderer())#renders in
anim_save('Temp_Dist_sampled.gif', hms_animated_TD)



####STOP HERE - BELOW IS UNCESSARY CODE####

###subsetting 
 pres<- presencePts[, c(1, 2, 3)]
 rand<- randomPts[, c(1, 2, 3)]
 pref<- PrefPts[, c(1, 2, 3)]
 dist<-DistPts[, c(1, 2, 3)]
 BY<- ByPts[, c(1, 2, 3)]
 prefsub<-pref[pref$year<1985,]
 
 Bydist<- ByDistPts[, c(1, 2, 3)]
  world<-map_data("world")
  USA<-map_data("state")
  west_coast <- subset(USA, region %in% c("california", "oregon", "washington"))
  wcm<-ggplot(data = world) + 
    geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "grey") + 
    coord_quickmap(xlim=c(-180,-50),ylim=c(19,75))+
    theme_classic()
  plot(wcm)
  
  
  PresenceMap<- wcm +
    geom_point(data=pres,aes(x = lon, y = lat), size = 0.5, color='green')+
    geom_point(data=pref,aes(x = lon, y = lat), size = 0.8, color='black')+
    #facet_wrap(~year)
    labs(title = 'Year: {frame_time}') +
    transition_time(year)
    ease_aes('linear')

print(PresenceMap)
Anim_presMap<-animate(PresenceMap, nframes = 121, fps=3)
anim_save("PresenceMap.gif", Anim_presMap)

## Using Steph's code
#abundance
hms_animate_map <- ggplot(data = data, aes(x=lon,y=lat))+
  geom_tile(aes(fill=abundance_t))+
  theme_classic() +  labs(y="", x="") +
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = cmocean("haline")(256),limits = c(0, max(data$pres))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  # scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  # transition_time(year)+
  # ease_aes("linear") +
  # labs(title="HAD: HMS Biomass {frame_time}") #takes some time
hms_animated <- gganimate::animate(hms_animate_map,nframes = 116, fps=3, renderer = gifski_renderer())#renders in
anim_save('AnimatedMap_HMS_abundance.gif', hms_animated)

#Presence-Absence & sampling locations
hms_animate_PA_map <- ggplot(data = data, aes(x=lon,y=lat))+
  geom_point(data=pres,aes(x = lon, y = lat), size = 1.0, color='green')+
  geom_point(data=pref,aes(x = lon, y = lat), size = 1.0, color='black')+
  theme_classic() +  labs(y="", x="") +
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-136,-115.8),ylim=c(28,49)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  transition_time(year)+
  ease_aes("linear") +
  labs(title="HAD: HMS Biomass {frame_time}") #takes some time
hms_animated <- gganimate::animate(hms_animate_PA_map,nframes = 116, fps=3, renderer = gifski_renderer())#renders in
anim_save('AnimatedMap_HMS_PA_2.gif', hms_animated)
