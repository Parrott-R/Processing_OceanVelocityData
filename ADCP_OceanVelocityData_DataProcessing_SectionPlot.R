# name: Importing, processing, and plotting of shipboard ADCP data
# version: 1.0.0
# created by: Ruan G. Parrott on 09/12/2019
# contributed to and edited by: Shantelle Smith
# last edited: 13/09/2020

#### Load required packages ####
require(oce)
require(tidyverse)
require(lubridate)
require(leaflet)
require(sf)
require(MBA)
require(scales)
require(marmap)
require(metR)
require(reshape2)
#### Read in x number of files containing ADCP data #####
### Read in ADCP data (multiple ADCP file formats accepted) into list for x number of files
setwd("/path/")

adp <- c(read.adp("file1.LTA"),
         read.adp("file2.LTA"),
         read.adp("file2.LTA"))

#### Process/clean raw data and create modified variables (velocity and distance) ####
# Choose date bounds for each file (can be omitted. see cruise.adcp below)
date = data.frame(as_datetime(c("2018-07-18 21:31:53", "2018-07-21 04:51:53")),
        as_datetime(c("2018-07-18 21:31:53", "2018-07-23 11:00:00")),
        as_datetime(c("2018-07-18 21:31:53", "2018-07-23 11:00:00")))
colnames(date) <- c("1","2","3")

# Initialise variables 
adp.bin <- c()
time <- c()
cruise.adcp <- c()
distance <- c()
lon <- c()
lat <- c()
v <- c()
u <- c()
vel <- c()
cruise.sf <- c()
ship.vel <- c()
vel.abs <- c()
loni <- c()
lati <- c()
abs.vel.profile <- c()
abs.vel.profile.long <- c()
percent_faulty <- c()

# Clean data (remove erroneous velocities) and generate distance/velocity variables
for (adcp in 1:length(adp)) {
  print(paste("Starting to process file",adcp))
  
  ### Bin-map adp object 
  # by interpolating velocities, backscatter amplitudes, etc., to uniform depth bins, 
  # compensating for the pitch and roll of the instrument.
  adp.bin[[adcp]] <- binmapAdp(adp[[adcp]])

  ### Extract variables and subset
  # Create subset of dataset between specific dates 
  time[[adcp]] <- adp[[adcp]][["time"]]
  # Let cruise.adcp = adp.bin if no time subsetting is required
  cruise.adcp[[adcp]] = subset(adp.bin[[adcp]], time >= date[1,adcp] & time <= date[2,adcp])
  
  # Extract subset variables
  time[[adcp]] <- cruise.adcp[[adcp]][["time"]]
  distance[[adcp]] <- cruise.adcp[[adcp]][["distance"]] 
  lon[[adcp]] <- cruise.adcp[[adcp]][["firstLongitude"]]
  lat[[adcp]] <- cruise.adcp[[adcp]][["firstLatitude"]]
  
  # Calculate velocity --> velocity in ms^−1 = √(U^2+V^2)
  v[[adcp]] <- cruise.adcp[[adcp]][["v"]][,,1] # v component of velocity
  u[[adcp]] <- cruise.adcp[[adcp]][["v"]][,,2] # u component of velocity
  vel[[adcp]] <- sqrt(v[[adcp]]^2 + u[[adcp]]^2) 
  
  ### Confirm sampling locations
  # Create simple features object
  cruise.sf[[adcp]] = data.frame(lon[[adcp]],lat[[adcp]],vel[[adcp]])%>%
                        st_as_sf(coords = c("lon..adcp..", "lat..adcp.."))%>% 
                        st_set_crs(4326) # WGS84 reference ellipsoid
  
  # Plot the location of ADCP measurements
  leaflet(data = cruise.sf[[adcp]])%>%
    addTiles()%>%
    addMarkers(popup = ~time[[adcp]])
  
  ### Calculate absolute velocity and compile other variables into new data frame
  # Calculate ship velocity 
  ship.vel[[adcp]] <- cruise.adcp[[adcp]][["avgSpeed"]]
  
  # Obtain absolute current velocity
  vel.abs[[adcp]] <- abs(vel[[adcp]]-ship.vel[[adcp]])
  
  # Create absolute current velocity profiles
  abs.vel.profile[[adcp]] <- data.frame(distance[[adcp]],t(vel.abs[[adcp]]))
  colnames(abs.vel.profile[[adcp]]) <- c("distance",paste("profile", 
                            2:length(abs.vel.profile[[adcp]])-1, sep = ""))
  
  # Make long data frame for absolute current velocity profiles 
  abs.vel.profile.long[[adcp]] <- abs.vel.profile[[adcp]] %>% 
                                gather(key = "profile", value = "velocity", 
                                       2:length(abs.vel.profile[[adcp]]))
  colnames(abs.vel.profile.long[[adcp]])[1] <- "depth"
  
  # Add lats & longs & distance, they repeat in patterns for all profiles
  abs.vel.profile.long[[adcp]]$lon <- rep(NA,length(abs.vel.profile.long[[adcp]]$profile))
  abs.vel.profile.long[[adcp]]$lat <- rep(NA,length(abs.vel.profile.long[[adcp]]$profile))
  abs.vel.profile.long[[adcp]]$distance <- rep(NA,length(abs.vel.profile.long[[adcp]]$profile))
  for (i in 1:length(unique(abs.vel.profile.long[[adcp]]$profile))){
    abs.vel.profile.long[[adcp]]$lon[which(
          abs.vel.profile.long[[adcp]]$profile 
          == paste("profile", i, sep = ""))] <- lon[[adcp]][i]
    abs.vel.profile.long[[adcp]]$lat[which(
          abs.vel.profile.long[[adcp]]$profile 
          == paste("profile", i, sep = ""))] <- lat[[adcp]][i]
  }
  for (j in 1:length(abs.vel.profile.long[[adcp]]$lon)) {
    loni[[adcp]] <- abs.vel.profile.long[[adcp]]$lon[j]
    lati[[adcp]] <- abs.vel.profile.long[[adcp]]$lat[j]
    abs.vel.profile.long[[adcp]]$distance[j] <- geodDist(longitude1 = first(abs.vel.profile.long[[1]]$lon),
                                                 latitude1 = first(abs.vel.profile.long[[1]]$lat),
                                                 longitude2 = loni[[adcp]], 
                                                 latitude2 = lati[[adcp]], 
                                                 alongPath = F)
  }
  
  # Find percentage of data that is faulty (data quality analysis)
  # (>2.5 ms^-1 chosen from literature knowledge --> to be edited to appropriate value for respective current)
  percent_faulty[[adcp]] <- length(abs.vel.profile.long[[adcp]]$velocity[which
              (abs.vel.profile.long[[adcp]]$velocity > 2.5)])/
              length(abs.vel.profile.long[[adcp]]$velocity)*100
  print(paste(round(percent_faulty[[adcp]], digits=1),"% erroneous data", sep=""))
  
  # Set faulty data to NA
  abs.vel.profile.long[[adcp]]$velocity[which(abs.vel.profile.long[[adcp]]$velocity > 2.5)] <- NA

} # END of for loop

#### Compile x number of abs velocity profiles into one data frame #####
# (if data is in multiple files for transect)
cruise.Vel <- rbind(abs.vel.profile.long[[1]],
                  abs.vel.profile.long[[2]],
                  abs.vel.profile.long[[3]])
cruise.Vel.sel <- select(cruise.Vel, c("distance","depth","velocity"))

#### Create interpolated data frame for section plot ####
# Remove NA's and subset for top 1000 m
dat.var <- na.omit(cruise.Vel.sel)
surfdata <- filter(dat.var, depth < 1000) # Optional: subset for surface data

# Interpolation
data_mba <- mba.surf(surfdata, no.X = 300, no.Y = 300,extend = T,n=1,m=1)
dimnames(data_mba$xyz.est$z) <- list(data_mba$xyz.est$x, 
                                     data_mba$xyz.est$y)
data_mba <- melt(data_mba$xyz.est$z, varnames = c('Distance', 'Depth'),
                 value.name = 'velocity')

#### Optional: Create station labels for section plot ####
# Read in table of station positions
# Paste lat/lon table in text="" or read in csv/text file (only containing columns for lat and lon)
cruiseLine <- read.table(header = T, text= "") 

# Calculate distance between stations
cruiseStns <- as.matrix(earth.dist(cruiseLine))[1,]

# Create data frame for station labels for section plot
stationnames <- data.frame(unique(as.character(rownames(cruiseLine))))
latitudenames <- data.frame(as.numeric(as.character(
                                unique(cruiseLine$lat))))
longitudenames <- data.frame(as.numeric(as.character(
                                unique(cruiseLine$long))))
distancenames <- data.frame(as.numeric(as.character(cruiseStns)))

stations <- data.frame(stationnames,
                       latitudenames, longitudenames,
                       c(2,2,2,2,2,2,2,2,2,2, # set depth to 2 m (surface)
                         2,2,2,2,2,2,2,2,2,2), distancenames)
colnames(stations) <- c('station','lat','long','depth', 'distance')

#### Bathymetry for section plot ####
depth.plot <- 400 # change to depth of section plot

# Read in GEBCO file for bathymetry data
setwd("/path/")
Mybathy <- readGEBCO.bathy("gebco_2020.nc") # insert own file

# Subset bathymetry data
NewLatituden <- na.omit(cruise.Vel$lat)
NewLongituden <- na.omit(cruise.Vel$lon)
lati <- seq(min(NewLatituden), max(NewLatituden), 0.01)
loni <- approx(NewLatituden, NewLongituden, lati, rule=2)$y
dist <- rev(geodDist(loni, lati, alongPath=TRUE))
bottom <- data.frame(get.depth(Mybathy, x=loni, y=lati, locator=FALSE), dist)
bottom$depth <- as.numeric(as.character(bottom$depth))

BathyFull_df <- data.frame(x = c(min(dist), last(dist),
                                 rev(dist[2:length(dist)]),
                                 max(dist)),
                           y= c(first(bottom$depth), last(bottom$depth),
                                rev(bottom$depth[2:length(dist)]),
                                first(bottom$depth)))

BathySurf_df <- data.frame(BathyFull_df[["x"]][which(
  (-BathyFull_df[["y"]] - depth.plot) <= 0)], BathyFull_df[["y"]][which(
    (-BathyFull_df[["y"]] - depth.plot) <= 0)])
names(BathySurf_df)[1] <- "x"
names(BathySurf_df)[2] <- "y"

# Create coords for the bathymetry
PolyX <- c(BathySurf_df[["x"]][1],
           BathySurf_df[["x"]][which(BathySurf_df[["y"]] > -depth.plot)],
           last(BathySurf_df[["x"]][which(BathySurf_df[["y"]] > -depth.plot)]))

PolyY <- c(depth.plot, 
           -BathySurf_df[["y"]][which(BathySurf_df[["y"]] > -depth.plot)],
           depth.plot)
rm(Mybathy) # remove large file to save environment memory

#### Section plot: ggplot2 ####
ggplot() +
  scale_y_reverse(limits=c(400,-10), expand = c(0,0))+
  scale_x_continuous(breaks = c(0,50,100,150,200,250), expand = c(0,0))+ 
  geom_tile(data = data_mba, aes(fill = velocity,
                                 x= Distance, y= Depth)) + 
  theme_classic() +
  geom_contour(aes(z = data_mba$velocity,
                   x = data_mba$Distance,
                   y = data_mba$Depth ), 
               binwidth = 0.15, 
               colour = "black", alpha = 0.3) + 
  geom_text_contour(aes(z = data_mba$velocity, 
                        x = data_mba$Distance,
                        y = data_mba$Depth ),
                    check_overlap = TRUE, binwidth = 0.3, alpha = 0.7) +
  scale_fill_gradientn(colours = viridis_pal()(10),
                       values = rescale(c(0,0.3,0.6)), 
                       guide = "colorbar", limits=c(0,2.5),oob=squish,
                       na.value = 'white') + 
  guides(fill = guide_colourbar(barwidth = 0.8, barheight = 15, 
                                nbin = 50,
                                draw.ulim = FALSE, draw.llim = FALSE,
                                label.theme = element_text(size = 15),
                                ticks.colour = "black")) +
  geom_polygon(aes(x = PolyX,
                   y = PolyY),
               fill="black")+
  labs(y = "Depth [m]", x = "Distance from Station 1 [km]", 
       fill = "Velocity [m/s]") + 
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))+
  coord_cartesian(expand = 0)

pathsave = "/path/" # insert path name

ggsave(plot = last_plot(),
       filename = "sample_image.png",
       path = pathsave,
       device = "png", dpi = 1200,
       units = "cm", width = 22, height = 14)
dev.off()