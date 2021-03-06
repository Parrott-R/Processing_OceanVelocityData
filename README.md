# Importing, processing, and plotting of shipboard ADCP data
[![DOI](https://zenodo.org/badge/295239289.svg)](https://zenodo.org/badge/latestdoi/295239289)

This code can be applied to any dataset of ADCP (Acoustic Doppler Current Profiler) data from shipboard instrumentation. There is capacity for the importing of 'x' number of input files. The ADCP data is used to create velocity and distance (from first sampling point; used as a substitute for latitude/longitude - it is ideal for a ship track that does not follow a longitude or latitude). The absolute velocity is calculated from the relative current velocity and the ship's average speed. 

We use the ggplot2 package, supplemented by the scales and metR packages, to create a section plot of absolute current velocity. 

## Getting Started

These instructions will get you a copy of the project up and running.

### Prerequisites

You will need to install the required packages as follows:

```
install.packages(oce)
install.packages(tidyverse)
install.packages(lubridate)
install.packages(leaflet)
install.packages(sf)
install.packages(MBA)
install.packages(scales)
install.packages(marmap)
install.packages(metR)
install.packages(reshape2)
```
The installed packages can then be loaded using library() e.g. library(oce). 

### Importing data

read.adp is used to import the ADCP data and can accept files with extensions used by most instruments. Some pre-processing of files may be required if read.adp is not able to read the file type. 

You may install multiple files at once using the following:

```
adp <- c(read.adp("file1.LTA"),
         read.adp("file2.LTA"),
         read.adp("file3.LTA"))
```
Larger ADCP data files may be pre-processed using adpEnsembleAverage() to reduce the number of profiles and to smooth the data. 

## Processing data

The variables for the processing step are first initialised. 

```
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
```
The data processing is done in a for loop to allow for the input of multiple ADCP files. For this reason, it is strongly recommended to use the ensemble averaging method mentioned above if the files are very large. 

The data is first bin-mapped by interpolating to uniform depths to compensate for the pitch and roll of the ship. There is an option to subset the data between specific dates/times if necessary. To confirm the locations of the sampling, the subsetted data is converted into a simple feature object and visualised using leaflet(). 

Velocity is calculated from the meridional and zonal components using: velocity (in ms^−1) = sqrt(u^2+v^2) and is then used to calculate aboslute velocity by subtracting the ship's average speed. The structure of the data frame is then modified to obtain depth, profile ID, longitude, latitude, distance (calculated to be distance from station 1 using geodDist()), velocity into a single data frame for plotting later on. Spurious data are removed from this data frame using an upper limit (e.g. 2.5 ms^-1) as a cut-off for reasonable values and by changing unreasonable values to 'NA'. A percentage of spurious data is calculated and is useful for assessing general data quality e.g. in the dataset used for testing, all data above 35 m was discarded, likely due to interference from the ship. 

```
for (adcp in 1:length(adp)) {
  print(paste("Starting to process file",adcp))
  
  ### Bin-map an ADP object 
  # by interpolating velocities, backscatter amplitudes, etc., to uniform depth bins, 
  # thus compensating for the pitch and roll of the instrument.
  adp.bin[[adcp]] <- binmapAdp(adp[[adcp]])

  ### Extract variables and subset
  time[[adcp]] <- adp[[adcp]][["time"]]
  
  # Create subset of dataset between specific dates 
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
  print(percent_faulty[[adcp]])
  
  # Set faulty data to NA
  abs.vel.profile.long[[adcp]]$velocity[which(abs.vel.profile.long[[adcp]]$velocity > 2.5)] <- NA

} # END of for loop
```

## Compiling data frames 

The data frames for each file, that are produced above, are combined and subsetted for plotting. The order of the columns in the subsetted data frame are important for the interpolation step below since the function mba.surf() uses data from the three columns in a pre-determined way, selecting each one for a specific dimension of the interpolation (i.e. x, y, z).

```
cruise.Vel <- rbind(abs.vel.profile.long[[1]],
                  abs.vel.profile.long[[2]],
                  abs.vel.profile.long[[3]])
cruise.Vel.sel <- select(cruise.Vel, c("distance","depth","velocity"))
```

## Interpolating the absolute current velocity using the MBA package

Prior to interpolating the data, 'NA' values are removed and, as an optional extra step, the data is subsetted to include only surface data. It's advised to minimise the size of the data frame, where possible, to maximise the efficiency of mba.surf(), which is used to perform the interpolation. The interpolated data is then melted into a long data frame for plotting. 

```
dat.var <- na.omit(cruise.Vel.sel)
surfdata <- filter(dat.var, depth < 1000) 

data_mba <- mba.surf(surfdata, no.X = 300, no.Y = 300,extend = T,n=1,m=1)
dimnames(data_mba$xyz.est$z) <- list(data_mba$xyz.est$x, 
                                     data_mba$xyz.est$y)
data_mba <- melt(data_mba$xyz.est$z, varnames = c('Distance', 'Depth'),
                 value.name = 'velocity')
```

## GEBCO bathymetry - import and subsetting data

The bathymetry data overlaid on the section plot is prepared as shown in the following block. Bathymetry data are obtained from GEBCO, of which smaller files can be created at [https://download.gebco.net/](https://download.gebco.net/) for a specific area defined by latitude and longtiude. 
The bathymetry data is subsetted to the exact latitude and longitude range in the ADCP data and to the top 1000 m. The depth.plot variable needs to match the depth chosen for the section plot below, otherwise the polygon will not fill in correctly. The last step is to create variables PolyX and PolyY for the x-axis and y-axis points of the bathymetry polygon, respectively, setting the beginning and end points of the polygon as the edges of the plot. 

```
depth.plot <- 400 # change to depth of section plot

# Read in GEBCO file for bathymetry data
setwd("/path/")
Mybathy <- readGEBCO.bathy("gebco_2020.nc")

# Subset bathymetry data
NewLatituden <- na.omit(cruise.Vel$lat)
NewLongituden <- na.omit(cruise.Vel$lon)
lati <- seq(min(NewLatituden), max(NewLatituden), 0.01)
loni <- approx(NewLatituden, NewLongituden, lati, rule=2)$y
dist <- rev(geodDist(loni, lati, alongPath=TRUE))
bottom <- get.depth(Mybathy, x=loni, y=lati, locator=FALSE)
bottom <- data.frame(bottom, dist)
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
rm(Mybathy)
```

## Plotting the absolute current velocity using the ggplot2 package

We created a template for a ggplot2 section plot, with a classic theme, to plot the interpolated current velocity data. We recommend ggplot2 over other R plotting packages due to the customisation possibilities. On the x-axis, distance can be replaced with longitude or latitude depending on the layout of the sampling transect. This template uses the colour blind-friendly viridis palette with a customised legend. geom_polygon is used to overlay the bathymetry. We refer the reader to the optional section in the R code if other station/locations need to be plotted as overlaid points on the section plot. 

```
ggplot() +
  scale_y_reverse(limits=c(400,-10), expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+ 
  geom_tile(data = data_mba, aes(fill = velocity,
                                 x= Distance, y= Depth)) +
  theme_classic() +
  geom_contour(aes(z = data_mba$velocity,
                   x = data_mba$Distance,
                   y = data_mba$Depth ), 
               binwidth = 0.15, 
               colour = "black", alpha = 0.3) + #adds contours
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
```
![Absolute velocity in the current.](/sample_image.png)

## Built With

* [oce](https://cran.r-project.org/web/packages/oce/index.html) - Used to import and analyse ADCP data
* [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html) - Used to rearrange data frames
* [leaflet](https://cran.r-project.org/web/packages/leaflet/index.html) - Used to plot sampling locations
* [sf](https://cran.r-project.org/web/packages/sf/index.html) - Used to create simple feature object
* [MBA](https://cran.r-project.org/web/packages/MBA/index.html) - Used to interpolate spatially
* [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html) - Used to melt wide data frame
* [marmap](https://cran.r-project.org/web/packages/marmap/index.html) - Used to import GEBCO bathymetry data
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) - Used for plotting 
* [scales](https://cran.r-project.org/web/packages/scales/index.html) - Used for plotting 
* [metR](https://cran.r-project.org/web/packages/metR/index.html) - Used with ggplot for section plot

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/Parrott-R/Processing_OceanVelocityData/tags). 

## Authors

The authors below contributed equally.
* **Ruan G. Parrott** - *Initial work, research, and plotting* [https://parrott-r.github.io/Me/](https://parrott-r.github.io/Me/)
* **Shantelle Smith** - *Editing, restructuring, and plotting. README file compilation.* [https://shantellesmith.github.io/](https://shantellesmith.github.io/)

See also the list of [contributors](https://github.com/Parrott-R/Processing_OceanVelocityData/graphs/contributors) who participated in this project.

## Acknowledgments

* We would like to acknowledge Masumbuko Semba, whose [blog post](https://semba-blog.netlify.app/10/15/2018/processing-adcp-data-with-r/) was used to assist in the importing of the raw data and in the calculating of the absolute current velocity 
