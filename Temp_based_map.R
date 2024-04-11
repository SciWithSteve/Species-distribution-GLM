setwd("C:/Users/steph/Documents/PhD/Experiments/Temperature development and survival/Temp based map")

dev.temp = read.csv("data/DevelopmentSummary.csv")
dev.temp

plot(Incubation_time ~ Temperature, ylab = "development time (days)", xlab = "deg C", data=dev.temp)

dev.temp$rate = 1 / dev.temp$Incubation_time # get developmental rate (proportion per day)
dev.temp$lnrate = log(dev.temp$rate) # get log of developmental rate
dev.temp$invK = 1 / (dev.temp$Temperature + 273.15) # get inverse of temperature in Kelvin

dev.temp

# plot raw rate vs. temperature
plot(rate ~ Temperature, ylab = "development rate 1/d", xlab = "deg C", data=dev.temp)

plot(lnrate ~ invK, ylab = "ln(development rate)", xlab = "1/K", data=dev.temp)

model = lm(dev.temp$lnrate ~ dev.temp$invK) # fit the model
int=coef(model)[1] # get the estimated intercept
slp=coef(model)[2] # get the estimated slope
plot(lnrate ~ invK, ylab = "ln(development rate)", xlab = "1/K", data=dev.temp) # plot log rate vs. 1/temperature (in Kelvin)
points(dev.temp$invK, slp * dev.temp$invK + int ,type='l', col=1) # plot the prediction

T_A = slp * -1 # Arrhenius temperature
T_REF = 22 # reference temperature
kdot_ref = subset(dev.temp, Temperature == T_REF)$rate # reference rate
kdot_ref

# definition of the 1-parameter Arrhenius function
ArrFunc1 <- function(x, T_A, T_REF, kdot_ref){
  exp(T_A * (1 / (273.15 + T_REF) - 1 / (273.15 + x))) * kdot_ref
}

# plot the data
plot(rate ~ Temperature,  data=dev.temp, ylab = "development rate 1/d", xlab = "deg C")

curve(ArrFunc1(x, T_A, T_REF, kdot_ref),add=TRUE)




# definition of the 5-parameter Arrhenius function
ArrFunc5 <- function(x, T_A, T_AL, T_AH, T_L, T_H, T_REF, kdot_ref){
  exp(T_A * (1 / (273.15 + T_REF) - 1 / (273.15 + x))) / 
    (1 + 
       exp(T_AL * (1 /(273.15 + x) - 1 / T_L)) + 
       exp(T_AH * (1 / T_H - 1 / (273.15 + x)))
    ) * kdot_ref
}

# initial guesses for the other 4 parameters (already have an estimate of T_A)
T_AL = 30000 # Arrhenius temperature at lower temperature threshold
T_AH = 40000 # Arrhenius temperature at upper temperature threshold
T_L = 273.15 + 14 # lower temperature threshold
T_H = 273.15 + 38 # upper temperature threshold

# fit to data
library(minpack.lm)
x = dev.temp$Temperature
y = dev.temp$rate
fit <- nlsLM(y ~ ArrFunc5(x, T_A, T_AL, T_AH, T_L, T_H, T_REF, kdot_ref), 
           start=list(T_A = T_A, T_AL = T_AL, T_AH = T_AH, T_L = T_L, T_H = T_H))

# retrieve the estimated coefficients
coeffs<-coef(fit)
T_A = coeffs[1] # Arrhenius temperture
T_AL = coeffs[2] # Arrhenius temperature at lower temperature threshold
T_AH = coeffs[3] # Arrhenius temperature at upper temperature threshold
T_L = coeffs[4] # lower temperature threshold
T_H = coeffs[5] # upper temperature threshold 

# plot the data
plot(rate ~ Temperature, ylab = "development rate 1/d", xlab = "deg C", data=dev.temp,  ylim=c(0,0.2))

# plot the fitted 1-parameter model
curve(ArrFunc1(x, T_A, T_REF, kdot_ref),add=TRUE)

# plot the fitted 5-parameter model
curve(ArrFunc5(x, T_A, T_AL, T_AH, T_L, T_H, T_REF, kdot_ref), col = "blue", add=TRUE)



###Soil data
library(raster)
library(ncdf4)

soildata <- terra::rast("data/soil_2.5cm_0shade/soil2.5cm_0pctShade_2017.nc") #takes a few mins
soildata=soildata/10
plot(soildata[[1]], main="Soildata 1")
plot(soildata[[13]], main="Soildata 13")
plot(soildata)

#soil.folder = "data/environmental/" # specify the path to the soil files
#soil.files = list.files(soil.folder) # get a list all the files in that folder
#head(soil.files) # have a look at the first 5 values of soil.files

#soil.order = substr(soil.files, 16,18) # use the substr function to subsample the text from positions 16 to 18
#soil.order = gsub(pattern = 'n', x = soil.order, replacement = '')
#soil.order = gsub(pattern = '[.]', x = soil.order, replacement = '')

#soil.files.path = paste0(soil.folder,soil.files[1]) # paste the soil folder to the first file in our list
#soil = brick(soil.files.path) # read the 365 soil layers into 1 object called 'soil'

#soil = soil / 10

#plot(soil[[1]]) # plot the 1st layer in 'soil'
#plot(soil[[13]]) # plot the 13th layer in the soil

lon.lat = cbind(138,-32) # define the longitude and latitude
plot(soildata[[1]]) # plot the midnight soil temperature again
points(lon.lat, cex=1.5, pch=16, col='red') # note that cex specifies the size of the point and 'col' the colour (the default is black).

temp.hr = extract(soildata, lon.lat) # extract data for all layers in 'soil' at location 'lon.lat'
temp.hr = t(temp.hr) # transpose the values, in preparation for the plot command
temp.hr
hrs = seq(1,8760) # create a sequence of hours
plot(temp.hr ~ hrs, type = 'b') # plot temperature as a function of hour, as a line and point graph (using type = 'b')


develop.rate = ArrFunc5( # apply the 5-parameter Arrhenius temperature curve function to the extracted temperatures
  x = temp.hr, T_A = T_A, T_AL = T_AL, T_AH = T_AH, T_L = T_L, T_H = T_H, T_REF = T_REF, kdot_ref = kdot_ref
)
plot(develop.rate ~ hrs, type='b') # plot the results

development = cumsum(develop.rate) # do a cumulative sum of development
plot(development ~ hrs, type='b') # plot cumulative sum of development

develop.rate = develop.rate / 24 # convert to rate per hour
plot(develop.rate ~ hrs, type='b') # plot development rate per hour

development = cumsum(develop.rate) # do a cumulative sum of development
plot(development ~ hrs, type='b') # plot cumulative sum of development

#all.soil = stack(paste0(soil.folder,soil.files[1:365]))
#all.soil = all.soil / 10 # divide by 10, because the original grids were rounded to the nearest degree and multiplied by 10 to reduce file size





lon.lat = cbind(142.14, -34.20) # define the longitude and latitude
temp.hr = extract(soildata, lon.lat) # extract the data for this site, across all layers (24 hrs x 1 year)
temp.hr = t(temp.hr) # transpose the result from row to column

tzone = paste("Etc/GMT-",10,sep = "") # specify a time zone
dates = seq(ISOdate(2017,1,1,tz = tzone) - 3600 * 12, ISOdate(2018,1,1,tz = tzone) - 3600 * 13, by = "hours") # create a sequence of dates/times in hourly steps for all of 1990

library("clock")
start <- date_time_parse("2017-01-01 00:00:00", zone = "Australia/Sydney")
end <- date_time_parse(
  c("1990-01-01 00:00:00", "2018-01-01 00:00:00"),
  zone = "Australia/Sydney"
) #set end ()
date_count_between(start, end, "hour")
#8760 hours in a year

plot(temp.hr ~ dates, type = 'l') # plot the soil temperatures at 2.5 cm for this site




develop.rate = ArrFunc5( # compute the development rates  (rate per day)
  x = temp.hr, T_A = T_A, T_AL = T_AL, T_AH = T_AH, T_L = T_L, T_H = T_H, T_REF = T_REF, kdot_ref = kdot_ref
)

develop.rate = develop.rate / 24 # convert to development rate per hour

plot(develop.rate ~ dates, type = 'l') # plot the result



development = cumsum(develop.rate) # do a cumulative sum from the start of the dataset (1st Jan)

plot(development ~ dates, type = 'l',
     xlab = "Month",
     ylab = "Number of generations completed",
     ylim= c(0,9))

abline(1,0, col = 'red', lty=2) # show where development would be complete by drawing a line with the 'abline' function



ggplot()+
  geom_path(x=dates, y=development)



develop.rate = ArrFunc5( # compute the development rates  (rate per day)
  x = soildata, T_A = T_A, T_AL = T_AL, T_AH = T_AH, T_L = T_L, T_H = T_H, T_REF = T_REF, kdot_ref = kdot_ref
) #takes a few mins
develop.rate_perhour = develop.rate / 24 # convert to development rate per hour
plot(develop.rate_perhour[[1]])# plot development rate specified day of year





development = cumsum(develop.rate_perhour) # get cumulative sum of development across all rasters, i.e. each hour from midnight 1st June to 11pm 30th September 1990.
last.day.of.month=c(90,181,273,365) # vector of last day of each month in a sequence of days starting from the 1st of June
last.hr.of.month=last.day.of.month*24 # multiplying the previous vector by 24 to get the position in terms of hours of 11pm on the last day of each quarter
cum.month = development[[last.hr.of.month]] # make a new variable 'cum.month' which just contains the four raster layers in the development brick at the last hour of the last day of each quarter
plot(cum.month, main = c("Jan-March","April-June","Jul-Sep","Oct-Dec")) # plot these four raster layers and label them accordingly 
cum.year = development[[8760]]
plot(cum.year,  main="Yearly development (number cycles)")

cum.year
cumyear_df <- as.data.frame(cum.year, xy = TRUE)
cumyear_df
library(terra)
cumyearrast<-rast(cumyear_df)
cumyearrast2<-raster(cumyearrast)
cumyearrast
crs(cumyearrast) <- '+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs' #define coordinates reference system

library(ggplot2)
library(tidyterra)
library(raster)
ggplot()+
  geom_spatraster(data = cumyearrast)





library("sf")
Contour5<-rasterToContour(cumyearrast2, levels = 5) #extract contour line at #generations=5
Contour5
crs(Contour5) <- '+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs' #define coordinates reference system
plot(cum.year,  main="Yearly development (number cycles)")
plot(Contour5, add=TRUE, col="red", show.legend=FALSE) #add line to map of aus
Contour5 <- raster::disaggregate(Contour5)
Contour5 <- st_as_sf(Contour5)
Contour5
Contour5rast<-rast(Contour5) #now prepared for use in ggplot (PA model projection)

Contour4<-rasterToContour(cumyearrast2, levels = 4) #extract contour line at #generations=5
Contour4
crs(Contour4) <- '+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs' #define coordinates reference system
plot(cum.year,  main="Yearly development (number cycles)")
plot(Contour4, add=TRUE, col="red", show.legend=FALSE) #add line to map of aus
Contour4 <- raster::disaggregate(Contour4)
Contour4 <- st_as_sf(Contour4)
Contour4
Contour4rast<-rast(Contour4) #now prepared for use in ggplot (PA model projection)



###Mask if required
#cum.yearcopy<-cum.year
#cum.yearcopy[cum.yearcopy < 5] <- NA
#cum.yearcopy[cum.yearcopy > 20] <- NA
#plot(cum.yearcopy, main="cum.yearcopy")
#masked <- mask(cum.year, cum.yearcopy)
#plot(masked, main="masked")
#plot(aust_bound, add=TRUE)

cum.month[cum.month>5]<-5 # replace values greater than 10 in 'cum.month' with a value of 5
plot(cum.month, main = c("Jan-March","April-June","Jul-Sep","Oct-Dec")) # replot the results

decmaxtemp = raster("data/BOM/maxdec.txt")
plot(decmaxtemp) #try mask at >36?
julmintemp = raster("data/BOM/minjul.txt")
plot(julmintemp) #try mask at <3 degrees, maybe <2

decmaxtempmask<-decmaxtemp #create duplicate
decmaxtempmasklim<-decmaxtempmask>38
plot(decmaxtempmasklim, main="decmax")
maxmask<-terra::rast(decmaxtempmasklim)
plot(maxmask)
maxmask<-filter(maxmask, values(maxmask)==TRUE)
maxmaska<-aggregate(maxmask, 30)
plot(maxmask, col="red", add=TRUE)
maxmask.sub <- drop_na(maxmask)
maxmaskpoly<-as.polygons(maxmask.sub)
plot(maxmaskpoly)
maxmaskpoly

julmintempmask<-julmintemp #create duplicate
julmintempmasklim<-julmintempmask<2
plot(julmintempmasklim, main="julmin")
minmask<-terra::rast(julmintempmasklim)
plot(minmask)
minmask<-filter(minmask, values(minmask)==TRUE)
minmaska<-aggregate(minmask, 10)
minmaska<-filter(minmaska, values(minmaska)==TRUE)
plot(minmask, add=TRUE, col="blue")
minmask.sub <- drop_na(minmask)
minmaskpoly<-as.polygons(minmask.sub)
plot(minmaskpoly)
minmaskpoly




###FINAL MAP###
library(ozmaps)
plot(cum.year)
ozmap(add=TRUE)
plot(maxmask, col="red", add=TRUE, legend=FALSE)
plot(minmask, add=TRUE, col="blue", legend=FALSE)
plot(Contour5, add=TRUE, col="white", lwd=2, show.legend=FALSE) #add line to map of aus
crs(Contour5) <- '+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs' #define coordinates reference system

AllTrunc<-read.csv("Trunc_presence_only.csv")
head(AllTrunc)

basemap_vect <- ozmap_states %>% 
  vect()
library(ggspatial)

FinalMap<-ggplot()+
  geom_spatraster(data = cumyearrast)+
  geom_sf(data = basemap_vect, # vic state lines on top of raster
                                                fill = NA,
                                                col = "white",
                                                lwd = 0.5)+
  geom_spatvector(data=minmaskpoly, color="black", fill="blue")+

  geom_sf(data = Contour5, col="white", linetype=2)+
  scale_fill_viridis_c(option = "turbo", limits = c(0, 15),
                       breaks = c(0, 5, 10, 15, 20),
                       labels = c(0, 5, 10, 15, 20))+

  coord_sf(xlim = c(114, 153), ylim = c(-44.5, -12)) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Carpophilus truncatus development cycles") +
  labs(data=cumyearrast, fill = "Number of\ngenerations\nper year")

FinalMap



ggsave("Final_map_tidyterra.pdf", width=10, height=6)

  





###Working space###

library(maptools)
P4S <- CRS("+proj=longlat +datum=WGS84")
aust_bound<-readShapeLines("data/borders/ausborder_polyline.shp", verbose=TRUE, proj4string=P4S)
state_bound<-readShapeLines("data/borders/state_boundaries.shp", verbose=TRUE, proj4string=P4S)

fun <- function() {
  plot(aust_bound, col="black", lwd=1.0,add=TRUE)  
  plot(state_bound, col="black", lwd=1.0,add=TRUE)
}
plot(cum.month, main = c("June 30th","July 31st","August 31st","September 30th"), addfun = fun) # replot the results
