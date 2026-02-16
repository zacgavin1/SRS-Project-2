##### Import Data and Load Necessary Packages ----------------------------------
library(ncdf4)
library(lubridate)
library(RColorBrewer)
library(lattice)
library(fields)
library(terra)
library(data.table)

# Add an extra line for your personal file path and comment out as appropriate

setwd("C:\\Users\\Luke Egan\\OneDrive\\Desktop\\Statistical Research Skills\\Assignment 2")
# setwd("C:/Users/zgavi/Documents/Edinburgh Term 2/SRS")
# unzip("temp_data_SRS.nc.zip")
nc <- nc_open("temp_data_SRS.nc")


###### Initial Look at data ----------------------------------------------------
print(nc)
names(nc$dim)
names(nc$var)


# Primary Variable is temperature - measured in degrees celcius
temp <- ncvar_get(nc, "tmp")
temp_units <- ncatt_get(nc, "tmp", "units")
temp_units


# Extract spatial and temporal data
lon <- nc$dim[[1]]$vals
lat <- nc$dim[[2]]$vals
time <- nc$dim[[3]]$vals


# Dimensions of spatial and temporal variables
dim(lon)
dim(lat)
dim(time)


# Check temporal units
t_units <- ncatt_get(nc, "time", "units")
t_units
c(time[1]/365, time[1344]/365)
# ~112 years of data, however time step are not uniform (as month length varies)


# Change time format to DD/MM/YYYY
t_ustr <- strsplit(t_units$value, " ")
t_dstr <- strsplit(unlist(t_ustr)[3], "-")
date <- ymd(t_dstr) +ddays(time)

# c(date[1], date[1344])
# Data ranges from 16/01/1901 to 16/12/2012


# -----------------------------------------------------
# putting the data into a regular df in tidy form - don't use the 
# data in this form for a while

#### NOTE::: This doesn't work at the moment
rast_temp <- rast(
  aperm(temp, c(2, 1, 3)),
  extent=c(xmin = min(lon),
           xmax = max(lon),
           ymin = min(lat),
           ymax = max(lat)),
  crs = "EPSG:4326"
)
names(rast_temp) <- as.character(date)

# making weird dataframe
df <- terra::as.data.frame(rast_temp, xy = TRUE)

# Manually reshape 
coords <- df[, 1:2]
vals   <- df[, -c(1,2)]

# Seasonal Factor
month <- as.POSIXlt(date)$mon + 1

seasons <- factor((month %% 12) %/% 3 + 1,
                  labels = c("winter", "spring", "summer", "autumn"))

# EDIT LUKE: Changed "date" to maintain "Date" class/format
# EDIT LUKE: Inserted season factor column
df_temp <- data.frame(
         lon = rep(coords$x, times = ncol(vals)),
         lat = rep(coords$y, times = ncol(vals)),
         date = rep(date, each = nrow(vals)),
         # date = rep(names(vals), each = nrow(vals)),
         seasons = rep(seasons, each = nrow(vals)),
         temperature = as.vector(as.matrix(vals))
)
# df_temp is now in tidy form
dim(df_temp)
names(df_temp)
head(df_temp)

# -------------------------------------------------------
# Looking further at tidy data and putting seasonal factor column in
# Tried to put seasonal factor column in here but was very inefficient
# Should consider "data.table" R package. See link below for details
# https://r-datatable.com/

# Also consider "dplyr" package for dealing with tidy data
# -------------------------------------------------------
library(rlang)
library(dplyr)

# TO DO: Check if lat and lon is correct in tidy data format
range(lat)
range(df_temp$lat)




# -------------------------------------------------------
# Compare heatmaps same day 111 years apart
par(mfrow = c(2,1), mar = c(3,3,2,5))
image.plot(temp[,,1], col = rev(brewer.pal(10,"RdBu")), main = "Global Temp 01/1901")
image.plot(temp[,,1344], col = rev(brewer.pal(10,"RdBu")), main ="Global Temp 01/2012")



# Difference in temperature 111 years apart 
dev.off()
image.plot(temp[,,1344] - temp[,,1], col = rev(brewer.pal(10,"RdBu")), main ="Temp differences 111 years apart")


# Split data by seasons
ii_winter <- grep("-12-|-01-|-02-", date)
ii_spring <- grep("-03-|-04-|-05-", date)
ii_summer <- grep("-06-|-07-|-08-", date)
ii_autumn <- grep("-09-|-10-|-11-", date)


ii_warm <- sort(c(ii_spring, ii_summer))
ii_cold <- sort(c(ii_autumn, ii_winter))


# Consider one region for EDA - Ireland 
# Longitide (-10, 40), Latitude (35, 72)
ireland_lon_ii <- which(lon > -10 & lon < -5)
ireland_lat_ii <- which(lat > 50 & lat < 55)

temp_ire <- temp[ireland_lon_ii, ireland_lat_ii,]


par(mfrow=c(1,2))
boxplot(as.vector(temp_ire[,, ii_warm[1:112]]),
        as.vector(temp_ire[,, ii_warm[113:224]]),
        as.vector(temp_ire[,, ii_warm[225:336]]),
        as.vector(temp_ire[,, ii_warm[337:448]]),
        as.vector(temp_ire[,, ii_warm[449:560]]),
        as.vector(temp_ire[,, ii_warm[561:672]]))

boxplot(as.vector(temp_ire[,, ii_cold[1:112]]),
        as.vector(temp_ire[,, ii_cold[113:224]]),
        as.vector(temp_ire[,, ii_cold[225:336]]),
        as.vector(temp_ire[,, ii_cold[337:448]]),
        as.vector(temp_ire[,, ii_cold[449:560]]),
        as.vector(temp_ire[,, ii_cold[561:672]]))




############################################################
###### ----- Maximum Yearly Temperatures -------- ##########
############################################################

# Look at yearly maxima in Ireland - global maxima don't really mean that much 

# first just take maxima over each observation time
obs_max_ire <- apply(temp_ire, 3, max, na.rm=T) 

# get the year of each set of observations
year <- as.integer(gsub("-[0-9][0-9]-[0-9][0-9]","", date)) 
years <- unique(year)

yearly_max_ire <- rep(NA, length(years))
for (i in years){
  yearly_max_ire[i-1900] <- max(obs_max_ire[year==i])
}

yearly_max_ire
plot(years, yearly_max_ire)
# unfortunately this plot really doesn't show anything useful in terms of trend

# try another small area
eq_lon_ii <- which(lon > -100 & lon < -95)
eq_lat_ii <- which(lat > 40 & lat < 45)
temp_eq <- temp[eq_lon_ii, eq_lat_ii,]

obs_max_eq <- apply(temp_eq, 3, max, na.rm=T) 
yearly_max_eq <- rep(NA, length(years))
for (i in years){
  yearly_max_eq[i-1900] <- max(obs_max_eq[year==i])
}
yearly_max_eq
plot(years, yearly_max_eq, ylab="yearly max temp observation")

# notes: for -5<lat<0, 0<lon<10 have weird set of equal max and min from 1920-40

# For some regions there is some increase, but overall there is no real 
# upward trend in maxima here




############################################################
######## ------- Overall trend over time -------- ##########
############################################################

# we'll just take the global yearly mean, and plot it over time

obs_mean <- apply(temp, 3, mean, na.rm=T) # this takes a while
yearly_mean <- rep(NA, length(years))
for (i in years){
  yearly_mean[i-1900] <- max(obs_mean[year==i])
}


plot(years, yearly_mean - mean(yearly_mean[1:20]))
fit <- smooth.spline(years, yearly_mean - mean(yearly_mean[1:20]))
lines(fit, col="purple", lwd="2")
# wait, climate change is real??


##### 
# Does this differ for the northern and southern hemispheres?
#####
obs_meann <- apply(temp, 3, mean, na.rm=T) # this takes a while
yearly_mean <- rep(NA, length(years))
for (i in years){
  yearly_mean[i-1900] <- max(obs_mean[year==i])
}


plot(years, yearly_mean - mean(yearly_mean[1:20]))
fit <- smooth.spline(years, yearly_mean - mean(yearly_mean[1:20]))
lines(fit, col="purple", lwd="2")






