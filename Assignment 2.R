##### Import Data and Load Necessary Packages ----------------------------------
library(ncdf4)
library(lubridate)
library(RColorBrewer)
library(lattice)
library(fields)

setwd("C:\\Users\\Luke Egan\\OneDrive\\Desktop\\Statistical Research Skills\\Assignment 2")
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
c(time[1]/356, time[1344]/365)
# ~112 years of data, however time step are not uniform (as month length varies)


# Change time format to DD/MM/YYYY
t_ustr <- strsplit(t_units$value, " ")
t_dstr <- strsplit(unlist(t_ustr)[3], "-")
date <- ymd(t_dstr) +ddays(time)
c(date[1], date[1344])
# Data ranges from 16/01/1901 to 16/12/2012


# Visualise Day earliest and latest days 
temp_day1 <- temp[,,1]
temp_day1344 <- temp[,,1344]


# Compare heatmaps same day 111 years apart
par(mfrow = c(2,1))
image.plot(temp_day1, col = rev(brewer.pal(10,"RdBu")), main = "Global Temp 01/1901")
image.plot(temp_day1344, col = rev(brewer.pal(10,"RdBu")), main ="Global Temp 01/2012")


# Difference 111 years apart 
diff_temp = temp_day1344 - temp_day1
image.plot(diff_temp, col = rev(brewer.pal(10,"RdBu")), main ="Temp differences 111 years apart")


# Split data by seasons
ii_winter <- grep("-12-|-01-|-02-", date)
ii_spring <- grep("-03-|-04-|-05-", date)
ii_summer <- grep("-06-|-07-|-08-", date)
ii_autumn <- grep("-09-|-10-|-11-", date)


# winter_temp <- temp[,,ii_winter]
# spring_temp <- temp[,,ii_spring]
# summer_temp <- temp[,,ii_summer]
# autumn_temp <- temp[,,ii_autumn]


# Consider one region for eda - Ireland 
# Longitide (-10, 40), Latitude (35, 72)
ireland_lon_ii <- which(lon > -10.5 & lon < -5.5)
ireland_lat_ii <- which(lat > 51.5 & lat < 55.5)

boxplot(as.vector(temp[ireland_lon_ii, ireland_lat_ii, ii_winter[1:50]]),
        as.vector(temp[ireland_lon_ii, ireland_lat_ii, ii_winter[51:100]]),
        as.vector(temp[ireland_lon_ii, ireland_lat_ii, ii_winter[101:150]]),
        as.vector(temp[ireland_lon_ii, ireland_lat_ii, ii_winter[151:250]]),
        as.vector(temp[ireland_lon_ii, ireland_lat_ii, ii_winter[251:300]]),
        as.vector(temp[ireland_lon_ii, ireland_lat_ii, ii_winter[301:350]]))


# Consider sampling places in units to reduce computational load
# Split into yearly chuncks and check boxplots along
# Take yearly averages ?

