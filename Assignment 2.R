##### Import Data and Load Necessary Packages ----------------------------------
library(ncdf4)
library(lubridate)
library(RColorBrewer)
library(lattice)
library(fields)
library(terra)
library(data.table)

# Add an extra line for your personal file path and comment out as appropriate
#little change to check git push is working

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



############################################################
######## --------- Seasonal Periodicity --------- ##########
############################################################

# Monthly mean temp over Ireland coordinates
obs_ire_mean <- apply(temp_ire, 3, mean, na.rm = T)
obs_ire_mean
# Create data frame including relevant data
df_ire <- data.frame(
  date = date,
  mean_temps = obs_ire_mean
)
df_ire
df_ire$month <- as.numeric(format(df_ire$date, "%m"))
df_ire$year <- as.numeric(format(df_ire$date, "%Y"))

monthly_temps <- split(df_ire$mean_temps, df_ire$month)


#### As a start, look at seasonal data averaged over entire time period
yearly_temps_by_month <- sapply(monthly_temps, mean)

plot(1:12, yearly_temps_by_month, xlab = "Month" , ylab = "Mean Temp")
fit <- smooth.spline(1:12, yearly_temps_by_month)
lines(fit, col = "purple")
## Not really useful


#### Split data into 14 year sections and compute
split_data <- lapply(monthly_temps, function(x) {
  intervals <- ceiling(seq_along(x) / 14)
  tapply(x, intervals, mean)
} )

split_data_matrix = matrix(0, nrow = 12, ncol = 8)
for(i in 1:12){
  split_data_matrix[i,] <- split_data[[i]]
}

fits = list()
for(i in 1:8){
  fits[[i]] <- smooth.spline(1:12, split_data_matrix[,i])
}

plot(NULL,
     xlim = c(1,12),
     ylim = range(split_data_matrix),
     xaxt = "n",
     xlab = "Month",
     ylab = "Mean Temperature")
axis(1, at = 1:12, labels = month.abb)

# Manually put in values 1:8 corresponding to 14 year blocks
lines(fits[[8]], col = "brown")

# Amplitudes are increasing over blocks and overall larger mean temps each month


#### Fit Harmonic model to data 
# https://stats.stackexchange.com/questions/60994/fit-a-sinusoidal-term-to-data
# http://www-stat.wharton.upenn.edu/%7Estine/stat910/lectures/06_harmonic_regr.pdf
# A lot of sources recommend time series model 

sin_term <- sin((2*pi*1:12)/12)
cos_term <- cos((2*pi*1:12)/12)

harm_models = list()
for(i in 1:8){
  harm_models[[i]] <- lm(split_data_matrix[,i] ~ sin_term + cos_term)
}
fitted_vals <- lapply(harm_models, predict)

plot(NULL,
     xlim = c(1,12),
     ylim = range(split_data_matrix),
     xaxt = "n",
     xlab = "Month",
     ylab = "Mean Temperature")
axis(1, at = 1:12, labels = month.abb)
split_data_matrix
# Manually put in values 1:8 corresponding to 14 year blocks
plot(1:12, split_data_matrix[,8], col = "orange")
lines(fitted_vals[[1]], col = "orange")
lines(fitted_vals[[8]], col = "green")

# Calculate the amplitude and phase from harmonic fits. Intercept represents mean
# temperature over 14 year periods
amp = c()
phase = c()
for(i in 1:8){
  mod <- harm_models[[i]]
  amp[i] <- sqrt((coef(mod)[2])^2 + (coef(mod)[3])^2)
  phase[i] <- atan2(coef(mod)[3], coef(mod)[2])
}


#### Last two analyses were done using 14 year average. Instead fit harmonic
# models for each year and compare amplitudes and phase shifts 

year_month_temp_mat <- matrix(0, nrow = 112, ncol = 12)
for(i in 1:12){
  year_month_temp_mat[,i] <- monthly_temps[[i]]
}

beta0 <- c()
amp <- c()
phase <- c()

for(i in 1:nrow(year_month_temp_mat)){
  mod <- lm(year_month_temp_mat[i,] ~ sin_term + cos_term)
  coefs <- coef(mod)
  
  beta0[i] <- coefs[1]
  amp[i] <- sqrt( (coefs[2])^2 + (coefs[3])^2)
  phase[i] <- atan2(coefs[2], coefs[3])
  
}

# Plot the intercepts (Yearly mean temp) with a line representing the trend in 
# these fitted values
plot(beta0, type = "b" , col = "black")
abline(lm(beta0 ~ seq_along(beta0)), col = "red" )

# Note that models are of the form:
# Temp_m = beta0 + beta1*sin(2pi*month/12)+beta2*cos(2pi*month/12) + (error)
# But this can be rewritten in non-linear form as 
# Temp_m = beta0 + Amp*cos(2pi*m/12 - phase) + (error)
# So this is what is represented by amplitude and phase

# Plot amplitudes 
plot(amp, type = "b" , col = "black")
abline(lm(amp ~ seq_along(amp)), col = "red" )


# Plot phases
plot(phase %% 12 , type = "b" , col = "black")
abline(lm((phase %% 12) ~ seq_along((phase %% 12))), col = "red" )

# Use phases to find which month is temp highest and lowest to see if summer/winter
# seasons are starting earlier/later

# Highest temp occurs when 2*pi*m/12 - phase = 0 -> m = 12*phase/2*pi
# Because we're measuring months in rads, we can get monthly values not in [1,12]
# Modulo operator used to ensure correct values


highest_month <-  (phase * (12/(2*pi))) %% 12
plot(highest_month, type="b", main="Highest Temp Month")
abline(lm(highest_month ~ seq_along(highest_month)), col = "red" )

# Lowest temp occurs when 2*pi*m/12 - phase = pi -> m = (-12*phase + pi)/2*pi
lowest_month <- ((phase + pi) *12/(2*pi)) %% 12
plot(lowest_month, type="b", main="Lowest Temp Month")
abline(lm(lowest_month ~ seq_along(highest_month)), col = "red" )


# Extra Notes: 
# 1) Can be done using R-package:
# https://rdrr.io/cran/HarmonicRegression/man/harmonic.regression.html
# 2) Information/Intro for harmonic regression used above
# https://www.mdpi.com/1660-4601/17/4/1318
# This reference may also be useful for extensions



