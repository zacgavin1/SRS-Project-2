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

# setwd("C:\\Users\\Luke Egan\\OneDrive\\Desktop\\Statistical Research Skills\\Assignment 2")
setwd("C:/Users/zgavi/Documents/Edinburgh Term 2/SRS")
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
ireland_lon_ii <- which(lon > 26 & lon < 31)
ireland_lat_ii <- which(lat > 39 & lat <42 )

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
  yearly_mean[i-1900] <- mean(obs_mean[year==i])
}


plot(years, yearly_mean - mean(yearly_mean[1:20]), ylim=c(-.3,1.5))
fit <- smooth.spline(years, yearly_mean - mean(yearly_mean[1:20]))
lines(fit, col="purple", lwd="2")
# wait, climate change is real??


##### 
# Does this differ in the North Atlantic region?
#####

NA_lon_ii<-which(lon > -80 & lon < 30); NA_lat_ii<-which(lat > 0 & lat <80 )
east_lon_ii<-which(lon>50 & lon <180); east_lat_ii<-which(lat>-80&lat<80)

# setting up the yearly mean temperatures in each region for plotting
obs_meane <- apply(temp[east_lon_ii,east_lat_ii,], 3, mean, na.rm=T) # this takes a while
yearly_mean_east <- rep(NA, length(years))
for (i in years){
  yearly_mean_east[i-1900] <- mean(obs_meane[year==i], na.omit=T)
}

obs_meanNA <- apply(temp[NA_lon_ii,NA_lat_ii,], 3, mean, na.rm=T)
yearly_mean_NA <- rep(NA, length(years))
for (i in years){
  yearly_mean_NA[i-1900] <- mean(obs_meanNA[year==i], na.omit=T)
}

par(mfrow=c(1,2))
# plot where the baseline is average temp at start of 20th century 
plot(years, yearly_mean_east - mean(yearly_mean_east[1:20]), ylim=c(-.3,1.5), main="Eastern Hemisphere") 
fit <- smooth.spline(years, yearly_mean_east - mean(yearly_mean_east[1:20]))
lines(fit, col="purple", lwd="2") # this is eastern hem

plot(years, yearly_mean_NA - mean(yearly_mean_NA[1:20]), ylim=c(-.3,1.5), main="North Atlantic region") 
fit <- smooth.spline(years, yearly_mean_NA - mean(yearly_mean_NA[1:20]), penalty=1.5)
lines(fit, col="purple", lwd="2") # this is NA region

# There is a more pronounced increase in the North Atlantic region during 
# the 40s peak compared to the eastern hemisphere


############################################################ ###################
######## --------- Seasonal Periodicity --------- ########## ###################
############################################################ ###################

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
cols <- colorRampPalette(c("blue", "red"))(8)
lines(fits[[1]],col=cols[1] )
lines(fits[[2]],col=cols[2] )
lines(fits[[3]],col=cols[3] )
lines(fits[[4]],col=cols[4] )
lines(fits[[5]],col=cols[5] )
lines(fits[[6]],col=cols[6] )
lines(fits[[7]],col=cols[7] )
lines(fits[[8]],col=cols[8] )



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

plot(amp)
plot(phase)


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
modamp <- lm(amp~seq_along(amp)); summary(modamp)

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


# It seems that the idea here is 
# temp=(overall trend)+(seasonal effects)+(noise), and then we are 
# looking at how the seasonal affects change over time?


################################################
## - model checking & model generalisations - ##
################################################

par(mfrow=c(1,1))
monthly_mean <- sapply(monthly_temps, mean)
plot(monthly_mean, ylim=c(3,16))
mod_check <- lm(monthly_mean~sin_term+cos_term)
lines(mod_check$fitted.values)

plot(mod_check) # "U" shaped residual plot


# Try second (and third) order terms to remove errors

sin2_term <- sin((4*pi*1:12)/12); cos2_term <- cos((4*pi*1:12)/12)
sin3_term <- sin((6*pi*1:12)/12); cos3_term <- cos((6*pi*1:12)/12)
  
mod_2check <- lm(monthly_mean~sin_term+cos_term +sin2_term+cos2_term)
lines(mod_2check$fitted.values)
plot(mod_2check) # still not perfect, but gets rid of "U" errors

mod_3check <- lm(monthly_mean~sin_term+cos_term +sin2_term+cos2_term+
                 cos3_term+sin3_term)

# Model checking: compare AICs (recall: 2*p-2*log-lik )
AIC(mod_check) # 20.2
AIC(mod_2check) # -7.17
AIC(mod_3check) # -4.73
# suggesting that the model with the second order terms might be the
# best option for the split into 8/ yearly regressions
# Looking online, people comment that the temp peaks are too sharp to 
# be captured by the first order models, and so use second order


######################################
## - implement second order model - ##
######################################

sin_term <- sin((2*pi*1:12)/12); cos_term <- cos((2*pi*1:12)/12)
sin2_term <- sin((4*pi*1:12)/12); cos2_term <- cos((4*pi*1:12)/12)

## First for each of the 8 blocks 
harm_models2 = list()
for(i in 1:8){
  harm_models2[[i]] <- lm(split_data_matrix[,i] ~ sin_term + cos_term+
                            sin2_term+ cos2_term)
}
fitted_vals2 <- lapply(harm_models2, predict)


# plot to demonstrate the fit changing across time
plot(NULL,
     xlim = c(1,12),
     ylim = range(split_data_matrix),
     xaxt = "n",
     xlab = "Month",
     ylab = "Mean Temperature")
axis(1, at = 1:12, labels = month.abb)
# Manually put in values 1:8 corresponding to 14 year blocks
plot(1:12, split_data_matrix[,8], col = "orange")
lines(fitted_vals2[[1]], col = "orange")
lines(fitted_vals2[[8]], col = "green")
lines(fitted_vals[[8]], col="red") # comparison to first order


# amplitude and phase analysis changes a bit here
# NOT COMPLETE YET
amp = c()
phase = c()
for(i in 1:8){
  mod <- harm_models2[[i]]
  amp[i] <- sqrt((coef(mod)[2])^2 + (coef(mod)[3])^2)
  phase[i] <- atan2(coef(mod)[3], coef(mod)[2])
}

plot(amp)      # phase of 
plot(phase)
plot(sapply(harm_models2, function(x) coef(x)[1])) # intercept v time


## now for each year individually
year_month_temp_mat <- matrix(0, nrow = 112, ncol = 12)
for(i in 1:12){
  year_month_temp_mat[,i] <- monthly_temps[[i]]
}

beta0 <- c()
amp <- c()
phase <- c()

for(i in 1:nrow(year_month_temp_mat)){
  mod <- lm(year_month_temp_mat[i,] ~ sin_term + cos_term+
              sin2_term+cos2_term)
  coefs <- coef(mod)
  
  beta0[i] <- coefs[1]
  amp[i] <- (max(mod$fitted.values)-min(mod$fitted.values))*.5
  phase[i] <- which(mod$fitted.values==max(mod$fitted.values))
  
}

# fit models for each 
intmod <- lm(beta0~years); 
ampmod <- lm(amp~years)
phasemod <- lm(phase~years)

plot(years,beta0); lines(years, intmod$fitted.values); summary(intmod)
plot(years, amp); lines(years, ampmod$fitted.values); summary(ampmod)
plot(years, phase); lines(years, phasemod$fitted.values); summary(phasemod)




##################################################################
######## ---- Investigating seasonal ARIMA models ---- ###########
##################################################################
library(SWMPr)

# Steps done:
# 1) decompose time series into trend, seasonal and stationary parts
# 2) Then what?


ts_data <- ts(obs_ire_mean, start=c(1901,1), frequency=12) # making a class time series object
decomp <- decompose(ts_data, type = "additive")
plot(decomp)
plot(decomp$trend)      # get a similar overall shape to yearly mean over time
plot(decomp$seasonal[1:36]) # just a few cycles here

plot()

# looking at stationary part
acf(decomp$random, na.action=na.omit)
pacf(decomp$random, na.action=na.omit)


## I'm not convinced by any of this yet
# I don't quite know what Q I'm trying to answer here. 
# I'm not sure how this will let us get at a change in the seasonality



########################################################
########## ------ some other stuff ------- #############
########################################################

# Try to fit a hierarchical model of lm for amplitude -> harmonic regression for seasons
# for each point in space. Then plot colour map of the amplitude

sin_term <- sin((2*pi*1:12)/12); cos_term <- cos((2*pi*1:12)/12)
sin2_term <- sin((4*pi*1:12)/12); cos2_term <- cos((4*pi*1:12)/12)

test_lon_ii <- which(lon > 0 & lon < 41)
test_lat_ii <- which(lat > -40 & lat <80 )
temp_test <- temp[test_lon_ii, test_lat_ii,]

n_lat <- length(test_lat_ii)
n_lon <- length(test_lon_ii)

mods_start <- matrix(rep(0, n_lat * n_lon), nrow = n_lon, ncol = n_lat)
mods_end <- matrix(rep(0, n_lat * n_lon), nrow = n_lon, ncol = n_lat)

for (i in 1:n_lon){
  for (j in 1:n_lat){
    if (!is.na(temp_test[i,j,1])){
      mods_start[i,j] <- coef(lm(temp_test[i,j,1:12] ~ sin_term + 
                          cos_term+sin2_term+cos2_term))[1]
    } else mods[i,j] <- NA
  }
  print(i)
}

for (i in 1:n_lon){
  for (j in 1:n_lat){
    if (!is.na(temp_test[i,j,1])){
      mods_start[i,j] <- coef(lm(temp_test[i,j,1332:1344] ~ sin_term + 
                                   cos_term+sin2_term+cos2_term))[1]
    } else mods_end[i,j] <- NA
  }
  print(i)
}
# if we tried to do this 112 times this would take forever 
# (and would be 720*360*112=29030400 lms to fit)
# maybe we could average over small areas to reduce dimensionality?
# or fit a model every second year



image.plot(mods, col = rev(brewer.pal(10,"RdBu")))



