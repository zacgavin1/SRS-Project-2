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
image.plot(tmp[,,1344], col = rev(brewer.pal(10,"RdBu")), main ="Global Temp 01/2012")



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

test_lon_ii <- which(lon > -160 & lon < 150)
test_lat_ii <- which(lat > -80 & lat <80 )
temp_test <- temp[test_lon_ii, test_lat_ii,]

n_lat <- length(test_lat_ii)
n_lon <- length(test_lon_ii)

ints <- array(rep(0, n_lat * n_lon *56), c(n_lon, n_lat, 56))
amps <- array(rep(NA, n_lat * n_lon *56), c(n_lon, n_lat, 56))

for (y in 1:56){
  for (i in 1:n_lon){
    for (j in 1:n_lat){
      if (!is.na(temp_test[i,j,1])){
        mod <- lm(temp_test[i,j,(24*(y-1)+1):((24*(y-1))+12)] ~ sin_term + 
             cos_term+sin2_term+cos2_term)
        ints[i,j, y] <- coef(mod)[1]
        amps[i,j,y] <- .5*(max(mod$fitted.values)-min(mod$fitted.values))
      } else ints[i,j,y ] <- amps[i,j,y] <- NA
    }
    print(i)
  }
  print(y)
}


# if we tried to do this 112 times this would take forever 
# (and would be 720*360*112=29030400 lms to fit)
# maybe we could average over small areas to reduce dimensionality?
# or fit a model every second year
`


# image.plot(mods[,,1], col = rev(brewer.pal(10,"RdBu")))


intmod <- matrix(rep(0, n_lat*n_lon), n_lon, n_lat)
ampmod <- matrix(rep(0, n_lat*n_lon), n_lon, n_lat)

for (i in 1:n_lon){
  for (j in 1:n_lat){
    if (!is.na(temp_test[i,j,1])){
      intmod[i,j] <- coef(lm(ints[i,j,]~ seq_along(ints[i,j,])))[2]
      ampmod[i,j] <- coef(lm(amps[i,j,]~ seq_along(ints[i,j,])))[2]
    } else intmod[i,j] <- ampmod[i,j] <- NA
  }
}

#### Plotting the heatmaps
data <- read.csv("ints.csv")
my_matrix <- as.matrix(data)
ampmod <- as.matrix(read.csv("amps.csv"))

lim <- max(abs(intmod*56), na.rm = TRUE)
image.plot(intmod*56, col = rev(brewer.pal(10,"RdBu")), zlim=c(-lim,lim))

lim <- max(abs(ampmod*56), na.rm = TRUE)
image.plot(ampmod*56, col = rev(brewer.pal(10,"RdBu")), zlim=c(-lim,lim))

      


################################################################################
############## Comparing Raw Data to Model Approaches ##########################
################################################################################
sin1 <- sin((2*pi*1:12)/12); cos1 <- cos((2*pi*1:12)/12)
sin2 <- sin((4*pi*1:12)/12); cos2 <- cos((4*pi*1:12)/12)

n_lon <- length(lon)
n_lat <- length(lat)
n_years <- length(years)
year_month_index <- split(1:length(year), month)
raw_amp_array  <- array(NA, dim=c(n_lon, n_lat, n_years))
raw_peak_array <- array(NA, dim=c(n_lon, n_lat, n_years))
model_amp_array  <- array(NA, dim=c(n_lon, n_lat, n_years))
model_peak_array <- array(NA, dim=c(n_lon, n_lat, n_years))


for(i in 1:n_lon){
  for(j in 1:n_lat){
    
    cell_ts <- temp[i,j,]
    if(all(is.na(cell_ts))) next
    
    # Build 112 x 12 matrix
    year_mat <- matrix(NA, nrow=n_years, ncol=12)
    for(m in 1:12){
      year_mat[,m] <- cell_ts[ year_month_index[[m]] ]
    }
    
    for(y in 1:n_years){
      yvals <- year_mat[y,]
      raw_amp_array[i,j,y]  <- (max(yvals) - min(yvals))/2
      raw_peak_array[i,j,y] <- which.max(yvals)
      
      mod <- lm(yvals ~ sin1 + cos1 + sin2 + cos2)
      fitted_vals <- fitted(mod)
      
      model_amp_array[i,j,y]  <- (max(fitted_vals)-min(fitted_vals))/2
      model_peak_array[i,j,y] <- which.max(fitted_vals)
    }
  }
}

peak_diff_array <- model_peak_array - raw_peak_array
mean_peak_diff <- apply(peak_diff_array, c(1,2), mean, na.rm=TRUE)

image.plot(lon, lat, mean_peak_diff,
           main="Model - Raw Peak Month Difference")
################################################################################
#### Branch YZ #################################################################
################################################################################


library(ncdf4)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)

lon  <- ncvar_get(nc, "lon")
lat  <- ncvar_get(nc, "lat")
time <- ncvar_get(nc, "time")
tmp  <- ncvar_get(nc, "tmp")

## Data preprocessing
# Dealing with Missing values: Since the missing values have been assigned extreme values, 
# it is planned to replace the missing values with NA.
fv <- ncatt_get(nc, "tmp", "_FillValue")$value
if (is.null(fv) || is.na(fv)) fv <- ncatt_get(nc, "tmp", "missing_value")$value
tmp[tmp >= (fv * 0.9)] <- NA

# Standardized time
origin <- as.Date("1900-01-01")
data <- origin + time
# Extract the date and month
years  <- year(data)
months <- month(data)

# Convert to standard longitude and latitude
if (max(lon, na.rm = TRUE) > 180) {
  lon_new <- ifelse(lon > 180, lon - 360, lon)
  ord <- order(lon_new)
  lon <- lon_new[ord]
  tmp <- tmp[ord, , ]
}

## Limit the area to Ireland
lon_min <- -10.3; lon_max <- -5.3
lat_min <-  50.3; lat_max <- 55.3

lon_idx <- which(lon >= lon_min & lon <= lon_max)
lat_idx <- which(lat >= lat_min & lat <= lat_max)

tmp_ie <- tmp[lon_idx, lat_idx, , drop = FALSE]

# Calculate the average temperature of the Irish region
ie_month_temperature <- sapply(seq_along(time), function(ti) {
  mean(tmp_ie[, , ti], na.rm = TRUE)
})

# Monthly time series data frame
ie_table <- tibble(
  date  = data,
  year  = years,
  month = months,
  ie_mean_month_temp = ie_month_temperature
)
ie_month_temperature
# head(ie_table)
# summary(ie_table$ie_mean_month_temp)

## Calculate the temperature difference between the coldest and hottest months each year (There is no obvious finds)
dif_year <- summarise(
  group_by(ie_table, year),
  amp = max(ie_mean_month_temp, na.rm = TRUE) - min(ie_mean_month_temp, na.rm = TRUE),
  .groups = "drop"
)

plot(dif_year$year, dif_year$amp, type = "l",
     xlab = "Year", 
     ylab = "Temperature Difference (max month - min month, °C)",
     main = "Ireland Temperature Difference over time"
)
points(dif_year$year, dif_year$amp, pch = 16, cex = 0.4)

abline(lm(amp ~ year, data = dif_year), lwd = 2)

## Whether the warmest or coldest month is changing (There is no obvious finds)
wc_month <- summarise(
  group_by(ie_table, year),
  warm_month = month[which.max(ie_mean_month_temp)],
  cold_month = month[which.min(ie_mean_month_temp)],
  .groups = "drop"
)

# table(wc_month$warm_month)
# table(wc_month$cold_month)

plot_warm <- ggplot(wc_month, aes(x = year, y = warm_month)) +
  geom_point(alpha = 0.5, position = position_jitter(height = 0.1, width = 0)) +
  scale_y_continuous(breaks = 1:12) +
  labs(x = "Year", y = "Warmest month")

plot_cold <- ggplot(wc_month, aes(x = year, y = cold_month)) +
  geom_point(alpha = 0.6, position = position_jitter(height = 0.1, width = 0)) +
  scale_y_continuous(breaks = 1:12) +
  labs(x = "Year", y = "Coldest Month")

plot_warm
plot_cold

## Determine the impact of climate warming on the seasons 
## (Whether the warming occurs uniformly throughout the year or is mainly concentrated in certain seasons)
## There is some finds
ie_season <- mutate(
  ie_table,
  season = case_when(
    month %in% c(12, 1, 2) ~ "Winter",
    month %in% c(3, 4, 5)  ~ "Spring",
    month %in% c(6, 7, 8)  ~ "Summer",
    TRUE                   ~ "Autumn"
  )
)

season_trend <- summarise(
  group_by(ie_season, season, year),
  temp = mean(ie_mean_month_temp, na.rm = TRUE),
  .groups = "drop"
)

# Perform a regression comparison of the slope for each season with respect to year
by(season_trend, season_trend$season, \(d) summary(lm(temp ~ year, data = d)))

ggplot(season_trend, aes(x = year, y = temp)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ season, ncol = 2, scales = "free_y") +
  labs(x = "Year", y = "Seasonal Mean Temperature (°C)")

## Calculate the annual average temperatures of other countries at different longitudes
## Define a general function. Given the bbox, output the annual average temperature.
annual <- function(lon, lat, tmp, data, years, lon_min, lon_max, lat_min, lat_max, name = "Region") {
  lon_idx <- which(lon >= lon_min & lon <= lon_max)
  lat_idx <- which(lat >= lat_min & lat <= lat_max)
  
  if (length(lon_idx) == 0 || length(lat_idx) == 0) {
    stop(sprintf("There are no grid points within the bbox. Check the lon or lat range.", name))
  }
  
  tmp_sub <- tmp[lon_idx, lat_idx, , drop = FALSE]
  
  # Monthly time series data frame
  ie_table <- tibble(
    date  = data,
    year  = years,
    month = months,
    ie_mean_month_temp = ie_month_temperature
  )
  
  # Regional monthly average temperature
  month_temperature <- sapply(seq_len(dim(tmp_sub)[3]), function(ti) {
    mean(tmp_sub[, , ti], na.rm = TRUE)
  })
  
  re_month <- tibble(
    date = data,
    year = years,
    temp = month_temperature
  )
  
  # Annual average temperature
  re_annual <- mutate(
    summarise(
      group_by(re_month, year),
      temp_annual = mean(temp, na.rm = TRUE),
      .groups = "drop"
    ),
    region = name
  )
  
  return(re_annual)
}

# Ireland

ie_annual <- annual(
  lon = lon, lat = lat, tmp = tmp,
  data = data, years = years,
  lon_min = lon_min, lon_max = lon_max,
  lat_min = lat_min, lat_max = lat_max,
  name = "Ireland"
)

# Only change the longitude range
# Netherlands（The latitude range remains at 50.3 – 55.3）
nl_lon_min <- 3.0;  nl_lon_max <- 7.5
nl_lat_min <- lat_min; nl_lat_max <- lat_max

nl_annual <- annual(
  lon = lon, lat = lat, tmp = tmp,
  data = data, years = years,
  lon_min = nl_lon_min, lon_max = nl_lon_max,
  lat_min = nl_lat_min, lat_max = nl_lat_max,
  name = "Netherlands (same-lat band)"
)

# Belarus
by_lon_min <- 23.0; by_lon_max <- 32.0
by_lat_min <- lat_min; by_lat_max <- lat_max

by_annual <- annual(
  lon = lon, lat = lat, tmp = tmp,
  data = data, years = years,
  lon_min = by_lon_min, lon_max = by_lon_max,
  lat_min = by_lat_min, lat_max = by_lat_max,
  name = "Belarus (same-lat band)"
)

# Draw a comparison chart
annuals <- bind_rows(ie_annual, nl_annual, by_annual)

# print(head(annuals))
# print(tail(annuals))

ggplot(annuals, aes(x = year, y = temp_annual, group = region, linetype = region)) +
  geom_line(linewidth = 0.8) +
  labs(x = "Year",
       y = "Annual mean temperature (°C)",
       title = "Annual mean temperature: Ireland vs same-latitude bands in Netherlands & Belarus",
       linetype = "Region"
  )



################################################################################
#### Branch scx ################################################################
################################################################################


# warm/cold seasons
obs_ire_mean <- apply(temp_ire, 3, mean, na.rm = TRUE)

year  <- as.integer(format(date, "%Y"))
month <- as.integer(format(date, "%m"))
years <- sort(unique(year))

warm_months <- 3:8
cold_months <- c(9:12, 1:2)

# 1 number per year
warm_mean <- cold_mean <- amp_mean <- rep(NA, length(years))

for (i in seq_along(years)) {
  y <- years[i]
  ii_warm_y <- which(year == y & month %in% warm_months)
  ii_cold_y <- which(year == y & month %in% cold_months)
  
  warm_mean[i] <- mean(obs_ire_mean[ii_warm_y], na.rm = TRUE)
  cold_mean[i] <- mean(obs_ire_mean[ii_cold_y], na.rm = TRUE)
  amp_mean[i]  <- warm_mean[i] - cold_mean[i]
}

# Warm vs Cold plot
op <- par(no.readonly = TRUE)
layout(matrix(c(1,2), nrow=2), heights=c(4,1))
par(mar=c(3,4,2,1))

ylim_all <- range(c(warm_mean, cold_mean), na.rm=TRUE)
plot(years, warm_mean, type="l", ylim=ylim_all,
     ylab="Mean Temp (°C)", xlab="",
     main="Ireland: Warm (Mar–Aug) vs Cold (Sep–Feb) seasonal means")
lines(years, cold_mean, lty=2)

par(mar=c(0,0,0,0))
plot.new()
legend("center",
       legend=c("Warm (Mar–Aug)", "Cold (Sep–Feb)"),
       lty=c(1,2), bty="n", horiz=TRUE, cex=0.9)

layout(1)
par(op)

# (warm - cold) and linear trend
par(mar=c(3,4,2,1))
plot(years, amp_mean, type="l",
     ylab="Warm - Cold (°C)", xlab="Year",
     main="Ireland: Seasonal contrast (mean-based)")
abline(lm(amp_mean ~ years), lty=2)



# Trend summaries
summary(lm(warm_mean ~ years))
##Slope: 0.006595 °C per year, p-value: 5.76e-06
##Perhaps could think that the warm season in Ireland is experiencing significant warming.

summary(lm(cold_mean ~ years))
##Slope: 0.005395 °C per year, p-value: 0.000797, warms up more slowly

summary(lm(amp_mean  ~ years))
##Slope: 0.001200 °C per year, p-value: 0.505 (not significant)
##Both seasons are getting warmer, but the temperature difference doesn't change significantly.






#############################################
## Remove average seasonal cycle -> anomalies
clim_by_month <- tapply(obs_ire_mean, month, mean, na.rm = TRUE)
anom <- obs_ire_mean - clim_by_month[as.character(month)]

# the time series, distribution, QQ and ACF 
par(mfrow=c(4,1), mar=c(3,4,2,1))
plot(date, anom, type="l", xlab="", ylab="Temp anomaly (°C)",
     main="Ireland: monthly anomalies (seasonality removed)")
abline(lm(anom ~ as.numeric(date)), lty=2)

hist(anom, breaks=40, main="Anomaly distribution", xlab="°C")
qqnorm(anom); qqline(anom)
## the distribution is not completely normal, with heavier tails.
## it indicates that there is a high possibility of extreme temperatures.

acf(anom, na.action=na.omit, main="ACF of anomalies")
## There is autocorrelation
## I think, it is reasonable because the warm/cold anomalies will persist for a certain period of time.
## But does this mean that the simple model for months is unreliable?

par(mfrow=c(1,1))

##########################
##issue:
# Should the influence of the months with extreme temperatures on the average monthly temperatures in cold and warm seasons be taken into account?
# Does autocorrelation shown by the ACF violate independence and invalidate p-values?




#############################
## 1. Are extreme hot and cold become more frequent over time?
df_ext <- data.frame(
  date  = date,
  year  = year,
  month = month,
  anom  = anom
)

df_ext$season_wc <- ifelse(df_ext$month %in% warm_months, "warm",
                           ifelse(df_ext$month %in% cold_months, "cold", NA))
df_ext <- df_ext[!is.na(df_ext$season_wc), ]
df_ext$season_wc <- factor(df_ext$season_wc, levels = c("cold","warm"))

# Process data

df_ext$t <- df_ext$year - mean(df_ext$year, na.rm = TRUE)

thr_hot_warm  <- quantile(df_ext$anom[df_ext$season_wc=="warm"], 0.95, na.rm = TRUE) # Determine the threshold using 0.95
thr_cold_cold <- quantile(df_ext$anom[df_ext$season_wc=="cold"], 0.05, na.rm = TRUE)

df_ext$Ehot_warm <- NA_integer_
df_ext$Ecold_cold <- NA_integer_

df_ext$Ehot_warm[df_ext$season_wc=="warm"] <- as.integer(df_ext$anom[df_ext$season_wc=="warm"] > thr_hot_warm)
df_ext$Ecold_cold[df_ext$season_wc=="cold"] <- as.integer(df_ext$anom[df_ext$season_wc=="cold"] < thr_cold_cold)

df_warm <- df_ext[df_ext$season_wc=="warm", ]
df_cold <- df_ext[df_ext$season_wc=="cold", ]

#  Logistic regression - season + trend + s * t

# warm
g_hot_warm <- glm(Ehot_warm ~ t, data = df_warm, family = binomial())
summary(g_hot_warm)
## time coefficient t = 0.0189, p = 0.0017
## the probability of "extreme heat anomalies" during the warm season has significantly increased over time.

#cold
g_cold_cold <- glm(Ecold_cold ~ t, data = df_cold, family = binomial())
summary(g_cold_cold)
## the time coefficient t is approximately 0.00021, and p = 0.97
## there are no significant long-term changes.

# check
# group - 10 years
# abs - the actual proportion of extreme hot in a group
# fit - the probability of model fitting
df_warm$decade <- floor(df_warm$year/10)*10
df_cold$decade <- floor(df_cold$year/10)*10

cal_warm <- aggregate(cbind(obs = Ehot_warm, fit = fitted(g_hot_warm)) ~ decade,
                      data = df_warm, FUN = mean)
cal_cold <- aggregate(cbind(obs = Ecold_cold, fit = fitted(g_cold_cold)) ~ decade,
                      data = df_cold, FUN = mean)

cal_warm
cal_cold

## obs clearly shows pulsation
## the sample size is extremely small in itself. I think this situation is quite reasonable
## so, not intuitive

# plot
## observed: number of ex-hot/ex-cold months in the corresponding six-month period of that year( n/6 )
## smoothed obs: line by smoothing the observation points
## fitted: the logistic regression model provides

yr_warm <- aggregate(cbind(obs = Ehot_warm, fit = fitted(g_hot_warm)) ~ year,
                     data = df_warm, FUN = mean)
yr_cold <- aggregate(cbind(obs = Ecold_cold, fit = fitted(g_cold_cold)) ~ year,
                     data = df_cold, FUN = mean)

op <- par(no.readonly = TRUE)
par(mfrow=c(2,1), mar=c(3,4,2,1))

plot(yr_warm$year, yr_warm$obs, type="p",
     xlab="Year", ylab="P(extreme hot | warm)",
     main="Warm season: yearly observed vs fitted")
lines(yr_warm$year, yr_warm$fit, lty=2)
lines(smooth.spline(yr_warm$year, yr_warm$obs), lty=1)
legend("topleft",
       legend=c("Observed (yearly)", "Fitted (model)", "Smoothed obs"),
       pch=c(1, NA, NA), lty=c(NA,2,1), bty="n", cex=0.9)

plot(yr_cold$year, yr_cold$obs, type="p",
     xlab="Year", ylab="P(extreme cold | cold)",
     main="Cold season: yearly observed vs fitted")
lines(yr_cold$year, yr_cold$fit, lty=2)
lines(smooth.spline(yr_cold$year, yr_cold$obs), lty=1)
legend("topleft",
       legend=c("Observed (yearly)", "Fitted (model)", "Smoothed obs"),
       pch=c(1, NA, NA), lty=c(NA,2,1), bty="n", cex=0.9)

par(op)

##The extreme hot has increased over time, while the extreme cold has not.

###############################
## 2. sensitivity test for the influence of extreme values
## 3 models compare: median, trimmed mean(delete the top and bottom 10%), winsorized mean(<5% -- =5%, >95% -- =95%)

# winsorize
winsorize <- function(x, p = 0.05) {
  q <- quantile(x, probs = c(p, 1 - p), na.rm = TRUE)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  x
}

# model
warm_med  <- cold_med  <- amp_med  <- rep(NA, length(years))
warm_trim <- cold_trim <- amp_trim <- rep(NA, length(years))
warm_win  <- cold_win  <- amp_win  <- rep(NA, length(years))

for (i in seq_along(years)) {
  y <- years[i]
  ii_warm_y <- which(year == y & month %in% warm_months)
  ii_cold_y <- which(year == y & month %in% cold_months)
  
  xw <- obs_ire_mean[ii_warm_y]
  xc <- obs_ire_mean[ii_cold_y]
  
  # median
  warm_med[i] <- median(xw, na.rm = TRUE)
  cold_med[i] <- median(xc, na.rm = TRUE)
  amp_med[i]  <- warm_med[i] - cold_med[i]
  
  # trimmed mean
  warm_trim[i] <- mean(xw, trim = 0.10, na.rm = TRUE)
  cold_trim[i] <- mean(xc, trim = 0.10, na.rm = TRUE)
  amp_trim[i]  <- warm_trim[i] - cold_trim[i]
  
  # winsorized mean
  warm_win[i] <- mean(winsorize(xw, p = 0.05), na.rm = TRUE)
  cold_win[i] <- mean(winsorize(xc, p = 0.05), na.rm = TRUE)
  amp_win[i]  <- warm_win[i] - cold_win[i]
}

# Trend comparisons
summ_line <- function(mod) c(slope = coef(mod)[2], p = summary(mod)$coefficients[2,4])

res_tab <- rbind(
  warm_mean  = summ_line(lm(warm_mean ~ years)),
  warm_med   = summ_line(lm(warm_med  ~ years)),
  warm_trim  = summ_line(lm(warm_trim ~ years)),
  warm_win   = summ_line(lm(warm_win  ~ years)),
  
  cold_mean  = summ_line(lm(cold_mean ~ years)),
  cold_med   = summ_line(lm(cold_med  ~ years)),
  cold_trim  = summ_line(lm(cold_trim ~ years)),
  cold_win   = summ_line(lm(cold_win  ~ years)),
  
  amp_mean   = summ_line(lm(amp_mean  ~ years)),
  amp_med    = summ_line(lm(amp_med   ~ years)),
  amp_trim   = summ_line(lm(amp_trim  ~ years)),
  amp_win    = summ_line(lm(amp_win   ~ years))
)

res_tab

# plot
op <- par(no.readonly = TRUE)
par(mar=c(4,4,2,1))

plot(years, amp_mean, type="l",
     ylab="Warm - Cold (°C)", xlab="Year",
     main="Ireland: amplitude trend sensitivity (mean vs robust)")
lines(years, amp_med,  lty=2)
lines(years, amp_trim, lty=3)
lines(years, amp_win,  lty=4)

legend("topleft",
       legend=c("mean","median","trim(10%)","winsor(5%)"),
       lty=c(1,2,3,4), bty="n", cex=0.9)

par(op)

## Hard to see clearly, make a smooth line graph
op <- par(no.readonly = TRUE)
par(mar=c(4,4,2,1))

plot(years, amp_mean, type="n",
     ylab="Warm - Cold (°C)", xlab="Year",
     main="Ireland: amplitude trend sensitivity (smoothed)")

lines(smooth.spline(years, amp_mean), lty=1)
lines(smooth.spline(years, amp_med),  lty=2)
lines(smooth.spline(years, amp_trim), lty=3)
lines(smooth.spline(years, amp_win),  lty=4)

legend("topleft",
       legend=c("mean","median","trim(10%)","winsor(5%)"),
       lty=c(1,2,3,4), bty="n", cex=0.9)

par(op)

## In addition to the median, the slopes and p-values of the mean, trim, and winsor are all close.
## And the smooth lines of the three are almost touching each other.
## The median does not provide a strong explanation for the target.
## Perhaps we can consider that extreme values have little impact on the conclusion.
