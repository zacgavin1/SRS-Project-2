#### Summary -------------------------------------------------------------------

# We need to emphasise what our research question is and put more structure on 
# the code we already have.

# I think the overarching questions we're looking at is:
# "How has/does climate change effect the shape and timing of seasonal cycles"
# "How has the seasonal structure of temperature changed over time"
# "Are seasonal peaks becoming more extreme over time"

# In terms of applications we can relate this to biological phenology where 
# certain biological patterns rely on stable seasonal patterns

# In particular we want to look at:
# 1) Is the seasonal differences between cold and warm months changing (changes
#    in amplitude) or are they increasing at the same rate 
# 2) Is the timing of peak yearly temperatures shifting (looking at phase shift)
#    (only having monthly data for this question may be an issue) 
# 3) Are extreme warm/cold months becoming more frequent 

# Ideally we could do this globally and plot to see if particular regions are 
# effected more drastically. From testing, inland regions appear to be affected
# more than coastal areas

# How do we do this ?

# Question (1) and (2):
# - From raw data we can calculate the amp = (max(yearly_temp) - min(yearly_temp))/2
#   and the peak months = argmax/argmin(yearly_temp)
# - Look at trends in these quantities over time using linear model
# - Consider harmonic regression models which account for overall trend and
#   seasonal effects. Account for highly correlated errors with ARMA component
#   This way we have a model of the form:
#   T_t = trend + seasonal + ARMA errors + random noise

# Question (3)
# - Define what is meant by extreme temperature. Maybe rolling window -> 95th 
#   percentile. Or a certain deviation from expected value from models ?
# - Logistic regression ?








# Structured code --------------------------------------------------------------

################################################################################
######################## Data Preprocessing and EDA  ###########################
################################################################################
library(ncdf4)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)


# Load in Variables
setwd("C:\\Users\\Luke Egan\\OneDrive\\Desktop\\Statistical Research Skills\\Assignment 2")
nc <- nc_open("temp_data_SRS.nc")
lon  <- ncvar_get(nc, "lon")
lat  <- ncvar_get(nc, "lat")
time <- ncvar_get(nc, "time")
tmp  <- ncvar_get(nc, "tmp")



# Standardized time and extract year and month variables
origin <- as.Date("1900-01-01")
data <- origin + time
years  <- year(data)
months <- month(data)

year <- unique(years)

# Global temperature plots of first and last recording
# par(mfrow = c(2,1), mar = c(3,3,2,5))
# image.plot(temp[,,1], col = rev(brewer.pal(10,"RdBu")), main = "Global Temp 01/1901")
# image.plot(temp[,,1344], col = rev(brewer.pal(10,"RdBu")), main ="Global Temp 01/2012")
# 
# Global temperature plot of 112 year difference
# image.plot(temp[,,1344] - temp[,,1], col = rev(brewer.pal(10,"RdBu")), 
#            main ="Temp differences 111 years apart")



# Initial look at global mean temperature increase over the time period 
obs_mean <- apply(tmp, 3, mean, na.rm=T) # this take a while
yearly_mean <- rep(NA, length(years))
for (i in years){
  yearly_mean[i-1900] <- mean(obs_mean[year==i])
}

plot(years, yearly_mean - mean(yearly_mean[1:20]), ylim=c(-.3,1.5))
fit <- smooth.spline(years, yearly_mean - mean(yearly_mean[1:20]))
lines(fit, col="purple", lwd="2") # wait, climate change is real??






# Define a testing region. These are the coordinates for central Siberia as 
# seasonal changes are more drastic here. Easier for testing

lon_ii <- which(lon > 90 & lon < 110)
lat_ii <- which(lat > 55 & lat < 65)
temp_test <- tmp[lon_ii, lat_ii,]

dim(temp_test)

n_lon <- length(lon_ii)
n_lat <- length(lat_ii)
n_time <- dim(temp_test)[3]



# Now take mean over 3x3 grids in this region. Lower spatial dimension

# n_lon_block <- floor(n_lon / 3)
# n_lat_block <- floor(n_lat / 3)
# 
# temp_test_3x3 <- array(NA, dim = c(n_lon_block,
#                                   n_lat_block,
#                                   n_time))
# 
# for (i in 1:n_lon_block) {
#   for (j in 1:n_lat_block) {
#     
#     lon_range <- ((i-1)*3 + 1):(i*3)
#     lat_range <- ((j-1)*3 + 1):(j*3)
#     
#     block_vals <- temp_test[lon_range,
#                            lat_range,
#     ]
#     
#     temp_test_3x3[i, j, ] <-
#       apply(block_vals, 3, mean, na.rm = TRUE)
#   }
# }
# 
# dim(temp_test_3x3)

################################################################################
######################## Seasonal Data Analysis   ##############################
################################################################################

## Amplitude from data 


## Phase from data



################################################################################
######################## Harmonic Regression Model  ############################
################################################################################

## Fit harmonic models

# Using all siberian grid points to test on
# Model over entire time period
# Define time period and frequency of periodic cycles

t <- 1:n_time
w <- 2*pi/12

sin1 <- sin(w * t); sin2 <- sin(2 * w * t); sin3 <- sin(3 * w * t)
cos1 <- cos(w * t); cos2 <- cos(2 * w * t); cos3 <- cos(3 * w * t)

n_lon_block <- dim(temp_test)[1]
n_lat_block <- dim(temp_test)[2]



# Compare Harmonic Regression Models using relative measures of performance 
# (AIC and BIC)

AIC_order1 <- AIC_order2 <- AIC_order3<- BIC_order1<-BIC_order2<-BIC_order3 <- matrix(0,n_lon_block, n_lat_block)
best_order_AIC <- best_order_BIC <- matrix(0,n_lon_block, n_lat_block)
for (i in 1:n_lon_block) {
  for (j in 1:n_lat_block) {
    
    y <- temp_test[i, j, ]
    
    # Remove missing values if needed
    # if (all(is.na(y))) next
    
    # First order
    mod1 <- lm(y ~ t + sin1 + cos1 + t:sin1 + t:cos1)
    
    # Second order
    mod2 <- lm(y ~ t + sin1 + cos1 + sin2 + cos2 + t:sin1 + t:cos1
               + t:sin2 + t:cos2)
    
    # Third order
    mod3 <- lm(y ~ t + sin1 + cos1 + sin2 + cos2 +  sin3 + cos3 +
                 t:sin1 + t:cos1 + t:sin2 + t:cos2 + t:sin3 + t:cos3)
    
    AIC_order1[i,j] <- AIC(mod1)
    AIC_order2[i,j] <- AIC(mod2)
    AIC_order3[i,j] <- AIC(mod3)
    
    BIC_order1[i,j] <- BIC(mod1)
    BIC_order2[i,j] <- BIC(mod2)
    BIC_order3[i,j] <- BIC(mod3)
    
    
    best_order_AIC[i,j] <- which.min(
      c(AIC_order1[i,j],
        AIC_order2[i,j],
        AIC_order3[i,j])
    )
    
    best_order_BIC[i,j] <- which.min(
      c(BIC_order1[i,j],
        BIC_order2[i,j],
        BIC_order3[i,j])
    )
  }
}
summary(mod3)
c(length(which(best_order_AIC == 1)), length(which(best_order_AIC == 2)),
  length(which(best_order_AIC == 3)))

c(length(which(best_order_BIC == 1)), length(which(best_order_BIC == 2)),
  length(which(best_order_BIC == 3)))


# AIC values favour 3rd order harmonics but reduction in AIC is not huge between
# 2nd order and 3rd order. This is reflected by the BIC which penalises 
# the additional parameters more heavily and favours the second order model 




# Look at the 2nd order model in more detail. In particular, we want to look at the
# time and harmonic interaction terms and see if they're significant. Significance
# may indicate potential seasonal changes


# NOTE: The significance of these terms is looking at seasonality evolution (i.e.
#       both amplitude and phase shift simultaneously)


# To test amplitude and phase in isolation I think the yearly models and looking
# at trends in amp and phase is more appropriate

models_list <- list()
interaction_pvals <- matrix(NA, n_lon_block, n_lat_block)

for (i in 1:n_lon_block) {
  for (j in 1:n_lat_block) {
    
    y <- temp_test[i, j, ]
    
    if (all(is.na(y))) next
    
    mod2_reduced <- lm(y ~ t + sin1 + cos1 + sin2 + cos2)
    
    mod2 <- lm(y ~ t + sin1 + cos1 + sin2 + cos2 + 
                 t:sin1 + t:cos1 + t:sin2 + t:cos2)
    
    models_list[[paste0("lon", i, "_lat", j)]] <- mod2
    
    Ftest <- anova(mod2_reduced, mod2)
    
    interaction_pvals[i, j] <- Ftest$`Pr(>F)`[2]
    
  }
}


# Benjamini-Hochberg Correction for Multiple Testing. Maybe need to look
# at assumptions for this because I think the tests may be spatially dependent and
# independence may be an assumption
p_adjusted <- matrix(p.adjust(as.vector(interaction_pvals), method="BH"),
                     nrow=nrow(interaction_pvals))

p_adjusted_sig <- p_adjusted < 0.05
p_adjusted_sig


lon_region <- lon[lon_ii]
lat_region <- lat[lat_ii]

lon_region
lat_region

map_p_vals = matrix(NA, nrow = length(lon), ncol = length(lat))
map_p_vals[lat_region, lon_region] = interaction_pvals

library(maps)
image.plot(lon,
           lat,
           map_p_vals,
           main = "Seasonal Interaction p-values (Siberia)",
           xlab = "Longitude",
           ylab = "Latitude")

image.plot(lon_region,
           lat_region,
           p_adjusted,
           xlab = "Longitude",
           ylab = "Latitude",
           main = "Seasonal Interaction p-values (Siberia)",
           col = topo.colors(50))


map("world",
    add = TRUE,
    col = "black",
    lwd = 1)



# Model for each year and look at model coefficients


# month = 1:12
# n_years <- n_time / 12
# sin1 <- sin(2 * pi * month / 12); cos1 <- cos(2 * pi * month / 12)
# sin2 <- sin(4 * pi * month / 12); cos2 <- cos(4 * pi * month / 12)
# 
# trend_array <- array(rep(0, n_lat_block * n_lon_block *112), c(n_lon_block, n_lat_block, 112))
# amps <- array(rep(NA, n_lat_block * n_lon_block *112), c(n_lon_block, n_lat_block, 112))
# peak_month_array <- array(rep(NA, n_lat_block * n_lon_block *112), c(n_lon_block, n_lat_block, 112))
# amp_trend  <- matrix(NA, n_lon_block, n_lat_block)
# mean_trend <- matrix(NA, n_lon_block, n_lat_block)
# 
# 
# for (year_idx in 1:20) {
#   for (i in 1:n_lon_block) {
#     for (j in 1:n_lat_block) {
#       
#       y_year <- temp_test_3x3[i, j, ((year_idx - 1) * 12 + 1):(year_idx * 12)]  
#       
#       if (all(is.na(y_year))) next
#       
#       # Fit second-order harmonic model
#       mod <- lm(y_year ~ sin1 + cos1 + sin2 + cos2)
#       
#       trend_array[i, j, year_idx] <- coef(mod)[1]
#       amps[i,j,year_idx] <- 0.5*(max(mod$fitted.values)-min(mod$fitted.values))
#       peak_month_array[i, j, year_idx] <- which.max(mod$fitted.values)
# 
#     }
#   }
# }
# 
# amp_trend  <- matrix(NA, n_lon_block, n_lat_block)
# mean_trend <- matrix(NA, n_lon_block, n_lat_block)
# 
# 
# for (i in 1:n_lon_block) {
#   for (j in 1:n_lat_block) {
#     
#     amp_series  <- amps[i, j, ]
#     mean_series <- trend_array[i, j, ]
#     
#     if (all(is.na(amp_series))) next
#     
#     # Linear trend in amplitude
#     mod_amp <- lm(amp_series ~ year)
#     amp_trend[i, j] <- coef(mod_amp)[2]
#     
#     # Linear trend in mean temperature
#     mod_mean <- lm(mean_series ~ year)
#     mean_trend[i, j] <- coef(mod_mean)[2]
#   }
# }
# 
# amp_trend_total  <- amp_trend  * n_years
# mean_trend_total <- mean_trend * n_years
# lim_amp <- max(abs(amp_trend_total), na.rm = TRUE)
# 
# image.plot(amp_trend_total,
#            col = rev(heat.colors(50)),
#            zlim = c(-lim_amp, lim_amp),
#            main = "Total Change in Seasonal Amplitude",
#            xlab = "Lon blocks",
#            ylab = "Lat blocks")
# 
# 
# lim_mean <- max(abs(mean_trend_total), na.rm = TRUE)
# 
# image.plot(mean_trend_total,
#            col = rev(heat.colors(50)),
#            zlim = c(-lim_mean, lim_mean),
#            main = "Total Change in Mean Temperature",
#            xlab = "Lon blocks",
#            ylab = "Lat blocks")



################################################################################
################### Accounting for AR Errors ###################################
################################################################################
resids <- list()

for (i in 1:n_lon_block) {
  for (j in 1:n_lat_block) {
    
    y <- temp_test[i, j, ]
    
    if (all(is.na(y))) next
    
    mod2 <- lm(y ~ t + sin1 + cos1 + sin2 + cos2 + 
                 t:sin1 + t:cos1 + t:sin2 + t:cos2)
    
    models_list[[paste0("lon", i, "_lat", j)]] <- mod2
    
    resids[[paste0("lon", i, "_lat", j)]] <- residuals(mod2)
    
  }
}

acf(resids[[110]])
pacf(resids[[100]])

# Based on residual plots i think:

# An MA component may be needed as damped cosine pattern in PACF
# I don't think an AR component is needed
# Is all the seasonality captured by the harmonic regression ?


xreg <- cbind(
  t,
  sin1, cos1,
  sin2, cos2,
  t*sin1, t*cos1,
  t*sin2, t*cos2
)

models_list2 <- list()
resids2 <- list()
for(i in 1:n_lon_block) {
  for (j in 1:n_lat_block) {
    
    y <- temp_test[i, j, ]
    
    if (all(is.na(y))) next
    
    
    fit <- arima(
      y,
      order = c(1,0,1),   # example ARMA(1,1)
      xreg = xreg,
      method = "ML"
    )
    
    models_list2[[paste0("lon", i, "_lat", j)]] <- fit
    resids2[[paste0("lon", i, "_lat", j)]] <- residuals(fit)
    
  }
}

dim(temp_test)
mod11 <- models_list[[1]]
mod12 <- models_list2[[1]]

acf(resids2[[1]])
pacf(resids2[[1]])

fit1 <- fitted(mod11)
fit2 <- temp_test[1,1,] - resids2[[1]]


plot(1:140, temp_test[1,1,701:840])
lines(fit1[701:840], col = "red")
lines(fit2[701:840], col = "blue")

sum(resids[[420]]^2) 
sum(resids2[[220]]^2)


# The RSS is consistently smaller for the model that accounts for correlation from
# month to month. Accounting for correlated residuals is effective






################################################################################
######################## Extreme Temperature/Anomalies  ########################
################################################################################



plot(years, res2_vec[[1]])
plot(year,res2_vec[[2]])
plot(res2_vec[[300]], type="l")

threshs <- matrix(vector("list",n_lon_block*n_lat_block), nrow=n_lon_block, ncol=n_lat_block)
# finding quantile for threshold, and 
for(i in 1:(n_lat_block*n_lon_block)){ 
  
    threshs[[i]] <- c(quantile(resids2[[i]], 0.025), quantile(res2_vec[[i]], 0.975))
    
}


# plot with blue line lower extreme, red upper extreme threshholds
# change loc for different locations
loc <- 1
plot(years, resids2[[loc]], type="l")
abline(h=threshs[[loc]][1], col="blue")
abline(h=threshs[[loc]][2], col="red")

hot_exts_recent <- cold_exts_recent <- rec_hotext_rate<- rec_coldext_rate<- rep(NA, n_lat_block*n_lon_block)
for (i in 1:(n_lat_block*n_lon_block)){
  hot_exts_recent[i] <- sum(resids2[[i]][800:1344]>threshs[[i]][2])
  cold_exts_recent[i]<- sum(resids2[[i]][800:1344]<threshs[[i]][1])
  rec_hotext_rate <- hot_exts_recent/(1344-800)
  rec_coldext_rate<- cold_exts_recent/(1344-800)
}

mean(rec_hotext_rate)
mean(rec_coldext_rate)



