library(ncdf4)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)

nc <- nc_open("C:/Users/86187/Desktop/Assignment 2/temp_data_SRS/temp_data_SRS.nc")

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