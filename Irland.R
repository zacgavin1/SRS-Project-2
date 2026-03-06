###
# Extreme hot months are defined as months in which, after controlling for long term trend and seasonality,
# the temperature residual exceeds the 95th percentile of the residual distribution in the baseline period.
# res = observed temp - model predicted temp
# Using residuals removes the components explained by trend and seasonality,
# thereby isolating anomalies and reducing the influence of long term warming.
# A month is classified as an extreme hot month if its residual is above the baseline Q95 threshold, and is coded as 1
# Otherwise it is coded as 0. This gives each month a binary 0/1 label.
###
# The baseline period is 1901 to 1930.
# A sensitivity analysis was conducted using multiple 30 year baseline windows：1901–1930、1931–1960、1961–1990、1981–2010
# The results show that the per decade rate ratios are all close to 1, their 95% confidence intervals all include 1, 
# and all p values are greater than 0.05. 
# This indicates that the conclusions are not sensitive to the choice of baseline, supporting the use of 1901–1930.
###
# A Poisson regression is used to test for a temporal trend because the outcome is a count.
# The model is log(lamda_y) = alpha + beta*year +log(n_months)
# where lamda_y is the expected number of extreme hot months in year y.
# The coefficient beta represents the change associated with a one year increase. 
# Because this effect is typically small,
# it is converted to a per decade change to match the long term scale of the research question.
###
# Region: Ireland bbox
###
# Pipeline:
#  1) Read NetCDF by time slices
#  2) Compute area-weighted monthly mean temp
#  3) Fit 4 models:
#     M1: year + factor(month)
#     M2: year + 1st harmonic
#     M3: year + 2nd harmonic
#     M4: M3 + ARMA errors (arima with xreg)
#  4) Residual-based extremes using baseline thresholds
#  5) Count extreme months per year + Poisson trend test
###

library(ncdf4)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)

nc_path <- "C:/Users/86187/Desktop/Assignment 2/temp_data_SRS/temp_data_SRS.nc"
vn <- "tmp"

# Choose latitude and longitude
LON_MIN <- -11
LON_MAX <- -5
LAT_MIN <-  51
LAT_MAX <-  56

# Convert NetCDF "days since 1900-01-01" time to Date
trsdate <- function(nc, time_var = "time") {
  t <- ncvar_get(nc, time_var)
  as.Date("1900-01-01") + t
}

# Make harmonic regressors for month seasonality
# Map the months to the periodic characteristics of sin/cos
# Since seasons are cyclical, December and January should be consecutive
mharmonics <- function(monthv, K = 2, period = 12) {
  out <- list()
  for (i in 1:K) {
    out[[paste0("sin", i)]] <- sin(2 * pi * i * monthv / period)
    out[[paste0("cos", i)]] <- cos(2 * pi * i * monthv / period)
  }
  as.data.frame(out)
}

# Compute area-weighted monthly mean temperature time series
month_mean_tep <- function(nc_path, vn = "tmp", lon_min = -11, lon_max = -5, lat_min = 51,  lat_max = 56) {
  nc <- nc_open(nc_path)
  on.exit(nc_close(nc), add = TRUE)
  
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  dates <- trsdate(nc, "time")
  
  nt <- length(dates)
  
  # transform [0,360] to [-180,180]
  lon_1 <- lon
  lon_1[lon_1 > 180] <- lon_1[lon_1 > 180] - 360
  
  # Select region indices
  lon_idx <- which(lon_1 >= lon_min & lon_1 <= lon_max)
  lat_idx <- which(lat >= lat_min & lat <= lat_max)
  
  # Area weights ~ cos(lat)
  wlat <- cos(lat[lat_idx] * pi / 180)
  W <- matrix(rep(wlat, each = length(lon_idx)), nrow = length(lon_idx), ncol = length(lat_idx))
  
  ir_mean <- rep(NA_real_, nt)
  
  lon_start <- min(lon_idx)
  lat_start <- min(lat_idx)
  
  for (ti in 1:nt) {
    slice <- ncvar_get(
      nc, vn,
      start = c(lon_start, lat_start, ti),
      count = c(length(lon_idx), length(lat_idx), 1)
    )
    
    # Replace fill/missing with NA
    slice[slice > 1e30] <- NA
    
    u <- !is.na(slice)
    if (any(u)) {
      ir_mean[ti] <- sum(slice[u] * W[u]) / sum(W[u])
    }
  }
  
  tibble(
    date = dates,
    year = year(dates),
    month = month(dates),
    temp = ir_mean
  ) %>% filter(!is.na(temp))
}

# Compute baseline thresholds from residuals
thresholds <- function(residualv, years, bly) {
  base <- residualv[years %in% bly]
  base <- base[!is.na(base)]
  list(
    q95 = as.numeric(quantile(base, 0.95, na.rm = TRUE)),
    q05 = as.numeric(quantile(base, 0.05, na.rm = TRUE))
  )
}

# Make extreme indicator from residual and thresholds
ext_indicator <- function(residualv, th) {
  as.integer(residualv > th$q95)
}

# Fit Poisson trend for yearly extreme counts
fit_trend <- function(countv, years, exposure) {
  glm(countv ~ years, family = poisson(link = "log"), offset = log(exposure))
}

# Extract "per 10-year multiplier" from Poisson model
ext_trend <- function(glm_fit) {
  b <- coef(glm_fit)["years"]
  se <- sqrt(vcov(glm_fit)["years", "years"])
  tibble(
    beta_year = b,
    mult_per_10yr = exp(b * 10),
    ci_low = exp((b - 1.96 * se) * 10),
    ci_high = exp((b + 1.96 * se) * 10),
    p_value = summary(glm_fit)$coefficients["years", "Pr(>|z|)"]
  )
}

# Build monthly mean temperature series
df <- month_mean_tep(nc_path, vn, lon_min = LON_MIN, lon_max = LON_MAX, lat_min = LAT_MIN, lat_max = LAT_MAX)

# Center year for numerical stability
df$year_c  <- df$year - mean(df$year)
df$t_index <- seq_len(nrow(df))

# Add harmonic seasonal regressors to represent monthly cyclic seasonality
H1 <- mharmonics(df$month, K = 1)
H2 <- mharmonics(df$month, K = 2)

# Combine original data with harmonic regressors
df_m2 <- bind_cols(df, H1)
df_m3 <- bind_cols(df, H2)

## fit model M1~4
# M1: year + factor(month)
m1 <- lm(temp ~ year_c + factor(month), data = df)

# M2: year + 1st harmonic
m2 <- lm(temp ~ year_c + sin1 + cos1, data = df_m2)

# M3: year + 2nd harmonic
m3 <- lm(temp ~ year_c + sin1 + cos1 + sin2 + cos2, data = df_m3)

# M4: M3 + ARMA errors
xreg_m4 <- as.matrix(df_m3 %>% select(year_c, sin1, cos1, sin2, cos2))

# choose ARMA(p,q) for M4 by AIC (p,q <= 2)
# restrict the candidate orders to low values to avoid overfitting
# For each (p,q), fitting the model by maximum likelihood, 
# and then select the (p,q) with the smallest AIC as the error-term structure for Model M4
grid <- expand.grid(p = 0:2, q = 0:2)
grid$aic <- NA_real_

for (i in seq_len(nrow(grid))) {
  p <- grid$p[i]; q <- grid$q[i]
  fit_try <- try(
    arima(df_m3$temp, order = c(p, 0, q),
          xreg = xreg_m4, include.mean = TRUE, method = "ML"),
    silent = TRUE
  )
  if (!inherits(fit_try, "try-error")) {
    grid$aic[i] <- fit_try$aic
  }
}

grid <- grid[order(grid$aic), ]
print(grid)

best <- grid[1, ]
arma_order <- c(best$p, 0, best$q)
cat("Chosen ARMA order:", arma_order, "\n")

# fit M4 with chosen order
m4 <- arima(df_m3$temp, order = arma_order, xreg = xreg_m4, include.mean = TRUE, method = "ML")

# Model comparison (AIC/BIC)
comparison <- tibble(
  model = c("M1: year + factor(month)",
            "M2: year + 1st harmonic",
            "M3: year + 2nd harmonic",
            paste0("M4: M3 + ARMA(", arma_order[1], ",", arma_order[2], ",", arma_order[3], ") errors")),
  AIC = c(AIC(m1), AIC(m2), AIC(m3), m4$aic),
  BIC = c(BIC(m1), BIC(m2), BIC(m3), AIC(m4, k = log(nrow(df_m3))))
)
# M4 performs best because it captures the temporal autocorrelation in the residuals, 
# whereas M2 performs worst because a first-order harmonic is too coarse to represent the detailed seasonal pattern of temperature
# This indicates that, even after controlling for the long-term trend and seasonality, 
# the temperature series still exhibits substantial residual autocorrelation
print(comparison)

# Residuals
df$res_m1 <- resid(m1)
df$res_m2 <- resid(m2)
df$res_m3 <- resid(m3)
df$res_m4 <- as.numeric(residuals(m4))

# Sensitivity analysis for different baseline periods
baseline_sets <- list(
  "1901-1930" = 1901:1930,
  "1931-1960" = 1931:1960,
  "1961-1990" = 1961:1990,
  "1981-2010" = 1981:2010
)

# Q95 for hot extremes
threshold_q95 <- function(residualv, years, bly) {
  base <- residualv[years %in% bly]
  base <- base[!is.na(base)]
  list(q95 = as.numeric(quantile(base, 0.95, na.rm = TRUE)))
}

trend_one <- function(residual_vec, years_vec, baseline_years) {
  th <- threshold_q95(residual_vec, years_vec, baseline_years)
  ext <- as.integer(residual_vec > th$q95)
  
  yearly <- tibble(years = years_vec, ext = ext) %>%
    group_by(years) %>%
    summarise(n_months = n(), extreme_months = sum(ext, na.rm = TRUE), .groups = "drop")
  
  fit <- glm(extreme_months ~ years, family = poisson(link = "log"), offset = log(n_months), data = yearly)
  
  ext_trend(fit) %>%
    mutate(q95 = th$q95) %>%
    select(q95, beta_year, mult_per_10yr, ci_low, ci_high, p_value)
}

sens_table <- bind_rows(
  lapply(names(baseline_sets), function(bn) {
    bly_i <- baseline_sets[[bn]]
    
    bind_rows(
      M1 = trend_one(df$res_m1, df$year, bly_i),
      M2 = trend_one(df$res_m2, df$year, bly_i),
      M3 = trend_one(df$res_m3, df$year, bly_i),
      M4 = trend_one(df$res_m4, df$year, bly_i),
      .id = "model"
    ) %>% mutate(baseline = bn, .before = 1)
  })
)

print(sens_table)
# Choose baseline
bly <- baseline_sets[["1901-1930"]]

# baseline thresholds
th_m1 <- thresholds(df$res_m1, df$year, bly)
th_m2 <- thresholds(df$res_m2, df$year, bly)
th_m3 <- thresholds(df$res_m3, df$year, bly)
th_m4 <- thresholds(df$res_m4, df$year, bly)

# extreme indicators
df$ext_m1 <- ext_indicator(df$res_m1, th_m1)
df$ext_m2 <- ext_indicator(df$res_m2, th_m2)
df$ext_m3 <- ext_indicator(df$res_m3, th_m3)
df$ext_m4 <- ext_indicator(df$res_m4, th_m4)

# Yearly extreme month counts and Poisson trend tests
# aggregate monthly 0/1 extreme indicators into yearly counts
yearly_ext <- df %>%
  group_by(year) %>%
  summarise(
    n_months = n(),
    count_m1 = sum(ext_m1, na.rm = TRUE),
    count_m2 = sum(ext_m2, na.rm = TRUE),
    count_m3 = sum(ext_m3, na.rm = TRUE),
    count_m4 = sum(ext_m4, na.rm = TRUE),
    .groups = "drop"
  )

# fit Poisson regressions for trend over years with exposure offset
p_m1 <- fit_trend(yearly_ext$count_m1, yearly_ext$year, yearly_ext$n_months)
p_m2 <- fit_trend(yearly_ext$count_m2, yearly_ext$year, yearly_ext$n_months)
p_m3 <- fit_trend(yearly_ext$count_m3, yearly_ext$year, yearly_ext$n_months)
p_m4 <- fit_trend(yearly_ext$count_m4, yearly_ext$year, yearly_ext$n_months)

# Extract interpretable trend metrics
trend_table <- bind_rows(
  M1 = ext_trend(p_m1),
  M2 = ext_trend(p_m2),
  M3 = ext_trend(p_m3),
  M4 = ext_trend(p_m4),
  .id = "model"
)

# The estimated 10-year rate ratios are all close to 1 with 95% CIs spanning 1 and p-values > 0.05 across M1–M4,
# indicating no statistically detectable long-term trend in the frequency of extreme-hot months for Ireland under this specification
print(trend_table)

# Plots
# (1) monthly mean temperature
# Figure (1) shows a strong and persistent annual seasonal cycle in Ireland’s monthly mean temperature (summer peaks and winter troughs). 
# The seasonal amplitude dominates the variation in the series. 
# Consequently, any underlying long-term trend is not easily visible. 
# This motivates modelling trend and seasonality explicitly, before defining extreme months using residuals.
ggplot(df, aes(x = date, y = temp)) +
  geom_line() +
  labs(title = "Irland Monthly Mean Temperature",
       subtitle = paste0("Irland lon[", LON_MIN, ",", LON_MAX, "], lat[", LAT_MIN, ",", LAT_MAX, "]"),
       x = "Date", y = "Temp (°C)") +
  theme_minimal()

# (2) Residuals with thresholds
# After accounting for long term trends and seasonal patterns, 
# extreme hot months defined as residuals above the baseline 95th percentile occur throughout the entire study period,
# but the plot shows no obvious clustering over time or a sustained upward trend.
ggplot(df, aes(x = date, y = res_m4)) +
  geom_line() +
  geom_hline(yintercept = th_m4$q95, linetype = 2) +
  geom_hline(yintercept = th_m4$q05, linetype = 2) +
  labs(title = paste0("Ireland Residuals (M4) with Baseline Thresholds (", min(bly), "-", max(bly), ")"),
       subtitle = "Extreme Hot Month: residual > baseline Q95",
       x = "Date", y = "Residual (°C)") +
  theme_minimal()

# (3) Yearly extreme-hot months with poisson trend
# The annual count of extreme hot months fluctuates between 0 and 4, with most years recording 0 to 1.
# Under a Poisson regression, the expected number of extreme months remains broadly stable over time or shows a slight decline,
# and the trend lines are consistent across all four models.
# This suggests that, under this definition of extremes,
# there is no detectable long term increase in the frequency of extreme hot months in Ireland from 1900 to 2012.

# create prediction data for model
make_pred <- function(yearly_ext, fit_pois, model_name, count_col) {
  # Predict on link scale to get standard errors (log scale)
  pr_link <- predict(fit_pois, type = "link", se.fit = TRUE)
  eta <- as.numeric(pr_link$fit)
  se<- as.numeric(pr_link$se.fit)

  out <- data.frame(
    year = yearly_ext$year,
    model = model_name,
    extreme_months = yearly_ext[[count_col]],
    fit = as.numeric(predict(fit_pois, type = "response")),
    lo  = exp(eta - 1.96 * se),
    hi  = exp(eta + 1.96 * se)
  )
  out
}

# Combine predictions for all models into one long table
pred_all <- rbind(
  make_pred(yearly_ext, p_m1, "M1", "count_m1"),
  make_pred(yearly_ext, p_m2, "M2", "count_m2"),
  make_pred(yearly_ext, p_m3, "M3", "count_m3"),
  make_pred(yearly_ext, p_m4, "M4", "count_m4")
)

ext_hot <- ggplot(pred_all, aes(x = year)) +
  # highlight baseline period used to define Q95 threshold
  annotate("rect", xmin = min(bly), xmax = max(bly), ymin = -Inf, ymax = Inf, alpha = 0.12) +
  # 95% confidence interval band for fitted mean
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18) +
  # fitted Poisson mean curve
  geom_line(aes(y = fit), linewidth = 1) +
  # observed yearly counts
  geom_line(aes(y = extreme_months), alpha = 0.85) +
  geom_point(aes(y = extreme_months), size = 1.0) +

  facet_wrap(~ model, ncol = 2, scales = "free_y") +
  labs(title = "Ireland: Yearly Extreme Hot Months with Poisson Trend",
    x = "Year",
    y = "Extreme Hot Months Per Year",
  ) +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 12),
    axis.title   = element_text(size = 10),
    axis.text    = element_text(size = 9),
    strip.text   = element_text(size = 10),
    plot.caption = element_text(size = 8)
  )

print(ext_hot)