###
# Across the three regions, the monthly temperature series are dominated by seasonal variation.
# This motivates modelling trend and seasonality explicitly before defining extremes using residuals.
# We define an extreme-hot month as a residual exceeding the baseline Q95. 
# Under this definition, Siberia shows a clear and statistically significant increase in extreme-hot month frequency.
# The estimated decadal rate ratios are consistently above 1 across model specifications.
# The result is also broadly robust to the choice of baseline period.
# In contrast, Ireland and Singapore show no robust evidence of a long-run increase under the same definition.
# Their estimated decadal rate ratios are close to 1, and the confidence intervals typically span 1.
# Diagnostics further suggest strong residual autocorrelation in Singapore. 
# Adding an ARMA error structure therefore substantially improves model fit.
# For Siberia, the seasonal structure appears more complex than a simple sinusoid.
# This favours a flexible month-factor representation.
###
# Region: Singapore bbox

library(ncdf4)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)

nc_path <- "C:/Users/86187/Desktop/Assignment 2/temp_data_SRS/temp_data_SRS.nc"
vn <- "tmp"

# Choose latitude and longitude (Singapore bbox)
LON_MIN <- 103.6
LON_MAX <- 104.1
LAT_MIN <- 1.15
LAT_MAX <- 1.50

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
month_mean_tep <- function(nc_path, vn = "tmp", lon_min = 103.6, lon_max = 104.1, lat_min = 1.15,  lat_max = 1.50) {
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
# For Singapore, adding an ARMA error structure leads to a substantial improvement in fit.
# This indicates strong temporal correlation in the error term of the monthly temperature series.
# Without an ARMA component, this dependence would be left in the residuals, noticeably worsening model fit.
# Regarding seasonality, M1 and M3 have similar AIC values, whereas M2 performs worst,
# suggesting that a first-order harmonic is too crude to capture the seasonal pattern.
print(comparison)

# Residuals
df$res_m1 <- resid(m1)
df$res_m2 <- resid(m2)
df$res_m3 <- resid(m3)
df$res_m4 <- as.numeric(residuals(m4))

# Sensitivity analysis across alternative 30-year baselines shows that
# the substantive conclusion (no robust increasing trend) is not driven by this choice.
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

# It is found that no statistically detectable long-term increase in the frequency of extreme-hot months in Singapore from 1900–2012
# Estimated decadal rate ratios are close to 1, with 95% CIs generally spanning 1.
print(trend_table)

# Plots
# (1) monthly mean temperature
# The seasonal variation is very small (ranging roughly between 26–29), much smaller than that of Ireland and siberia.
ggplot(df, aes(x = date, y = temp)) +
  geom_line() +
  labs(title = "Singapore Monthly Mean Temperature",
       subtitle = paste0("Singapore lon[", LON_MIN, ",", LON_MAX, "], lat[", LAT_MIN, ",", LAT_MAX, "]"),
       x = "Date", y = "Temp (°C)") +
  theme_minimal()

# (2) Residuals with thresholds
# The occurrence of the extreme hot month does not increase significantly over time.
ggplot(df, aes(x = date, y = res_m4)) +
  geom_line() +
  geom_hline(yintercept = th_m4$q95, linetype = 2) +
  geom_hline(yintercept = th_m4$q05, linetype = 2) +
  labs(title = paste0("Singapore Residuals (M4) with Baseline Thresholds (", min(bly), "-", max(bly), ")"),
       subtitle = "Extreme Hot Month: residual > baseline Q95",
       x = "Date", y = "Residual (°C)") +
  theme_minimal()

# (3) Yearly extreme-hot months with poisson trend
# In the Singapore region, the annual frequency does not exhibit a clear long-term increasing trend.
# Across models, the Poisson trend lines are broadly flat (and in some cases slightly decreasing),
# with substantial year-to-year variability. 
# This suggests that there is no robust upward trend in the frequency of extreme hot months in Singapore.

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
  labs(title = "Singapore: Yearly Extreme Hot Months with Poisson Trend",
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