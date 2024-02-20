rm(list=(ls()))

# Load the required packages
library(fpp3)
library(seasonal)
library(urca)
library(strucchange)
library(forecast)
library(ncdf4)
library(ncdf4.helpers)
library(zoo)
library(sandwich)
library(gets)


####################
### DATA LOADING ###
####################

# Open the netCDF file
climate_filepath <- "/Users/MiraMetzger/Documents/CBS_Ms_Data_Science/FS2023/Predictive_Analytics/exam/artic/netcdf/ARCTIC_OMI_SI_extent_obs_19790115_P20220616_R19792021.nc"
climate_output <- nc_open(climate_filepath)
# read all values of that variable into an R object
# Monthly averaged Sea Ice Extent (km^2), Lake ice is not included
sea_ice_extent <- ncvar_get(climate_output, varid = "siextent")
ice_time <- nc.get.time.series(climate_output, time.dim.name = "time")

# convert to a tsibble object
sea_ice <- data_frame(time = ice_time, sea_ice_extent = as.vector(sea_ice_extent)) %>%
  mutate(time = as.Date(format(time, "%Y-%m-%d"))) %>%
  mutate(time = yearmonth(time)) %>% # time column is being converted from text to a monthly time object
  as_tsibble(index = time)

# see which months have no data
sea_ice %>% 
  filter(is.na(sea_ice_extent))
### there are missing values in the year of 1986 and 1987 --> only use data after the missing values

# only use data after Dec 1987 
sea_ice <- sea_ice %>%
  filter_index("1988 Jan"~"2021 Dec")


############################
#### VISUAL INVESTIGATION ##
############################

# plot of the data
sea_ice %>%
  autoplot(sea_ice_extent) +
  labs(y = expression("Sea Ice Extent ("~km^2~")"), x = "Time", title="Northern Hemisphere Sea Ice Extent (1988 to 2021)")

# autocorrelation
#--- data is not stationary (mean and variance are not constant)
#--- strong seasonality
#--- seems to have a decreasing trend
sea_ice %>%
  ACF(sea_ice_extent) %>%
  autoplot() +
  labs(title = "ACF")

# partial autocorrelation
sea_ice %>%
  PACF(sea_ice_extent) %>%
  autoplot() +
  labs(title = "PACF")

### not stationary (mean is not constant)
### seems to have a downward trend
### ACF: significant lags, strong seasonal pattern visible
### in all plots, seasonality is dominant that it's hard to see anything else.
### strong seasonal factors affect the time series

# seasonal plots
sea_ice %>%
  gg_subseries(sea_ice_extent) +
  labs(y = expression("Sea Ice Extent ("~km^2~")"), x = "Time")
sea_ice %>%
  gg_season(sea_ice_extent) +
  labs(y = expression("Sea Ice Extent ("~km^2~")"), x = "Time")
sea_ice %>%
  gg_lag(sea_ice_extent)


#################################
### TRANSFORM & Decomposition ###
#################################

### perform mathematical transformation
### need to transform back after forecasts, e.g to calculate accuracy

### LOG transformation
### reasons: growth rate, less skewed distribution, linear (better for modelling)

# log transformation --> too strong, variance not constant over time
sea_ice %>%
  autoplot(log(sea_ice_extent)) +
  labs(title = "Log Transformation", x = "time", y = "log transformed sea ice extent")

### log transformation (Box_Cox with lambda=0) is too strong (variance not constant)
### this is why we need to use the box cox transformation instead where we find the optimal lamda

# As log is not optimal (too strong) --> Find optimal value of lambda
lambda <- sea_ice %>%
  features(sea_ice_extent, features = guerrero) %>%
  pull(lambda_guerrero)

# Create new column with transformed CPI
sea_ice <- sea_ice %>%
  mutate(sea_ice_extent.boxcox = box_cox(sea_ice_extent, lambda))

# plot the transformed data
sea_ice %>%
  autoplot(sea_ice_extent.boxcox) + 
  labs(title = "Box Cox Transformation (λ=2)", x = "time", y = "box cox transformed sea ice extent")

# autocorrelation 
sea_ice %>%
  ACF(sea_ice_extent.boxcox) %>%
  autoplot()

# partial autocorrelation
sea_ice %>%
  PACF(sea_ice_extent.boxcox) %>%
  autoplot()

# Decomposition Additive (with boxcox as variation is now stable over time)
sea_ice %>%
  model(classical_decomposition(sea_ice_extent.boxcox, type = "additive")) %>%
  components() %>%
  autoplot() +
  labs(x = "Time", title = "Classical Additive Decomposition")

### random part looks like white noise


######################
### STATIONARITY #####
######################

### test for stationarity
### the data needs to be stationary before applying ARIMA and ets

# visual inspection, looks like one seasonal differencing is enough to make the data stationary
sea_ice %>% autoplot(sea_ice_extent.boxcox %>% difference(lag = 12)) + labs(y = "Sea Ice Extent (Seasonal Difference)", x = "Time")
sea_ice %>% autoplot(sea_ice_extent.boxcox %>% difference())
sea_ice %>% autoplot(sea_ice_extent.boxcox %>% difference(lag = 12) %>% difference())

# KPSS -------------------------------------------------------------------------

# data shows clear trend --> start with KPSS with trend (type “tau” in KPSS).
summary(ur.kpss(sea_ice$sea_ice_extent.boxcox, type = "tau")) # Trend + drift: tau
### t < cv --> not reject null, data might be stationary and fluctuates around the trend

# KPSS with drift (type "mu" in KPSS)
summary(ur.kpss(sea_ice$sea_ice_extent.boxcox, type = "mu")) 
### t > cv (at 5% significance level)
### reject null, data is not stationary with fluctuations around the mean

# ADF --------------------------------------------------------------------------

# clear trend -> ADF test for a random walk with trend and drift (type "trend" in ADF)
summary(ur.df(sea_ice$sea_ice_extent.boxcox, type = "trend", selectlags = "AIC", lags = 24))
### tau3 -> t > cv -> not reject null
### phi2 -> t < cv -> not reject null
### phi3 -> t < cv -> not reject null
### cannot reject the null hypothesis of non-stationarity, hence, our data might be non-stationary

# For this reason, we move on to test for a random walk with drift
summary(ur.df(sea_ice$sea_ice_extent.boxcox, type = "drift", selectlags = "AIC", lags = 24))
### tau2 -> t > cv -> not reject null
### phi1 -> t < cv -> not reject null
### cannot reject the null hypothesis of non-stationarity, hence, our data might be non-stationary

# ADF test for a random walk without drift and trend
summary(ur.df(sea_ice$sea_ice_extent.boxcox, type = "none", selectlags = "AIC", lags = 24)) 
### tau1: t > cv -> not reject null
### cannot reject the null hypothesis of non-stationarity, hence, our data might be non-stationary


# ADF seasonal diff ------------------------------------------------------------
# apply seasonal differencing to make the data stationary
# the previous analysis indicaties a period of 12 months

# take the difference and start again, ADF with a deterministic trend, drift + unit root
summary(ur.df(diff(sea_ice$sea_ice_extent.boxcox, lag=12), type = "trend", selectlags = "AIC", lags = 24))
### tau3 -> t < cv -> reject null 
### phi2 -> t > cv -> reject null
### phi3 -> t > cv -> reject null
### rejects the null hypothesis of non-stationarity 

# ADF test random walk with drift
summary(ur.df(diff(sea_ice$sea_ice_extent.boxcox, lag = 12), type = "drift", selectlags = "AIC", lags = 24)) # Drift + unit root 
#tau2 -> t < cv -> reject null
#phi1 -> t > cv -> reject null

# ADF for a random walk without drift and trend
summary(ur.df(diff(sea_ice$sea_ice_extent.boxcox, lag = 12), type = "none", selectlags = "AIC", lags = 24)) # Drift + unit root 
#tau1 -> t < cv -> reject null
### data is stationariy --> we have a unit root without trend and drift

# KPSS seasonal diff -----------------------------------------------------------

# KPSS with trend, HO: data is stationary
summary(ur.kpss(diff(sea_ice$sea_ice_extent.boxcox, lag = 12), type = "tau")) # Trend + drift
### t < cv -> accept null at all sig. level -> no unit root, time series is stationary

# KPPS with drift
summary(ur.kpss(diff(sea_ice$sea_ice_extent.boxcox, lag = 12), type = "mu")) # Drift
### t < cv -> accept null at all sig. level -> no unit root, time series is stationary


######################
### STRUC. BREAK #####
######################

# QLR --------------------------------------------------------------------------

# QLR test to test for endogenous breaks
sea_ice <- sea_ice %>%
  mutate(l0 = sea_ice_extent.boxcox,
         l12 = lag(sea_ice_extent.boxcox, 12))

sea_ice %>%
  features(l0, features = shift_level_max )

# Calculate F statistics
qlr_sea_ice <- Fstats(l0 ~ l12, data = sea_ice, from = 0.15)

test_sea_ice <- sctest(qlr_sea_ice, type = "supF")
test_sea_ice

# check where a breakpoint would be
breakpoints_sea_ice <- breakpoints(qlr_sea_ice, alpha = 0.05)

# plot F statistics and upper bound
plot(qlr_sea_ice, alpha = 0.05, main = "F Statistics")
lines(breakpoints_sea_ice)

# Conclusion: no structural break visible
# with a p-value of 0.4208 --> hypothesis of no structural break cannot be rejected
# There is likely to be no structural break in the data

# SIS --------------------------------------------------------------------------

#sea_ice.zoo <- as.zoo(as.ts(sea_ice$sea_ice_extent))
#sis <- isat(sea_ice.zoo, ar = c(1:4), mc = TRUE)
#plot(sis)


###########################
#### TRAIN/ TEST SPLIT ####
###########################

# train test split (20% of the data in the test set)
length_train <- (length(sea_ice$sea_ice_extent.boxcox)*0.8) # 298.4 in train set
sea_ice[floor(length_train),] # have until dec 2014 as the train set the rest in the test

train <- sea_ice %>% 
  slice(0:324)
# 7 years in the test set
test <- sea_ice %>%
  slice(325:408)


###########################
######### SNAIVE ##########
###########################

# fit a seasonal naive model
fit_naive <- train %>%
  model (
    SNAIVE = SNAIVE(box_cox(sea_ice_extent, lambda))
  )

# report for the snaive model
report(fit_naive)

# Residuals analysis
fit_naive %>%
  gg_tsresiduals(type="innovation")


###########################
######## ETS MODELS #######
###########################

# trend: linear component (additive)
# seasonality: constant amplitude and fixed pattern of seasonal variation (additive)
# Error: additive

# fit a manual and auto ETS model
fit_ets <- train %>%
  model(
    auto = ETS(box_cox(sea_ice_extent, lambda)),
    manual = ETS(box_cox(sea_ice_extent, lambda) ~ error("A") + trend("A") + season("A")),
  )

# get an overview of the model performance
glance(fit_ets) %>%
  arrange(AICc)

# report for the ETS models
report(fit_ets%>%select(auto))
report(fit_ets%>%select(manual))

# Residuals analysis
fit_ets %>%
  select(auto) %>%
  gg_tsresiduals(type="innovation")
fit_ets %>%
  select(manual) %>%
  gg_tsresiduals(type="innovation")

# Ljung Box test to check if the residuals have no autocorrelation
fit_ets %>%
  select(auto) %>%
  residuals()  %>%
  features(.resid, features = ljung_box, lag=24)
fit_ets %>%
  select(manual) %>%
  residuals()  %>%
  features(.resid, features = ljung_box, lag=24)

# forecast and plot the best model
fit_ets %>%
  select(auto) %>%
  forecast(h = "5 years") %>%
  autoplot(sea_ice, level = NULL)


###########################
######### ARIMA ###########
###########################

# look at the autocorrelation (ACF) and partial autocorrelation (PACF)
# to determine the order of the AR and MA part
# choose an ARIMA model for the series and write down the equation of the model in the form ARIMA(p,d,q)(P,D,Q)[m].

# ACF for differenced and seasonal differentiated data
sea_ice %>%
  gg_tsdisplay(difference(sea_ice_extent.boxcox, lag = 12), plot_type = "partial", lag_max = 36)
# ACF
sea_ice %>%
  ACF(difference(sea_ice_extent.boxcox, lag = 12), lag_max = 36) %>%
  autoplot()
# PACF
sea_ice %>%
  PACF(difference(sea_ice_extent.boxcox, lag = 12), lag_max = 36) %>%
  autoplot()

### (p,d,q)(P,D,Q)12

### non-seasonal part ###
### d(0) --> no differencing
### AR(2) --> PACF two significant spikes at the beginning of the PACF before it enters the non-significant boundary)
### MA(3) --> ACF 3 significant spikes at the begining of the ACF plot before it goes into the non-significant part
### (2,0,3)

### seasonal part ###
### D(1) --> one seasonal differencing
### MA(1) ACF: only significant at lag 12, no significant spike at lag 24
### lag 12 in the ACF is negative, indicating a MA(1) component.
### (0,1,1)

# fit a auto and manual ARIMA model with the identified components
fit_arima <- train %>%
  model(
    auto = ARIMA(box_cox(sea_ice_extent, lambda), stepwise = FALSE, approx = FALSE),
    manual1 = ARIMA(box_cox(sea_ice_extent, lambda) ~ 0 + pdq(0,0,3) + PDQ(0,1,1)),
    manual2 = ARIMA(box_cox(sea_ice_extent, lambda) ~ 0 + pdq(2,0,0) + PDQ(0,1,1)),
    manual3 = ARIMA(box_cox(sea_ice_extent, lambda) ~ 0 + pdq(2,0,3) + PDQ(0,1,1))
  )

# get an overview of the model performance
glance(fit_arima) %>%
  arrange(AICc)

# report for the ARIMA models
report(fit_arima%>%select(auto))
report(fit_arima%>%select(manual1))
report(fit_arima%>%select(manual2))
report(fit_arima%>%select(manual3))

# Inspect the residuals --------------------------------------------------------

# goal: obtain residuals that behave like white noise and follow a normal distribution
# residuals should be uncorrelated, have zero mean, constant variance, normal distribution

# plot the residuals
fit_arima %>%
  select(auto) %>%
  gg_tsresiduals(type = "innovation")
fit_arima %>%
  select(manual3) %>%
  gg_tsresiduals(type = "innovation")
### residuals look well distributed and fluctuate around 0
### one significant spike in the ACF at lag 17
### distribution close to a normal distribution with a few outliers

# ljung-box test on the residuals, test wheter first h autocorrelation are different than 0
fit_arima %>%
  select(auto) %>%
  residuals()  %>%
  features(.resid, features = ljung_box, lag = 24, dof = 3) # lag = 2m -> 2*12 = 24, dof = degrees of freedom p + q
### p-value of 0.310 --> no autocorrelation

fit_arima %>%
  select(manual3) %>%
  residuals()  %>%
  features(.resid, features = ljung_box, lag = 24, dof = 5) # lag = 2m -> 2*12 = 24, dof = degrees of freedom p + q
### p-value of 0.505 --> no autocorrelation


# shapiro-wilk test
shapiro.test(fit_arima %>%
               select(auto) %>%
               residuals() %>%
               select(.resid) %>%
               as.ts())
### reject H0 p-value of 0.01614
shapiro.test(fit_arima %>%
               select(manual3) %>%
               residuals() %>%
               select(.resid) %>%
               as.ts())
### reject H0 p-value of 0.03399

# boxpierceTest??
# reject null if p close to 0 --> residuals are not independent
#Box.test(residuals, lag=10, fitdf = 0)

# Forecast ---------------------------------------------------------------------

fit_arima %>%
  select(auto) %>%
  forecast(test) %>%
  autoplot(sea_ice) + 
  labs(title = "Forecasts auto model")

fit_arima %>%
  select(manual3) %>%
  forecast(test) %>%
  autoplot(sea_ice) + 
  labs(title = "Forecasts manual model")


#################################
########## FORECAST ALL #########
#################################

# fit an ARIMA model with the SNAIVE, ETS and ARIMA model
fit_all <- train %>%
  model (
    SNAIVE = SNAIVE(box_cox(sea_ice_extent, lambda)),
    ETS = ETS(box_cox(sea_ice_extent, lambda) ~ error("A") + trend("N") + season("A")),
    ARIMA = ARIMA(box_cox(sea_ice_extent, lambda) ~ 0 + pdq(2,0,3) + PDQ(0,1,1))
  )

# plot the forecast for the ARIMA model
fit_all %>%
  select(ARIMA) %>%
  forecast(h="7 years") %>%
  autoplot(test) +
  labs(y = expression("Sea Ice Extent ("~km^2~")"), x = "Time")

# plot the forecast for the ETS model
fit_all %>%
  select(ETS) %>%
  forecast(h="7 years") %>%
  autoplot(test) +
  labs(y = expression("Sea Ice Extent ("~km^2~")"), x = "Time")

# plot the forecast for the seasonal naive model
fit_all %>%
  select(SNAIVE) %>%
  forecast(h="7 years") %>%
  autoplot(test) +
  labs(y = expression("Sea Ice Extent ("~km^2~")"), x = "Time")

# plot the forecast for all models together
fit_all %>%
  forecast(test) %>%
  autoplot(test, level = NULL) +
  labs(y = expression("Sea Ice Extent ("~km^2~")"), x = "Time")

# assess the accuracy of all models
fabletools::accuracy(fit_all %>% select(ARIMA) %>% forecast(test), sea_ice)
fabletools::accuracy(fit_all %>% select(ETS) %>% forecast(test), sea_ice)
fabletools::accuracy(fit_all %>% select(SNAIVE) %>% forecast(test), sea_ice)
