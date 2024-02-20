# Forecasting_Sea_Ice_Extent
This GitHub repository hosts code and documentation for a comprehensive study on forecasting sea ice extent in the Northern Hemisphere. The diminishing sea ice extent has critical implications for ecosystems and inhabitants in these regions. The study employs advanced time series analysis techniques, specifically comparing Autoregressive Integrated Moving Average (ARIMA) and Exponential Smoothing (ETS) models to forecast sea ice extent based on historical data.

## Exploratory Data Analysis and Preprocessing
The dataset, sourced from the E.U. Copernicus Marine Service, covers monthly mean sea ice extent measurements from 1979 to 2021. The data analysis involves mathematical transformations such as the Box-Cox transformation to stabilize variance. Classical additive decomposition reveals a strong seasonal influence and potential fluctuations attributed to varying temperatures.

## Methodology
Stationarity tests, structural break analysis, and exploratory data analysis guide the selection of suitable models. The study introduces a benchmark model, Seasonal Naïve, and compares it with ARIMA and ETS models. For ETS models, both manually selected and automated models are derived. Similarly, ARIMA models are chosen based on significant spikes in autocorrelation and partial autocorrelation functions.

## Results
The chosen ARIMA model outperforms both ETS and Seasonal Naïve models in terms of AIC, AICc, and BIC metrics. The GitHub repository includes code for model identification, training, and evaluation. The results suggest the importance of systematically assessing different forecasting approaches, emphasizing that the most complex model does not always yield the best performance.
