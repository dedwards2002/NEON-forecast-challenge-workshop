## install.packages('remotes')
## install.packages('tidyverse') # collection of R packages for data manipulation, analysis, and visualisation
## install.packages('lubridate') # working with dates and times
## remotes::install_github('eco4cast/neon4cast') # package from NEON4cast challenge organisers to assist with forecast building and submission

# ------ Load packages -----
library(tidyverse)
library(lubridate)
#--------------------------#

# Change this for your model ID
# Include the word "example" in my_model_id for a test submission
# Don't include the word "example" in my_model_id for a forecast that you have registered (see neon4cast.org for the registration form)
my_model_id <- 'example_ID'

# --Model description--- #

# Add a brief description of your modeling approach

# -- Uncertainty representation -- #

# Describe what sources of uncertainty are included in your forecast and how you estimate each source.

#------- Read data --------
# read in the targets data
targets <- read_csv("https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D/aquatics-targets.csv.gz")

# read in the sites data
aquatic_sites <- read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-ci/refs/heads/main/neon4cast_field_site_metadata.csv") |>
  dplyr::filter(aquatics == 1)

focal_sites <- aquatic_sites |> 
  filter(field_site_subtype == 'Lake') |> 
  pull(field_site_id)

# Filter the targets
targets <- targets %>%
  filter(site_id %in% focal_sites,
         variable == 'temperature')
#--------------------------#



# ------ Weather data ------
met_variables <- c("air_temperature")

# Past stacked weather -----
weather_past_s3 <- neon4cast::noaa_stage3()

weather_past <- weather_past_s3  |> 
  dplyr::filter(site_id %in% focal_sites,
                datetime >= ymd('2017-01-01'),
                variable %in% met_variables) |> 
  dplyr::collect()

# aggregate the past to mean values
weather_past_daily <- weather_past |> 
  mutate(datetime = as_date(datetime)) |> 
  group_by(datetime, site_id, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  # convert air temperature to Celsius if it is included in the weather data
  mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
  pivot_wider(names_from = variable, values_from = prediction)

# Future weather forecast --------
# New forecast only available at 5am UTC the next day
forecast_date <- Sys.Date() 
noaa_date <- forecast_date - days(1)

weather_future_s3 <- neon4cast::noaa_stage2(start_date = as.character(noaa_date))

weather_future <- weather_future_s3 |> 
  dplyr::filter(datetime >= forecast_date,
                site_id %in% focal_sites,
                variable %in% met_variables) |> 
  collect()

weather_future_daily <- weather_future |> 
  mutate(datetime = as_date(datetime)) |> 
  # mean daily forecasts at each site per ensemble
  group_by(datetime, site_id, parameter, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  # convert air temperature to Celsius if it is included in the weather data
  mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
  pivot_wider(names_from = variable, values_from = prediction) |> 
  select(any_of(c('datetime', 'site_id', met_variables, 'parameter')))

#--------------------------#



# ----- Fit model & generate forecast----

# 1. Prepare historical data with the Lag and Derived Covariate
targets_lm <- targets |> 
  pivot_wider(names_from = 'variable', values_from = 'observation') |> 
  left_join(weather_past_daily, by = c("datetime","site_id")) |>
  arrange(site_id, datetime) |>
  group_by(site_id) |>
  mutate(
    # The 1-day lag in water temperature
    wtemp_yday = lag(temperature, 1),
    
    # The Derived Covariate: Prior 7-day mean air temperature
    prior_week_mean_airt = (lag(air_temperature, 1) + lag(air_temperature, 2) + 
                              lag(air_temperature, 3) + lag(air_temperature, 4) + 
                              lag(air_temperature, 5) + lag(air_temperature, 6) + 
                              lag(air_temperature, 7)) / 7
  ) |>
  ungroup() |>
  # Drop NAs created by the 7-day lag so the model fits cleanly
  filter(complete.cases(temperature, wtemp_yday, air_temperature, prior_week_mean_airt))

# Loop through each site to fit the model
forecast_df <- NULL

for(i in 1:length(focal_sites)) {  
  
  curr_site <- focal_sites[i]
  
  site_target <- targets_lm |> filter(site_id == curr_site)
  noaa_future_site <- weather_future_daily |> filter(site_id == curr_site)
  
  # Fit linear model based on past data
  fit <- lm(temperature ~ wtemp_yday + air_temperature + prior_week_mean_airt, data = site_target)
  fit_summary <- summary(fit)
  
  # --- UNCERTAINTY EXTRACTION ---
  coeffs <- fit$coefficients
  params_se <- fit_summary$coefficients[, 2]
  sigma <- sd(fit$residuals, na.rm = TRUE)
  
  # Assignment parameters
  n_members <- 310
  forecast_start_date <- forecast_date
  forecasted_dates <- seq(from = forecast_start_date, to = max(noaa_future_site$datetime), by = "day")
  
  # --- INITIAL CONDITIONS UNCERTAINTY ---
  curr_wt <- tail(site_target$temperature, 1) # Last known observation
  ic_sd <- 0.2 # Estimated observation error
  ic_uc <- rnorm(n = n_members, mean = curr_wt, sd = ic_sd)
  
  ic_df <- tibble(
    forecast_date = rep(forecast_start_date, times = n_members),
    ensemble_member = 1:n_members,
    forecast_variable = "temperature",
    value = ic_uc,
    uc_type = "total"
  )
  
  # --- PARAMETER UNCERTAINTY ---
  param_df <- data.frame(
    beta1 = rnorm(n_members, coeffs[1], params_se[1]),
    beta2 = rnorm(n_members, coeffs[2], params_se[2]),
    beta3 = rnorm(n_members, coeffs[3], params_se[3]),
    beta4 = rnorm(n_members, coeffs[4], params_se[4])
  )
  
  # Set up empty dataframe and insert Initial Conditions
  forecast_total_unc <- tibble(
    forecast_date = rep(forecasted_dates, times = n_members),
    ensemble_member = rep(1:n_members, each = length(forecasted_dates)),
    forecast_variable = "temperature",
    value = as.double(NA),
    uc_type = "total"
  ) |> 
    rows_update(ic_df, by = c("forecast_date","ensemble_member","forecast_variable","uc_type")) 
  
  # --- UPDATED: Setup rolling memory for derived covariate ---
  initial_past_temps <- weather_past_daily |> 
    filter(site_id == curr_site, datetime <= forecast_start_date) |> 
    tail(7) |> 
    pull(air_temperature)
  
  # Create a 310 (rows) x 7 (columns) matrix.
  recent_air_temps_matrix <- matrix(
    rep(initial_past_temps, each = n_members),
    nrow = n_members,
    ncol = 7
  )
  
  # --- FORECAST GENERATION (Mod 6 Time Loop) ---
  for(d in 2:length(forecasted_dates)) {
    
    # pull dataframes for relevant dates
    temp_pred <- forecast_total_unc |> filter(forecast_date == forecasted_dates[d])
    temp_lag <- forecast_total_unc |> filter(forecast_date == forecasted_dates[d-1])
    temp_driv <- noaa_future_site |> filter(datetime == forecasted_dates[d])
    
    # THE CAROUSEL: Expand 31 weather ensembles to 310
    temp_driv_310 <- temp_driv[rep(1:nrow(temp_driv), length.out = n_members), ]
    
    # Calculate derived covariate for today (row-wise means for all 310 members)
    curr_prior_week_mean <- rowMeans(recent_air_temps_matrix, na.rm = TRUE)
    
    # run model using param_df, temp_lag, and adding process noise
    temp_pred$value <- param_df$beta1 + 
      temp_lag$value * param_df$beta2 + 
      temp_driv_310$air_temperature * param_df$beta3 + 
      curr_prior_week_mean * param_df$beta4 + 
      rnorm(n = n_members, mean = 0, sd = sigma)
    
    # insert values back into the forecast dataframe
    forecast_total_unc <- forecast_total_unc |> 
      rows_update(temp_pred, by = c("forecast_date","ensemble_member","forecast_variable","uc_type"))
    
    # Update rolling memory for tomorrow
    # Drop the oldest day (column 1) and bind the new forecast day (column 7)
    recent_air_temps_matrix <- cbind(
      recent_air_temps_matrix[, -1], 
      temp_driv_310$air_temperature
    )
  }
  
  curr_site_df <- forecast_total_unc |>
    filter(forecast_date > forecast_start_date) |>
    rename(datetime = forecast_date, parameter = ensemble_member, prediction = value) |>
    mutate(site_id = curr_site, variable = "temperature") |>
    select(datetime, site_id, parameter, prediction, variable)
  
  forecast_df <- dplyr::bind_rows(forecast_df, curr_site_df)
  message(curr_site, ' 310-member forecast run complete')
  
}

#--------------------------#
#---- Covert to EFI standard ----

# Make forecast fit the EFI standards
forecast_df_EFI <- forecast_df %>%
  filter(datetime > forecast_date) %>%
  mutate(model_id = my_model_id,
         reference_datetime = forecast_date,
         family = 'ensemble',
         duration = 'P1D',
         parameter = as.character(parameter),
         project_id = 'neon4cast') %>%
  select(datetime, reference_datetime, duration, site_id, family, parameter, variable, prediction, model_id, project_id)
#---------------------------#



# ----- Submit forecast -----
# Write the forecast to file
theme <- 'aquatics'
date <- forecast_df_EFI$reference_datetime[1]
forecast_name <- paste0(forecast_df_EFI$model_id[1], ".csv")
forecast_file <- paste(theme, date, forecast_name, sep = '-')

write_csv(forecast_df_EFI, forecast_file)

neon4cast::forecast_output_validator(forecast_file)


neon4cast::submit(forecast_file =  forecast_file, ask = FALSE) # if ask = T (default), it will produce a pop-up box asking if you want to submit

#--------------------------#

forecast_df_EFI |> 
  ggplot(aes(x=datetime, y=prediction, group = parameter)) +
  geom_line() +
  facet_wrap(~site_id) +
  labs(title = paste0('Forecast generated for ', forecast_df_EFI$variable[1], ' on ', forecast_df_EFI$reference_datetime[1]))

plot_file_name <- paste0("Submit_forecast/", forecast_df_EFI$variable[1], '-', forecast_df_EFI$reference_datetime[1], ".png")
ggsave(plot_file_name)


# In class assignment  ####
# all_results <- duckdbfs::open_dataset("s3://bio230014-bucket01/challenges/forecasts/bundled-parquet/project_id=neon4cast/duration=P1D/variable=temperature
# ", s3_endpoint = "sdsc.osn.xsede.org", anonymous = TRUE)
# 
# my_forecast <- all_results |> 
#   filter(reference_datetime == as_datetime("2025-03-24 00:00:00"),
#          site_id == "BARC",
#          model_id %in% c("climatology", "persistenceRW")) |> 
#   collect()
# 
# single_forecast <- my_forecast |> 
#   filter( model_id == "persistenceRW",
#           datetime == as_datetime("2025-03-27 00:00:00"))
# 
# 
# 
# 
# 
# my_forecast |> 
#   filter(model_id == "climatology") |> 
#   pivot_wider(names_from = parameter, values_from = prediction ) |> 
#   mutate(lower_2.5 = mu - 1.96*sigma,
#          upper_97.5 = mu +1.96*sigma) |> 
#   ggplot(aes(x = datetime)) +
#   geom_ribbon(aes(ymin = lower_2.5, ymax= upper_97.5), fill = "lightblue")+
#   geom_line(aes(y = mu))
# 
# 
# my_forecast |> 
#   filter(model_id == "persistenceRW") |> 
#   ggplot(aes(x = datetime, y= prediciton, group = parameter)) +
#   geom_line()
# 
# 
# s3://anonymous@bio230014-bucket01/challenges/forecasts/bundled-parquet/project_id=neon4cast/duration=P1D/variable=temperature/model_id=climatology?endpoint_override=sdsc.osn.xsede.org