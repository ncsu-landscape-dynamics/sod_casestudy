library(PoPS)
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_den100m_2015.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_max100m.tif"
temp <- TRUE
precip <- FALSE
precipitation_coefficient_file <- ""
time_step <- "week"
season_month_start <- 1
season_month_end <- 12
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatments_file <- ""
treatment_method <- "all infected"
natural_kernel_type <- "exponential"
anthropogenic_kernel_type <- "cauchy"
natural_dir <- "N"
anthropogenic_dir <- "NONE"
natural_kappa <- 2
anthropogenic_kappa <- 0
pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
mask <- NULL
success_metric = "quantity and configuration"
output_frequency <- "year"
movements_file = ""  ## ignore - for pigs
use_movements = FALSE
num_iterations = 10000
number_of_cores = 40
params_to_estimate <- c(T, T, T, T, T, T)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4th: anthropogenic distance

## set data files for individual years
# infected_file_2016 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/Cumulative Infections/cum_inf_2016eu.tif"
infected_file_2017 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/Cumulative Infections/cum_inf_2017eu.tif"
infected_file_2018 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/Cumulative Infections/cum_inf_2018eu.tif"
infected_file_2019 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/Cumulative Infections/cum_inf_2019eu.tif"

infected_file_2016s <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/end_of_year/end_inf_2016eu.tif"
infected_file_2017s <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/end_of_year/end_inf_2017eu.tif"
infected_file_2018s <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/end_of_year/end_inf_2018eu.tif"

temperature_coefficient_file_2017 <- "H:/Shared drives/APHIS  Projects/shared resources/data/DaymetUS/Oregon/ForLEMMA/100m/Small/weather_coef_2017.tif"
temperature_coefficient_file_2018 <- "H:/Shared drives/APHIS  Projects/shared resources/data/DaymetUS/Oregon/ForLEMMA/100m/Small/weather_coef_2018.tif"
temperature_coefficient_file_2019 <- "H:/Shared drives/APHIS  Projects/shared resources/data/DaymetUS/Oregon/ForLEMMA/100m/Small/weather_coef_2019.tif"

posterior_means_2017 <- as.matrix(t(read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/ABC_Calibration/2016_2017_ratio_means.csv")))[1,1:6]
posterior_cov_matrix_2017 <- as.matrix(read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/ABC_Calibration/2016_2017_ratio_cov_matrix.csv"))

posterior_means_2018 <- as.matrix(t(read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/ABC_Calibration/2017_2018_ratio_means.csv")))[1,1:6]
posterior_cov_matrix_2018 <- as.matrix(read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/ABC_Calibration/2017_2018_ratio_cov_matrix.csv"))
## compare all calibrations to hindcast in 2017
start_date <- '2017-01-01'
end_date <- '2017-12-31'
treatment_dates <- c('2017-12-24')

val_2017_params2017 <- abc_validate(infected_years_file = infected_file_2017, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2017,
                                    posterior_cov_matrix= posterior_cov_matrix_2017,
                                    params_to_estimate,
                                    infected_file = infected_file_2016s, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2017, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2017_params2018 <- abc_validate(infected_years_file = infected_file_2017, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2018,
                                    posterior_cov_matrix= posterior_cov_matrix_2018,
                                    params_to_estimate,
                                    infected_file = infected_file_2016s, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2017, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2017_params2017$year <- "2017"
val_2017_params2018$year <- "2017"

val_2017_params2017$cal_year <- "2017"
val_2017_params2018$cal_year <- "2018"

## compare all calibrations to hindcast in 2018
start_date <- '2018-01-01'
end_date <- '2018-12-31'
treatment_dates <- c('2018-12-24')

val_2018_params2017 <- abc_validate(infected_years_file = infected_file_2018, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2017,
                                    posterior_cov_matrix= posterior_cov_matrix_2017,
                                    params_to_estimate,
                                    infected_file = infected_file_2017s, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2018, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2018_params2018 <- abc_validate(infected_years_file = infected_file_2018, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2018,
                                    posterior_cov_matrix= posterior_cov_matrix_2018,
                                    params_to_estimate,
                                    infected_file = infected_file_2017s, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2018, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2018_params2017$year <- "2018"
val_2018_params2018$year <- "2018"

val_2018_params2017$cal_year <- "2017"
val_2018_params2018$cal_year <- "2018"

## compare all calibrations to hindcast in 2018
start_date <- '2019-01-01'
end_date <- '2019-12-31'
treatment_dates <- c('2019-12-24')

val_2019_params2017 <- abc_validate(infected_years_file = infected_file_2019, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2017,
                                    posterior_cov_matrix= posterior_cov_matrix_2017,
                                    params_to_estimate,
                                    infected_file = infected_file_2018s, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2018, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2019_params2018 <- abc_validate(infected_years_file = infected_file_2019, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2018,
                                    posterior_cov_matrix= posterior_cov_matrix_2018,
                                    params_to_estimate,
                                    infected_file = infected_file_2018s, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2018, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2019_params2017$year <- "2019"
val_2019_params2018$year <- "2019"

val_2019_params2017$cal_year <- "2017"
val_2019_params2018$cal_year <- "2018"

## plot out statistics

val_total_params2017 <- val_2017_params2017[,1:10] + val_2018_params2017[,1:10] + val_2019_params2017[,1:10]
val_total_params2017$year <- "All"
val_total_params2017$cal_year <- "2017"

val_total_params2018 <- val_2017_params2018[,1:10] + val_2018_params2018[,1:10] + val_2019_params2018[,1:10]
val_total_params2018$year <- "All"
val_total_params2018$cal_year <- "2018"

vals <- rbind( val_2017_params2017, val_2017_params2018, val_2018_params2017, val_2018_params2018, val_2019_params2017, val_2019_params2018)
vals_total <- rbind(val_total_params2017, val_total_params2018)
vals_total$configuration_disagreement <- vals_total$configuration_disagreement/4
total_infs <- 67+88+112
vals_total$configuration_disagreement_weighted[vals_total$cal_year == "2017"] <- val_2017_params2017$configuration_disagreement * (67/total_infs) + val_2018_params2017$configuration_disagreement * (88/total_infs) + val_2019_params2017$configuration_disagreement * (112/total_infs)
vals_total$configuration_disagreement_weighted[vals_total$cal_year == "2018"] <- val_2017_params2018$configuration_disagreement * (67/total_infs) + val_2018_params2018$configuration_disagreement * (88/total_infs) + val_2019_params2018$configuration_disagreement * (112/total_infs)
# vals_total <- vals_total[vals_total$quantity_disagreement < 500,]

ggplot(vals, aes(x=year, y=quantity_disagreement, fill=cal_year)) + 
  geom_boxplot() +scale_y_continuous(limits = quantile(vals$quantity_disagreement, c(0.025, 0.975)))

ggplot(vals, aes(x=year, y=configuration_disagreement, fill=cal_year)) + 
  geom_boxplot()

ggplot(vals_total, aes(x=year, y=quantity_disagreement, fill=cal_year)) + 
  geom_boxplot() + scale_y_continuous(limits = quantile(vals_total$quantity_disagreement, c(0.01, 0.975)))

# ggplot(vals_total, aes(x=year, y=configuration_disagreement, fill=cal_year)) + 
#   geom_boxplot()

ggplot(vals_total, aes(x=year, y=configuration_disagreement_weighted, fill=cal_year)) + 
  geom_boxplot() + xlab("All years combined") + ylab("Configuration Disagreement") + 
  labs(fill ="Calibration Year") 

#check total values
inf_2019 <- raster(infected_file_2019)
sum(inf_2019[inf_2019 > 0] > 0)

inf_2018 <- raster(infected_file_2018)
sum(inf_2018[inf_2018 > 0] > 0)

inf_2017 <- raster(infected_file_2017)
sum(inf_2017[inf_2017 > 0] > 0)

inf_2016 <- raster(infected_file_2016)
sum(inf_2016[inf_2016 > 0] > 0)

inf_2015 <- raster(infected_file_2015)
sum(inf_2015[inf_2015 > 0] > 0)
