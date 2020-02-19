
## CALIBRATION FOR EU1

## number of infections in 2016: 8
## number of infections in 2017: 67
## number of infections in 2018: 75

library(PoPS)
# source('C:/Users/dagaydos/Documents/GitHub/rpops/R/uncertainty_propogation.R')
# source('C:/Users/dagaydos/Documents/GitHub/rpops/R/helpers.R')
# source('C:/Users/dagaydos/Documents/GitHub/rpops/R/checks.R')
# source("C:/Users/dagaydos/Documents/GitHub/rpops/R/abc_calibration.R")

infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/Cumulative Infections/cum_inf_2017eu.tif"
number_of_observations <- 67     ### This is the number of infected cells - just make sure it's consistent across years 
prior_number_of_observations <- 0
prior_means <- c(0, 0, 0, 0, 0, 0)    ### leave as 0 for now, means that you are giving them no weight
prior_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
params_to_estimate <- c(T, T, T, T, T, T)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4th: anthropogenic distance
number_of_generations <- 8
generation_size <- 1000
checks = c(100, 11000, 1000, 10000)   ### 1: difference between number of infected cells in different years * 2
                                      ### 2: difference in distance in different years *2
                                      ### 3 and 4 don't matter if you keep the metric of success the same
infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/Cumulative Infections/cum_inf_2016eu.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_den100m_2015.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_max100m.tif"
temp <- TRUE
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/shared resources/data/DaymetUS/Oregon/ForLEMMA/100m/Small/weather_coef_2017.tif"
precip <- FALSE
precipitation_coefficient_file <- ""
time_step <- "week"
season_month_start <- 1
season_month_end <- 12
start_date <- '2017-01-01'
end_date <- '2017-12-31'
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2017-12-24')
treatments_file <- ""
treatment_method <- "all infected"
natural_kernel_type <- "exponential"
anthropogenic_kernel_type <- "cauchy"
natural_dir <- "N"
natural_kappa <- 2
anthropogenic_dir <- "NONE"
anthropogenic_kappa <- 0
pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
mask <- NULL
success_metric <- "number of locations and total distance"  ### keep this the same
output_frequency <- "year"
movements_file = ""  ## ignore - for pigs
use_movements = FALSE

eu_sod_2016_2017 <- abc_calibration(infected_years_file, 
                                     number_of_observations, prior_number_of_observations,
                                     prior_means, prior_cov_matrix, params_to_estimate,
                                     number_of_generations,
                                     generation_size,
                                     checks,
                                     infected_file, host_file, total_plants_file, 
                                     temp, temperature_coefficient_file, 
                                     precip, precipitation_coefficient_file, 
                                     time_step, 
                                     season_month_start, season_month_end, 
                                     start_date, end_date, 
                                     use_lethal_temperature, temperature_file,
                                     lethal_temperature, lethal_temperature_month,
                                     mortality_on, mortality_rate, mortality_time_lag, 
                                     management, treatment_dates, treatments_file,
                                     treatment_method,
                                     natural_kernel_type, anthropogenic_kernel_type,
                                     natural_dir, natural_kappa, 
                                     anthropogenic_dir, anthropogenic_kappa,
                                     pesticide_duration, pesticide_efficacy,
                                     mask, success_metric, output_frequency,
                                     movements_file, use_movements)


infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/Cumulative Infections/cum_inf_2018eu.tif"
number_of_observations <- 75
prior_number_of_observations <- eu_sod_2016_2017$total_number_of_observations
prior_means <- eu_sod_2016_2017$posterior_means
prior_cov_matrix <- eu_sod_2016_2017$posterior_cov_matrix
params_to_estimate <- c(T, T, T, T, T, T)
number_of_generations <- 8
generation_size <- 1000
checks = c(100,10000, 1000, 10000)
infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/Cumulative Infections/cum_inf_2017eu.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_den100m_2015.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_max100m.tif"
temp <- TRUE
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/shared resources/data/DaymetUS/Oregon/ForLEMMA/100m/Small/weather_coef_2017.tif"
precip <- FALSEC
precipitation_coefficient_file <- ""
time_step <- "week"
season_month_start <- 1
season_month_end <- 12
start_date <- '2018-01-01'
end_date <- '2018-12-31'
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2018-12-24')
treatments_file <- ""
treatment_method <- "ratio"
natural_kernel_type <- "cauchy"
anthropogenic_kernel_type <- "cauchy"
natural_dir <- "NONE"
natural_kappa <- 0
anthropogenic_dir <- "NONE"
anthropogenic_kappa <- 0
pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
mask <- NULL
success_metric <- "number of locations and total distance"
output_frequency <- "year"
movements_file = ""
use_movements = FALSE


eu_sod_2017_2018 <- abc_calibration(infected_years_file, 
                                     number_of_observations, prior_number_of_observations,
                                     prior_means, prior_cov_matrix, params_to_estimate,
                                     number_of_generations,
                                     generation_size,
                                     checks,
                                     infected_file, host_file, total_plants_file, 
                                     temp, temperature_coefficient_file, 
                                     precip, precipitation_coefficient_file, 
                                     time_step, 
                                     season_month_start, season_month_end, 
                                     start_date, end_date, 
                                     use_lethal_temperature, temperature_file,
                                     lethal_temperature, lethal_temperature_month,
                                     mortality_on, mortality_rate, mortality_time_lag, 
                                     management, treatment_dates, treatments_file,
                                     treatment_method,
                                     natural_kernel_type, anthropogenic_kernel_type,
                                     natural_dir, natural_kappa, 
                                     anthropogenic_dir, anthropogenic_kappa,
                                     pesticide_duration, pesticide_efficacy,
                                     mask, success_metric, output_frequency,
                                     movements_file, use_movements)





