
## CALIBRATION FOR NA1  - USING RATIO MANAGEMENT

## number of infections end of year 2001: 13
## number of infections cumulative 2002: 36 
## diff over 2002: 23

## number of infections end of year 2002: 36
## number of infections cumulative 2003: 58
## diff over 2003: 22

## number of infections end of year 2003: 58
## number of infections cumulative 2004: 73
## diff over 2004: 15

## number of infections end of year 2004: 73
## number of infections cumulative 2005: 91
## diff over 2005: 18

## number of infections end of year 2005: 78
## number of infections cumulative 2006: 142
## diff over 2006: 64

## number of infections end of year 2006: 102
## number of infections cumulative 2007: 181
## diff over 2007: 79

## number of infections end of year 2007: 152
## number of infections cumulative 2008: 228
## diff over 2008: 76

## number of infections end of year 2008: 201
## number of infections cumulative 2009: 276
## diff over 2009: 75

## number of infections end of year 2009: 256
## number of infections cumulative 2010: 342
## diff over 2010: 86


library(PoPS)
# source('C:/Users/dagaydos/Documents/GitHub/rpops/R/uncertainty_propogation.R')
# source('C:/Users/dagaydos/Documents/GitHub/rpops/R/helpers.R')
# source('C:/Users/dagaydos/Documents/GitHub/rpops/R/checks.R')
# source("C:/Users/dagaydos/Documents/GitHub/rpops/R/abc_calibration.R")


### 2001 to 2002

infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/Cumulative Infections/cum_inf_2002.tif"
number_of_observations <- ### ADDDDDDDDDDD     ### This is the number of infected cells - just make sure it's consistent across years 
  prior_number_of_observations <- 0
prior_means <- c(0, 0, 0, 0, 0, 0)    ### leave as 0 for now, means that you are giving them no weight
prior_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
params_to_estimate <- c(T, T, T, T, T, T)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4th: anthropogenic distance
number_of_generations <- 8
generation_size <- 1000
checks = c(60,20000,1000,10000)   #### CHECKED USED FOR EU1 : c(60,20000,1000,10000)
infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/end_of_year/end_inf_2001.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_den100m.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_max100m.tif"
temp <- TRUE
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/shared resources/data/DaymetUS/Oregon/ForLEMMA/100m/Small/weather_coef_2002.tif"
precip <- FALSE
precipitation_coefficient_file <- ""
time_step <- "week"
season_month_start <- 1
season_month_end <- 12
start_date <- '2002-01-01'
end_date <- '2002-12-31'
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2016-12-24')
treatments_file <- ""
treatment_method <- "ratio"
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

na_sod_2001_2002 <- abc_calibration(infected_years_file, 
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


### 2002 to 2003

infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/Cumulative Infections/cum_inf_2003.tif"
number_of_observations <- ### ADDDDDDDDDDD     ### This is the number of infected cells - just make sure it's consistent across years 
  prior_number_of_observations <- na_sod_2001_2002$total_number_of_observations
prior_means <- na_sod_2001_2002$posterior_means
prior_cov_matrix <- na_sod_2001_2002$posterior_cov_matrix
params_to_estimate <- c(T, T, T, T, T, T)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4th: anthropogenic distance
number_of_generations <- 8
generation_size <- 1000
checks = c(60,20000,1000,10000)  #### CHECKED USED FOR EU1 : c(60,20000,1000,10000)
infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/end_of_year/end_inf_2002.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_den100m.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_max100m.tif"
temp <- TRUE
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/shared resources/data/DaymetUS/Oregon/ForLEMMA/100m/Small/weather_coef_2003.tif"
precip <- FALSE
precipitation_coefficient_file <- ""
time_step <- "week"
season_month_start <- 1
season_month_end <- 12
start_date <- '2003-01-01'
end_date <- '2003-12-31'
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2003-12-24')
treatments_file <- ""
treatment_method <- "ratio"
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

na_sod_2002_2003 <- abc_calibration(infected_years_file, 
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



### 2003 to 2004

infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/Cumulative Infections/cum_inf_2004.tif"
number_of_observations <- ### ADDDDDDDDDDD     ### This is the number of infected cells - just make sure it's consistent across years 
  prior_number_of_observations <- na_sod_2002_2003$total_number_of_observations
prior_means <- na_sod_2002_2003$posterior_means
prior_cov_matrix <- na_sod_2002_2003$posterior_cov_matrix
params_to_estimate <- c(T, T, T, T, T, T)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4th: anthropogenic distance
number_of_generations <- 8
generation_size <- 1000
checks = c(60,20000,1000,10000)   #### CHECKED USED FOR EU1 : c(60,20000,1000,10000)
infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/end_of_year/end_inf_2003.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_den100m.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_max100m.tif"
temp <- TRUE
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/shared resources/data/DaymetUS/Oregon/ForLEMMA/100m/Small/weather_coef_2004.tif"
precip <- FALSE
precipitation_coefficient_file <- ""
time_step <- "week"
season_month_start <- 1
season_month_end <- 12
start_date <- '2004-01-01'
end_date <- '2004-12-31'
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2004-12-24')
treatments_file <- ""
treatment_method <- "ratio"
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

na_sod_2003_2004 <- abc_calibration(infected_years_file, 
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



### 2004 to 2005

infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/Cumulative Infections/cum_inf_2005.tif"
number_of_observations <- ### ADDDDDDDDDDD     ### This is the number of infected cells - just make sure it's consistent across years 
  prior_number_of_observations <- na_sod_2003_2004$total_number_of_observations
prior_means <- na_sod_2003_2004$posterior_means
prior_cov_matrix <- na_sod_2003_2004$posterior_cov_matrix
params_to_estimate <- c(T, T, T, T, T, T)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4th: anthropogenic distance
number_of_generations <- 8
generation_size <- 1000
checks = c(60,20000,1000,10000)   #### CHECKED USED FOR EU1 : c(60,20000,1000,10000)
infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/end_of_year/end_inf_2004.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_den100m.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_max100m.tif"
temp <- TRUE
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/shared resources/data/DaymetUS/Oregon/ForLEMMA/100m/Small/weather_coef_2005.tif"
precip <- FALSE
precipitation_coefficient_file <- ""
time_step <- "week"
season_month_start <- 1
season_month_end <- 12
start_date <- '2005-01-01'
end_date <- '2005-12-31'
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2005-12-24')
treatments_file <- ""
treatment_method <- "ratio"
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

na_sod_2004_2005 <- abc_calibration(infected_years_file, 
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



### 2005 to 2006

infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/Cumulative Infections/cum_inf_2006.tif"
number_of_observations <- ### ADDDDDDDDDDD     ### This is the number of infected cells - just make sure it's consistent across years 
  prior_number_of_observations <- na_sod_2004_2005$total_number_of_observations
prior_means <- na_sod_2004_2005$posterior_means
prior_cov_matrix <- na_sod_2004_2005$posterior_cov_matrix
params_to_estimate <- c(T, T, T, T, T, T)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4th: anthropogenic distance
number_of_generations <- 8
generation_size <- 1000
checks = c(60,20000,1000,10000)   #### CHECKES USED FOR EU1 : c(60,20000,1000,10000)
infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/end_of_year/end_inf_2005.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_den100m.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_max100m.tif"
temp <- TRUE
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/shared resources/data/DaymetUS/Oregon/ForLEMMA/100m/Small/weather_coef_2006.tif"
precip <- FALSE
precipitation_coefficient_file <- ""
time_step <- "week"
season_month_start <- 1
season_month_end <- 12
start_date <- '2006-01-01'
end_date <- '2006-12-31'
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2006-12-24')
treatments_file <- ""
treatment_method <- "ratio"
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

na_sod_2005_2006 <- abc_calibration(infected_years_file, 
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

### 2006 to 2007

infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/Cumulative Infections/cum_inf_2007.tif"
number_of_observations <- ### ADDDDDDDDDDD     ### This is the number of infected cells - just make sure it's consistent across years 
  prior_number_of_observations <- na_sod_2005_2006$total_number_of_observations
prior_means <- na_sod_2005_2006$posterior_means
prior_cov_matrix <- na_sod_2005_2006$posterior_cov_matrix
params_to_estimate <- c(T, T, T, T, T, T)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4th: anthropogenic distance
number_of_generations <- 8
generation_size <- 1000
checks = c(60,20000,1000,10000)   #### CHECKES USED FOR EU1 : c(60,20000,1000,10000)
infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/end_of_year/end_inf_2006.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_den100m.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_max100m.tif"
temp <- TRUE
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/shared resources/data/DaymetUS/Oregon/ForLEMMA/100m/Small/weather_coef_2007.tif"
precip <- FALSE
precipitation_coefficient_file <- ""
time_step <- "week"
season_month_start <- 1
season_month_end <- 12
start_date <- '2007-01-01'
end_date <- '2007-12-31'
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2007-12-24')
treatments_file <- ""
treatment_method <- "ratio"
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

na_sod_2006_2007 <- abc_calibration(infected_years_file, 
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



### 2007 to 2008

infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/Cumulative Infections/cum_inf_2008.tif"
number_of_observations <- ### ADDDDDDDDDDD     ### This is the number of infected cells - just make sure it's consistent across years 
  prior_number_of_observations <- na_sod_2006_2007$total_number_of_observations
prior_means <- na_sod_2006_2007$posterior_means
prior_cov_matrix <- na_sod_2006_2007$posterior_cov_matrix
params_to_estimate <- c(T, T, T, T, T, T)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4th: anthropogenic distance
number_of_generations <- 8
generation_size <- 1000
checks = c(60,20000,1000,10000)   #### CHECKES USED FOR EU1 : c(60,20000,1000,10000)
infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/end_of_year/end_inf_2007.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_den100m.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_max100m.tif"
temp <- TRUE
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/shared resources/data/DaymetUS/Oregon/ForLEMMA/100m/Small/weather_coef_2008.tif"
precip <- FALSE
precipitation_coefficient_file <- ""
time_step <- "week"
season_month_start <- 1
season_month_end <- 12
start_date <- '2008-01-01'
end_date <- '2008-12-31'
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2008-12-24')
treatments_file <- ""
treatment_method <- "ratio"
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

na_sod_2007_2008 <- abc_calibration(infected_years_file, 
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


### 2008 to 2009

infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/Cumulative Infections/cum_inf_2009.tif"
number_of_observations <- ### ADDDDDDDDDDD     ### This is the number of infected cells - just make sure it's consistent across years 
  prior_number_of_observations <- na_sod_2007_2008$total_number_of_observations
prior_means <- na_sod_2007_2008$posterior_means
prior_cov_matrix <- na_sod_2007_2008$posterior_cov_matrix
params_to_estimate <- c(T, T, T, T, T, T)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4th: anthropogenic distance
number_of_generations <- 8
generation_size <- 1000
checks = c(60,20000,1000,10000)  #### CHECKES USED FOR EU1 : c(60,20000,1000,10000)
infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/end_of_year/end_inf_2008.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_den100m.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_max100m.tif"
temp <- TRUE
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/shared resources/data/DaymetUS/Oregon/ForLEMMA/100m/Small/weather_coef_2009.tif"
precip <- FALSE
precipitation_coefficient_file <- ""
time_step <- "week"
season_month_start <- 1
season_month_end <- 12
start_date <- '2009-01-01'
end_date <- '2009-12-31'
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2009-12-24')
treatments_file <- ""
treatment_method <- "ratio"
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

na_sod_2008_2009 <- abc_calibration(infected_years_file, 
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


### 2009 to 2010

infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/Cumulative Infections/cum_inf_2010.tif"
number_of_observations <- ### ADDDDDDDDDDD     ### This is the number of infected cells - just make sure it's consistent across years 
  prior_number_of_observations <- na_sod_2008_2009$total_number_of_observations
prior_means <- na_sod_2008_2009$posterior_means
prior_cov_matrix <- na_sod_2008_2009$posterior_cov_matrix
params_to_estimate <- c(T, T, T, T, T, T)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4th: anthropogenic distance
number_of_generations <- 8
generation_size <- 1000
checks = c(60,20000,1000,10000)   #### CHECKES USED FOR EU1 : c(60,20000,1000,10000)
infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/100mLEMMA/Small/ratio/end_of_year/end_inf_2009.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_den100m.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Host Data - Proportion/SmallerLEMMA/lemma_max100m.tif"
temp <- TRUE
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/shared resources/data/DaymetUS/Oregon/ForLEMMA/100m/Small/weather_coef_2010.tif"
precip <- FALSE
precipitation_coefficient_file <- ""
time_step <- "week"
season_month_start <- 1
season_month_end <- 12
start_date <- '2010-01-01'
end_date <- '2010-12-31'
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2010-12-24')
treatments_file <- ""
treatment_method <- "ratio"
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

na_sod_2009_2010 <- abc_calibration(infected_years_file, 
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


#### IF WE WANT TO CONTINUE AFTER 2010, WE WILL NEED TO MASK OUT GIA INFECTIONS






