#############################################################################
#                                                                           #  
# HEALTHY LIFESPAN INEQUALITY: MORBIDITY COMPRESSION FROM                   #
#  A GLOBAL PERSPECTIVE                                                     #
#                                                                           #
# Inaki PERMANYER - Francisco VILLAVICENCIO - Sergi TRIAS-LLIMOS            #
#                                                                           #
# CODE TO ESTIMATE INEQUALITY MEASURES                                      #    
# Requires: Functions.R                                                     #
#                                                                           #
#                                                             February 2023 #
#############################################################################


############################################################
### WORKING DIRECTORY AND PACKAGES                       ###  
############################################################

# Clear working directory
rm(list = ls())

# Prevent scientific notation
options(scipen = 999)

# Packages
# install.packages('EnvStats')
library(EnvStats)               # Truncated log-normal distribution
# install.packages('msm')
library(msm)                    # Truncated normal distribution

# Functions
source('Functions.R')


############################################################
### DATA: LIFE TABLES AND HALE ESTIMATES FROM GBD 2019   ###  
############################################################

# COUNTRY AND REGIONAL CLASSIFICATION
countryClass <- read.csv('GBD2019/GBD-CountryClass.csv')
dim(countryClass)
head(countryClass)

# YEARS AND SEX-GROUPS OF INTEREST
Years <- 1990:2019
Sex <- c('Female', 'Male', 'Both')


#-------------#
# LIFE TABLES #
#-------------#

# Source: Global Burden of Disease (GBD) Study 2019 
#  https://ghdx.healthdata.org/record/ihme-data/gbd-2019-life-tables-1950-2019
#  (data downloaded on 11 July 2022)

# Data frame to store life tables
datLT <- data.frame()

# Upload data
for (i in c(5:20, 28, 30:33, 44, 45, 148)) {
  
  # Auxiliary data set
  dat1 <- read.csv(paste0('GBD2019/LifeTables/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_',
                          i, '_WSHOCK_Y2020M11D13.CSV'))
  
  # Select countries
  dat1 <- dat1[dat1$location_id %in% unique(c(1, countryClass$location_id, 
                                              countryClass$super_region_id,
                                              countryClass$region_id)), ]
  
  # Select years
  names(dat1)[names(dat1) == 'year_id'] <- 'year'
  dat1 <- dat1[dat1$year %in% Years, ]
  
  # Select sex
  names(dat1)[names(dat1) == 'sex_name'] <- 'sex'
  dat1$sex[dat1$sex == 'female'] <- Sex[1]
  dat1$sex[dat1$sex == 'male'] <- Sex[2]
  dat1$sex[dat1$sex == 'both'] <- Sex[3]
  dat1 <- dat1[dat1$sex %in% Sex, ]
  
  # Select variables
  dat1 <- dat1[, c('location_id', 'location_name', 'sex', 'age_group_name', 
                   'year', 'measure_name', 'lower', 'val', 'upper')]
  
  # Add data set
  datLT <- rbind(datLT, dat1)
  rm(dat1, i)
  
}

# Re-label ages
names(datLT)[names(datLT) == 'age_group_name'] <- 'age'
datLT$age[datLT$age == '<1 year'] <- 0
datLT$age[datLT$age == '1 to 4'] <- 1
datLT$age[datLT$age == '5 to 9'] <- 5
datLT$age[datLT$age == '100 to 104'] <- 100
datLT$age[datLT$age == '105 to 109'] <- 105
datLT$age[datLT$age == '110 plus'] <- 110
datLT$age[!datLT$age %in% c('<1 year', '1 to 4', '5 to 9', '100', '105', '110')] <- 
  substr(datLT$age[!datLT$age %in% c('<1 year', '1 to 4', '5 to 9', '100', '105', '110')], 1, 2)
datLT$age <- as.numeric(datLT$age)

# Sex
datLT$sex <- factor(datLT$sex, levels = Sex)

# Add ISO3 code
datLT <- merge(datLT, countryClass[, c('ISO3', 'location_id')],
               by = 'location_id', all.x = T)
datLT <- datLT[, c(1, 2, ncol(datLT), 3:(ncol(datLT) - 1))]

# Tidy up
dat1 <- datLT[which(is.na(datLT$ISO3)), ]
dat2 <- datLT[which(!is.na(datLT$ISO3)), ]
datLT <- rbind(dat1[order(dat1$location_id, dat1$year, dat1$sex, dat1$age), ], 
               dat2[order(dat2$ISO3, dat2$year, dat2$sex, dat2$age), ])
rm(dat1, dat2)
rownames(datLT) <- NULL
head(datLT)


#----------------------------------------------#
# AGE-SPECIFIC HEALTH-ADJUSTED LIFE EXPECTANCY #
#----------------------------------------------#

# Source: Global Burden of Disease (GBD) Study 2019 
#  http://ghdx.healthdata.org/gbd-results-tool 
#  (data downloaded on 11 July 2022)

# Upload data
datHALE <- rbind(read.csv(file = 'GBD2019/HALE/IHME-GBD_2019_DATA-ebdcb389-1.csv'),
                 read.csv(file = 'GBD2019/HALE/IHME-GBD_2019_DATA-ebdcb389-2.csv'),
                 read.csv(file = 'GBD2019/HALE/IHME-GBD_2019_DATA-ebdcb389-3.csv'))
datHALE <- datHALE[, c('location_id', 'location_name', 'sex_name', 'age_name', 
                       'year', 'lower', 'val', 'upper')]

# Select super-regions, regions, and countries
datHALE <- datHALE[datHALE$location_id %in% unique(c(1, countryClass$location_id, 
                                                     countryClass$super_region_id,
                                                     countryClass$region_id)), ]

# Add ISO3 code
datHALE <- merge(datHALE, countryClass[, c('ISO3', 'location_id')],
                 by = 'location_id', all.x = T)
datHALE <- datHALE[, c(1, 2, ncol(datHALE), 3:(ncol(datHALE) - 1))]

# Re-label ages
names(datHALE)[names(datHALE) == 'age_name'] <- 'age'
datHALE <- datHALE[datHALE$age != 'All ages', ]
datHALE$age[datHALE$age == '<1 year'] <- 0
datHALE$age[datHALE$age == '1-4 years'] <- 1
datHALE$age[datHALE$age == '5-9 years'] <- 5
datHALE$age[!datHALE$age %in% c(0, 1, 5)] <- 
  substr(datHALE$age[!datHALE$age %in% c(0, 1, 5)], 1, 2)
datHALE$age <- as.numeric(datHALE$age)
table(datHALE$age)

# Sex
names(datHALE)[names(datHALE) == 'sex_name'] <- 'sex'
datHALE <- datHALE[datHALE$sex %in% Sex, ]
datHALE$sex <- factor(datHALE$sex, levels = Sex)

# Tidy up
dat1 <- datHALE[which(is.na(datHALE$ISO3)), ]
dat2 <- datHALE[which(!is.na(datHALE$ISO3)), ]
datHALE <- rbind(dat1[order(dat1$location_id, dat1$year, dat1$sex, dat1$age), ], 
                 dat2[order(dat2$ISO3, dat2$year, dat2$sex, dat2$age), ])
rm(dat1, dat2)
rownames(datHALE) <- NULL
head(datHALE)

# SAVE
SAVE <- F
if (SAVE) {
  
  # Life tables
  write.csv(datLT, row.names = F,
            paste0('GBD2019/', format(Sys.Date(), "%Y%m%d"), 
                   '-GBD-LT-', min(Years), '-', max(Years), '.csv'))
  
  # Health-adjusted life expectancy
  write.csv(datHALE, row.names = F,
            paste0('GBD2019/', format(Sys.Date(), "%Y%m%d"), 
                   '-GBD-HALE-', min(Years), '-', max(Years), '.csv'))
  
}


############################################################
### LIFESPAN INEQUALITY AND HEALTHY LIFESPAN INEQUALITY  ###  
############################################################

#---------------# 
# RATIO ax / ex #
#---------------# 

# Object to store ax/ex ratios
rho <- c()

# Calculate ax/ex ratios at older ages
for (age in seq(90, 105, 5)) {
  
  # Life expectancy at age x
  ex1 <- datLT$val[datLT$age == age & datLT$measure_name == 'Life expectancy']
  
  # Life expectancy at age x+n
  ex2 <- datLT$val[datLT$age == age + 5 & datLT$measure_name == 'Life expectancy']
  
  # Probability of death at age x
  qx <- datLT$val[datLT$age == age & datLT$measure_name == 'Probability of death']
  
  # Estimate ax
  ax <- (ex1 - (ex2 + 5)*(1-qx)) / qx
  
  # Average ratio
  rho <- c(rho, round(mean(ax/ex1), 3))
  
}
rm(age, ax, ex1, ex2, qx)


#---------------------#
# INEQUALITY MEASURES #
#---------------------#

# Locations
locations <- unique(datLT$location_id)

# Ages for which we calculate inequality measures
Ages <- c(0, 65)

# Calculate uncertainty intervals (UI)? (T/F)
UI <- T
if (UI) {
  # Uncertainty interval to report (by default set to 80%)
  UI <- .8
  # Number of random draws to calculate UI
  ndraw <- 10000
} else {
  UI <- NULL
  ndraw <- NULL
}

# Parallel computing? (T/F)
PARALLEL <- T

#---------------------------------------------------------------------#
# NOTE: Expected execution time on a 4 CPU laptop                     #
#                                                                     #
# Parallel computing without uncertainty: less than 1 min             #
# No parallel computing without uncertainty: 30 min                   #
# Parallel computing with uncertainty: around 2.5 hours               #
# No parallel computing with uncertainty: 10+ hours (not recommended) #
#---------------------------------------------------------------------#

# Time control
Start <- Sys.time()

# Calculate lifespan inequality measures
if (PARALLEL) {
  
  # Upload necessary packages
  # install.packages('doParallel')
  library(doParallel)
  
  # Number of cores
  ncores <- min(20, detectCores() - 1, length(locations))
  
  # Start cluster
  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  # Parallel computing
  out <- foreach(i = 1:length(locations), .packages = c('msm', 'EnvStats')) %dopar% {
    CalcHLI(id = locations[i], 
            lt = datLT[datLT$location_id == locations[i], ], 
            hale = datHALE[datHALE$location_id == locations[i], ], 
            countries = countryClass, years = Years, ages = Ages,
            sexGr = Sex, UI = UI, n = ndraw)
  }
  
  # Stop cluster
  stopCluster(cl)

  # Re-arrange data
  results <- data.frame()
  ratios <- data.frame()
  for (i in 1:length(out)) {
    results <- rbind(results, out[[i]]$results)
    ratios <- rbind(ratios, out[[i]]$ratios)
  }
  
} else {
  
  # No parallel computing
  out <- CalcHLI(id = locations, lt = datLT, hale = datHALE, 
                 countries = countryClass, years = Years,
                 ages = Ages, sexGr = Sex, UI = UI, n = ndraw)
  results <- out$results
  ratios <- out$ratios
  
}

# Time control
End <- Sys.time()
print(End - Start)

# Save estimates
SAVE <- F
if (SAVE) {
  
  # LIFE EXPECTANCY AND INEQUALITY INDICATORS
  write.csv(results, row.names = F, 
            file = paste0('Results/EstimatesHLI-LI-', 
                          min(Years), '-', max(Years), '.csv'))
  
  # HLI / LI RATIOS
  write.csv(ratios, row.names = F, 
            file = paste0('Results/RatiosHLI-LI-', 
                          min(Years), '-', max(Years), '.csv'))
  
}

