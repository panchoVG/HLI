#############################################################################
#                                                                           #  
# HEALTHY LIFESPAN INEQUALITY: MORBIDITY COMPRESSION FROM                   #
#  A GLOBAL PERSPECTIVE                                                     #
#                                                                           #
# Inaki PERMANYER - Francisco VILLAVICENCIO - Sergi TRIAS-LLIMOS            #
#                                                                           #
# CODE TO RECONSTRUCT LIFE TABLES AND ESTIMATE INEQUALITY MEASURES          #    
# Requires: Functions.R                                                     #
#                                                                           #
#                                                                 July 2022 #
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
### DATA: LIFE TABLES AND HALE ESTIMATSE FROM GBD 2019   ###  
############################################################

# COUNTRY AND REGIONAL CLASSIFCIATION #
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

# Data frame to store data
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
  dat1$sex[dat1$sex == 'female'] <- 'Female'
  dat1$sex[dat1$sex == 'male'] <- 'Male'
  dat1$sex[dat1$sex == 'both'] <- 'Both'
  dat1 <- dat1[dat1$sex %in% Sex, ]
  
  # Select variables
  dat1 <- dat1[, c("location_id", "location_name", "sex", "age_group_name", 
                   "year", "measure_name", "lower", "val", "upper")]
  
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
datLT$sex <- factor(datLT$sex, levels = c('Female', 'Male', 'Both'))

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
#  http://ghdx.healthdata.org/gbd-results-tool (data downloaded on 11 July 2022)

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
datHALE$sex <- factor(datHALE$sex, levels = c('Female', 'Male', 'Both'))

# Tidy up
dat1 <- datHALE[which(is.na(datHALE$ISO3)), ]
dat2 <- datHALE[which(!is.na(datHALE$ISO3)), ]
datHALE <- rbind(dat1[order(dat1$location_id, dat1$year, dat1$sex, dat1$age), ], 
                 dat2[order(dat2$ISO3, dat2$year, dat2$sex, dat2$age), ])
rm(dat1, dat2)

datHALE <- datHALE[order(datHALE$ISO3, datHALE$location_id, 
                         datHALE$year, datHALE$sex, datHALE$age), ]
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
### BUILD LIFE TABLES AND CALCULATE HLI                  ###  
############################################################

# Data frame to store LI and HLI estimates
results <- data.frame()

# Data frame to store HLI / LI ratios
ratios <- data.frame()

# Ages for which we calculate inequality measures
Ages <- c(0, 65)

# Number of random draws
ndraw <- 10000

# Save life tables (TRUE/FALSE)? (FALSE reduces computation time substantially)
saveLT <- F
if (saveLT) {
  ltMortality <- data.frame()
  ltMorbidity <- data.frame()
}


#---------------# 
# RATIO ax / ex #
#---------------# 

# Object to store ratio ax / ex ratios
rho <- c()

# Calculate ax / ex ratio at older ages
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


#-------------------------------------------------------------#
# MORTALITY AND MORBIDITY LIFE TABLES AND INEQUALITY MEASURES #
#-------------------------------------------------------------#

# Time control
Start <- Sys.time()

# BUILD SEX-SPECIFIC LIFE TABLES FOR ALL COUNTRY-YEARS 
for (loc in unique(datLT$location_id)) {
  
  # IDENTIFYING VARIABLES
  if (loc %in% c(1, countryClass$super_region_id)) {
    if (loc == 1) {
      locName <- 'Global'
      type <- 'Global'
      supReg <- NA
    } else {
      locName <- unique(countryClass$super_region_name[countryClass$super_region_id == loc])
      type <- 'Super-region'
      supReg <- locName
    }
    iso3 <- NA
    reg <- NA
  }
  if (loc %in% countryClass$region_id) {
    locName <- unique(countryClass$region_name[countryClass$region_id == loc])
    type <- 'Region'
    # North Africa and South Asia
    if (loc %in% c(137, 158)) type <- 'Super-region/Region'
    iso3 <- NA
    supReg <- unique(countryClass$super_region_name[countryClass$region_name == locName])
    reg <- locName
  }
  if (loc %in% countryClass$location_id) {
    locName <- unique(countryClass$location_name[countryClass$location_id == loc])
    type <- 'Country/Territory'
    iso3 <- countryClass$ISO3[countryClass$location_id == loc]
    supReg <- countryClass$super_region_name[countryClass$location_id == loc]
    reg <- countryClass$region_name[countryClass$location_id == loc]
  }
  
  # LIFE TABLES FOR ALL YEARS AND SEX GROUPS
  for (year in Years) {
    
    for (sex in Sex) {
      
      # Vector of ages
      x <- sort(unique(datLT$age))
      nage <- length(x)
      
      # Death probabilities input data
      qx <- datLT$val[datLT$location_id == loc & datLT$year == year &
                        datLT$sex == sex & datLT$measure_name == 'Probability of death']
      qxLow <- datLT$lower[datLT$location_id == loc & datLT$year == year &
                             datLT$sex == sex & datLT$measure_name == 'Probability of death']
      qxUp <- datLT$upper[datLT$location_id == loc & datLT$year == year &
                            datLT$sex == sex & datLT$measure_name == 'Probability of death']
      
      # Amend inconsistencies in GBD input data
      idAdj <- which(qx < qxLow)
      if (length(idAdj) > 0) {
        qxAux <- qx[idAdj]
        qx[idAdj] <- qxLow[idAdj]
        qxLow[idAdj] <- qxAux
        rm(qxAux)
      }
      idAdj <- which(qxUp < qx)
      if (length(idAdj) > 0) {
        qxAux <- qx[idAdj]
        qx[idAdj] <- qxUp[idAdj]
        qxUp[idAdj] <- qxAux
        rm(qxAux)
      }
      idAdj <- which(round(qx[-length(qx)], 6) <= round(qxLow[-length(qxLow)], 6) | 
                       round(qx[-length(qx)], 6) >= round(qxUp[-length(qxUp)], 6))
      if (length(idAdj) > 0) qx[idAdj] <- (qxLow[idAdj] + qxUp[idAdj]) / 2
      rm(idAdj)
      
      # Death probabilities (Random draws)
      qxMat <- t(apply(cbind(qx[-nage], 
                             (qxUp[-nage] - qxLow[-nage]) / 3.92), 1,
                       function(z) {
                         rtnorm(ndraw, mean = z[1], sd = z[2], lower = 0, upper = 1)   
                       }))
      qxMat <- rbind(qxMat, 1)
      
      # Life expectancy input data
      ex <- datLT$val[datLT$location_id == loc & datLT$year == year &
                        datLT$sex == sex & datLT$measure_name == 'Life expectancy']
      exLow <- datLT$lower[datLT$location_id == loc & datLT$year == year &
                             datLT$sex == sex & datLT$measure_name == 'Life expectancy']
      exUp <- datLT$upper[datLT$location_id == loc & datLT$year == year &
                            datLT$sex == sex & datLT$measure_name == 'Life expectancy']
      
      # Amend inconsistencies in GBD input data
      idAdj <- which(ex < exLow)
      if (length(idAdj) > 0) {
        exAux <- ex[idAdj]
        ex[idAdj] <- exLow[idAdj]
        exLow[idAdj] <- exAux
        rm(exAux)
      }
      idAdj <- which(exUp < ex)
      if (length(idAdj) > 0) {
        exAux <- ex[idAdj]
        ex[idAdj] <- exUp[idAdj]
        exUp[idAdj] <- exAux
        rm(exAux)
      }
      idAdj <- which(round(ex, 6) <= round(exLow, 6) | 
                       round(ex, 6) >= round(exUp, 6))
      if (length(idAdj) > 0) ex[idAdj] <- (exLow[idAdj] + exUp[idAdj]) / 2
      rm(idAdj)
      
      # Life expectancy (Random draws)
      exMat <- t(apply(cbind(ex, (exUp - exLow) / 3.92), 1,
                       function(z) {
                         rtnorm(ndraw, mean = z[1], sd = z[2], lower = 0)   
                       }))
      
      # Average person-years lived by those dying between x and x+n
      ax <- CalcAx(x, qx, ex)
      
      # Two types of life tables
      for (measure in c('Mortality', 'Morbidity')) {
        
        # Morbidity
        if (measure == 'Morbidity') {
          
          # Store existing LI estimates to calculate HLI / LI ratio
          LI <- sdDx
          LImat <- sdMat
        
          # Health-adjusted life expectancy input data
          ex <- datHALE[datHALE$location_id == loc & datHALE$year == year & 
                          datHALE$sex == sex, 'val']
          exLow <- datHALE[datHALE$location_id == loc & datHALE$year == year &
                             datHALE$sex == sex, 'lower']
          exUp <- datHALE[datHALE$location_id == loc & datHALE$year == year &
                            datHALE$sex == sex, 'upper']
          
          # Amend inconsistencies in GBD input data
          idAdj <- which(ex < exLow)
          if (length(idAdj) > 0) {
            exAux <- ex[idAdj]
            ex[idAdj] <- exLow[idAdj]
            exLow[idAdj] <- exAux
            rm(exAux)
          }
          idAdj <- which(exUp < ex)
          if (length(idAdj) > 0) {
            exAux <- ex[idAdj]
            ex[idAdj] <- exUp[idAdj]
            exUp[idAdj] <- exAux
            rm(exAux)
          }
          idAdj <- which(round(ex, 6) <= round(exLow, 6) | 
                           round(ex, 6) >= round(exUp, 6))
          if (length(idAdj) > 0) ex[idAdj] <- (exLow[idAdj] + exUp[idAdj]) / 2
          rm(idAdj)
          
          # Health-adjusted life expectancy (Random draws)
          exMat <- t(apply(cbind(ex, (exUp - exLow) / 3.92), 1,
                           function(z) {
                             rtnorm(ndraw, mean = z[1], sd = z[2], lower = 0)   
                           }))
          
          # Re-define age and ax vectors
          nage <- length(ex)
          x <- x[1:nage]
          ax <- ax[1:nage]
          ax[nage] <- ex[nage]
          
          # Second-last ax: Decreasing linear trend
          if (ax[nage - 1] > ex[nage - 1]) {
            axNew <- (ax[nage - 2] + ax[nage]) / 2  
            if (axNew < ex[nage - 1]) ax[nage - 1] <- axNew
            rm(axNew)
          }
          
          # Adjust ax (needed for morbidity at higher ages in some cases)
          idAdj <- which(ex < ax)
          if (length(idAdj) > 0) ax[idAdj] <- ex[idAdj] * rho[idAdj - nage + length(rho) + 1]
          rm(idAdj)
          
          # Health-loss probabilities (Point estimates)
          qx <- CalcQx(x, ax, ex)
          
          # Health-loss probabilities (Random draws)
          qxMat <- apply(exMat, 2, function(z) {CalcQx(x, ax, z)})
          
          # Re-sample to ensure qx in [0, 1] from a TRUNCATED LOG-NORMAL DISTRIBUTION
          mu <- rowMeans(qxMat[-nage, ])
          qxMat[which(mu > 0.5), ] <- 1 - qxMat[which(mu > 0.5), ]
          qxMat[-nage, ] <- t(apply(qxMat[-nage, ], 1, function(z) {
            SampleLogNorm(n = ndraw, mean(z), quantile(z, c(0.025, 0.975)))
          }))
          qxMat[which(mu > 0.5), ] <- 1 - qxMat[which(mu > 0.5), ]
          rm(mu)
          
        }
        
        # Survival probabilities (Point estimates)
        lx <- c(1, cumprod(1 - qx[-nage]))
        
        # Survival probabilities (Random draws)
        lxMat <- apply(qxMat, 2, function(z) {cumprod(1 - z[-length(z)])})
        lxMat <- rbind(1, lxMat)
        
        # Distribution of deaths / health-loss (Point estimates)
        dx <- lx[-nage] - lx[-1]
        dx <- c(dx, lx[nage])
        
        # Distribution of deaths / health-loss (Random draws)
        dxMat <- lxMat[-nage, ] -  lxMat[-1, ]
        dxMat <- rbind(dxMat, lxMat[nage, ])
        
        
        #--------------------#
        # INEQUALITY MEASURE #
        #--------------------#
        
        # SD distribution of deaths / health-loss (Point estimates)
        sdDx <- CalcSD(x, ax, dx, ex, ages = Ages)
        
        # SD distribution of deaths / health-loss (Random draws)
        sdMat <- t(apply(array(data = c(dxMat, exMat), dim = c(dim(dxMat), 2)), 2,
                         function(z) {CalcSD(x, ax, z[, 1], z[, 2], Ages)}))
        
        # Uncertainty intervals
        sdLow <- apply(sdMat, 2, quantile, .1)
        sdUp <- apply(sdMat, 2, quantile, .9)
        
        
        #---------------#
        # STORE RESULTS #
        #---------------#
        
        # Output data frame
        results <- rbind(results,
                         data.frame(location_id = loc, location_name = locName,
                                    ISO3 = iso3, Type = type,
                                    SuperRegion = supReg, Region = reg, Year = year,
                                    Sex = sex, Age = Ages, Measure = measure,
                                    ex = ex[which(x %in% Ages)], 
                                    sdLow = sdLow, sd = sdDx, sdUp = sdUp))
        
        # Save life tables (avoiding this reduces computation time substantially)
        if (saveLT) {
          
          # Uncertainty intervals
          if (measure == 'Morbidity') {
            qxLow <- apply(qxMat, 1, quantile, 0.025)
            qxUp <- apply(qxMat, 1, quantile, 0.975)
          }
          lxLow <- apply(lxMat, 1, quantile, 0.025)
          lxUp <- apply(lxMat, 1, quantile, 0.975)
          dxLow <- apply(dxMat, 1, quantile, 0.025)
          dxUp <- apply(dxMat, 1, quantile, 0.975)
          
          # Life tables
          lifeTable <- data.frame(location_id = loc, location_name = locName,
                                  ISO3 = iso3, Type = type,
                                  SuperRegion = supReg, Region = reg, 
                                  Year = year, Sex = sex, Measure = measure,
                                  x = x, ax = ax,
                                  qx = qx, qxLow = qxLow, qxUp = qxUp,
                                  lx = lx, lxLow = lxLow, lxUp = lxUp,
                                  dx = dx, dxLow = dxLow, dxUp = dxUp, 
                                  ex = ex, exLow = exLow, exUp = exUp)
          
          # Store results
          if (measure == 'Mortality') ltMortality <- rbind(ltMortality, lifeTable)
          if (measure == 'Morbidity') ltMorbidity <- rbind(ltMorbidity, lifeTable)
          
          # Remove objects
          if (measure == 'Morbidity') rm(qxLow, qxUp)
          rm(lifeTable, lxLow, lxUp, dxLow, dxUp)
          
        }
        
        # HLI / LI ratios
        if (exists('sdDx') & exists('LI')) {
          
          # Uncertainty intervals
          ratioLow <- apply(sdMat / LImat, 2, quantile, .1)
          ratioUp <- apply(sdMat / LImat, 2, quantile, .9)
          
          # Output data frame
          ratios <- rbind(ratios,
                          data.frame(location_id = loc, location_name = locName,
                                     ISO3 = iso3, Type = type,
                                     SuperRegion = supReg, Region = reg, 
                                     Year = year, Sex = sex, Age = Ages,
                                     RatioLow = ratioLow, Ratio = sdDx / LI, 
                                     RatioUp = ratioUp))
          
          # Remove objects
          rm(LI, LImat, ratioLow, ratioUp, sdDx, sdMat)

        }
        
        # Remove objects
        if (measure == 'Mortality') rm(qxLow, qxUp)
        rm(qx, qxMat, lx, lxMat, dx, dxMat, ex, exLow, exUp, exMat, sdLow, sdUp)
        
      }
      
      # Remove objects
      rm(ax, nage, x)
      
    }
    
  }
  
  # Reduce size of input data
  datLT <- datLT[which(datLT$location_id != loc), ]
  datHALE <- datHALE[which(datHALE$location_id != loc), ]
  
  # Remove objects
  rm(iso3, loc, locName, measure, reg, sex, supReg, type, year)
  
}

End <- Sys.time()
print(End - Start)


# SAVE ESTIMATES
SAVE <- T
if (SAVE) {
  
  # LIFE EXPECTANCY AND INEQUALITY INDICATORS
  write.csv(results, row.names = F, 
            file = paste0('Results/EstimatesHLI-LI-', min(Years), '-', max(Years), '.csv'))
  
  # HLI / LI RATIOS
  write.csv(ratios, row.names = F, 
            file = paste0('Results/RatiosHLI-LI-', min(Years), '-', max(Years), '.csv'))
  
  # LIFE TABLES
  if (saveLT) {
    
    # Mortality
    write.csv(ltMortality, row.names = F, 
              file = paste0('Results/ltMortality-', min(Years), '-', max(Years), '.csv'))
    
    # Morbidity
    write.csv(ltMorbidity, row.names = F, 
              file = paste0('Results/ltMorbidity-', min(Years), '-', max(Years), '.csv'))
    
  }
  
}

