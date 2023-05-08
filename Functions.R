#############################################################################
#                                                                           #  
# HEALTHY LIFESPAN INEQUALITY: MORBIDITY COMPRESSION FROM                   #
#  A GLOBAL PERSPECTIVE                                                     #
#                                                                           #
# Inaki PERMANYER - Francisco VILLAVICENCIO - Sergi TRIAS-LLIMOS            #
#                                                                           #
# SET OF FUNCTIONS TO REPRODUCE THE RESULTS FROM THE PAPER                  #   
# Sourced by: CalcHLI.R                                                     #
#                                                                           #
#                                                             February 2023 #
#############################################################################


############################################################
### FUNCTIONS                                            ###  
############################################################

# CALCULATE ax VALUES FROM qx AND ex
CalcAx <- function(x, qx, ex) {
  
  ## x         Vector with ages (beginning of the age interval)
  ## qx        Age-specific probabilities of death
  ## ex        Age-specific remaining life expectancy
  
  # Length of the age intervals
  nx <- x[-1] - x[-length(x)]
  
  # Estimate ax (average person-years lived by those dying in the interval)
  ax <- (ex[-length(ex)] - (ex[-1] + nx)*(1 - qx[-length(qx)])) / qx[-length(qx)]
  
  # Last age group
  ax <- c(ax, ex[length(ex)])
  
  # Return ax values
  return(ax)
  
}


# CALCULATE DEATH PROBABILITIES FROM AGE-SPECIFIC LIFE EXPECTANCY DATA
CalcQx <- function(x, ax, ex) {
  
  ## x        Vector with ages (beginning of the age interval)
  ## ax       Average person-years lived by those dying between x and x+n
  ## ex       Age-specific remaining (healthy) life expectancy
  
  # Length of the age intervals
  nx <- x[-1] - x[-length(x)]
  
  # Unconditional age-specific probabilities of death
  qx <- 1 - (ex[-length(ex)] - ax[-length(ax)]) / (ex[-1] + nx - ax[-length(ax)])
  
  # Last age group
  qx <- c(qx, 1)
  
  # Output
  return(qx)
  
}


# FUNCTION TO SAMPLE FROM A (TRUNCATED) LOG-NORMAL DISTRIBUTION
SampleLogNorm <- function(n, mu, lim) {
  
  ## n        Number of random values
  ## mu       Mean value (reported point estimate)
  ## lim      Lower and upper bounds of the 95% UI
  
  # Estimate normal SD
  s <- (lim[2] - lim[1]) / 3.92
  
  # Log-normal mean
  muLog <- log(mu^2 / sqrt(s^2 + mu^2))
  
  # Log-normal SD
  sLog <- sqrt(log(1 + (s^2 / mu^2)))
  
  # Draw random values (TRUNCATED LOG-NORMAL DISTRIBUTION)
  randVal <- rlnormTrunc(n, meanlog = muLog, sdlog = sLog, max = 1)
  
  # Output
  return(randVal)
  
}


# CALCULATE SD of the life table DISTRIBUTION OF DEATHS
CalcSD <- function(x, ax, dx, ex, ages = c(0, 65)) {
  
  ## x        Vector with ages (beginning of the age interval)
  ## ax       Average person-years lived by those dying between x and x+n
  ## dx       Age-specific deaths / health-loss counts
  ## ex       Age-specific remaining (healthy) life expectancy
  ## ages     Vector with ages above which SD is calculated (OPTIONAL)
  
  # Output vector
  sdOut <- c()
  
  # Calculate SD for the desired ages
  for (i in ages) {
    y <- (x[x >= i] + ax[x >= i] - ex[x == i] - i)^2
    sdOut <- c(sdOut, (crossprod(y, dx[x >= i]) / sum(dx[x >= i]))^0.5)
  }
  
  # Output
  return(sdOut)
  
}


# FUNCTION TO CALCULATE (HEALTHY) LIFESPAN INEQUALITY MEASURES
CalcHLI <- function(id, lt, hale, countries, years, ages, sexGr, UI = NULL, n = NULL) {
  
  ## id           Country/regions ID (from GBD)
  ## lt           Life table data (qx and ex) from GBD
  ## hale         Health-adjusted life expectancy data from GBD
  ## countries    List of countries and their regional classification
  ## years        Years of interest
  ## ages         Ages for which inequality measures are calculated
  ## sexGr        Sex groups
  ## UI           Uncertainty interval to report
  ## n            Number of random draws to assess uncertainty
  
  
  # Bounds of the uncertainty interval
  if (!is.null(UI)) {
    if (UI <= 0 | UI >= 1) {
      UI <- c(.1, .9)
    } else UI <- .5 + c(-UI/2, UI/2) 
  }
  
  # Data frame to store LI and HLI estimates
  results <- data.frame()
  
  # Data frame to store HLI / LI ratios
  ratios <- data.frame()
  
  # LOOP ACROSS LOCATIONS
  for (loc in id) {
    
    # IDENTIFYING VARIABLES
    if (loc %in% c(1, countries$super_region_id)) {
      if (loc == 1) {
        locName <- 'Global'
        type <- 'Global'
        supReg <- NA
      } else {
        locName <- unique(countries$super_region_name[countries$super_region_id == loc])
        type <- 'Super-region'
        supReg <- locName
      }
      iso3 <- NA
      reg <- NA
    }
    if (loc %in% countries$region_id) {
      locName <- unique(countries$region_name[countries$region_id == loc])
      type <- 'Region'
      # North Africa and South Asia
      if (loc %in% c(137, 158)) type <- 'Super-region/Region'
      iso3 <- NA
      supReg <- unique(countries$super_region_name[countries$region_name == locName])
      reg <- locName
    }
    if (loc %in% countries$location_id) {
      locName <- unique(countries$location_name[countries$location_id == loc])
      type <- 'Country/Territory'
      iso3 <- countries$ISO3[countries$location_id == loc]
      supReg <- countries$super_region_name[countries$location_id == loc]
      reg <- countries$region_name[countries$location_id == loc]
    }
    
    # LOOP ACROSS YEARS
    for (year in years) {
      
      # LOOP ACROSS SEX GROUPS
      for (sex in sexGr) {
        
        # Vector of ages
        x <- sort(unique(lt$age))
        nage <- length(x)
        
        # Input data: Death probabilities
        qx <- lt$val[lt$location_id == loc & lt$year == year &
                       lt$sex == sex & lt$measure_name == 'Probability of death']
        qxLow <- lt$lower[lt$location_id == loc & lt$year == year &
                            lt$sex == sex & lt$measure_name == 'Probability of death']
        qxUp <- lt$upper[lt$location_id == loc & lt$year == year &
                           lt$sex == sex & lt$measure_name == 'Probability of death']
        
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
        if (!is.null(UI)) {
          qxMat <- t(apply(cbind(qx[-nage], (qxUp[-nage] - qxLow[-nage]) / 3.92), 1,
                           function(z) {
                             rtnorm(n, mean = z[1], sd = z[2], lower = 0, upper = 1)   
                           }))
          qxMat <- rbind(qxMat, 1)  
        }
        
        # Input data: Life expectancy
        ex <- lt$val[lt$location_id == loc & lt$year == year &
                       lt$sex == sex & lt$measure_name == 'Life expectancy']
        exLow <- lt$lower[lt$location_id == loc & lt$year == year &
                            lt$sex == sex & lt$measure_name == 'Life expectancy']
        exUp <- lt$upper[lt$location_id == loc & lt$year == year &
                           lt$sex == sex & lt$measure_name == 'Life expectancy']
        
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
        if (!is.null(UI)) {
          exMat <- t(apply(cbind(ex, (exUp - exLow) / 3.92), 1,
                           function(z) {
                             rtnorm(n, mean = z[1], sd = z[2], lower = 0)   
                           }))  
        }
        
        # Average person-years lived by those dying between x and x+n
        ax <- CalcAx(x, qx, ex)
        
        
        # Two measures: MORTALITY and MORBIDITY
        for (measure in c('Mortality', 'Morbidity')) {
          
          #-----------------------------------#
          # LIFE TABLE DISTRIBUTION OF DEATHS #
          #-----------------------------------#
          
          # Morbidity
          if (measure == 'Morbidity') {
            
            # Store existing LI estimates to calculate HLI / LI ratio
            LI <- SDdx
            if (!is.null(UI)) LImat <- SDmat
            
            # Input data: Health-adjusted life expectancy
            ex <- hale[hale$location_id == loc & hale$year == year &
                         hale$sex == sex, 'val']
            exLow <- hale[hale$location_id == loc & hale$year == year &
                            hale$sex == sex, 'lower']
            exUp <- hale[hale$location_id == loc & hale$year == year &
                           hale$sex == sex, 'upper']
            
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
            if (!is.null(UI)) {
              exMat <- t(apply(cbind(ex, (exUp - exLow) / 3.92), 1,
                               function(z) {
                                 rtnorm(n, mean = z[1], sd = z[2], lower = 0)   
                               }))  
            }
            
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
            if (!is.null(UI)) {
              qxMat <- apply(exMat, 2, function(z) {CalcQx(x, ax, z)})
              # Re-sample to ensure qx in [0, 1] from a TRUNCATED LOG-NORMAL DISTRIBUTION
              mu <- rowMeans(qxMat[-nage, ])
              qxMat[which(mu > 0.5), ] <- 1 - qxMat[which(mu > 0.5), ]
              qxMat[-nage, ] <- t(apply(qxMat[-nage, ], 1, function(z) {
                SampleLogNorm(n = n, mean(z), quantile(z, c(0.025, 0.975)))
              }))
              qxMat[which(mu > 0.5), ] <- 1 - qxMat[which(mu > 0.5), ]
              rm(mu)
            }
            
          }
          
          # Survival probabilities (Point estimates)
          lx <- c(1, cumprod(1 - qx[-nage]))
          
          # Survival probabilities (Random draws)
          if (!is.null(UI)) {
            lxMat <- apply(qxMat, 2, function(z) {cumprod(1 - z[-length(z)])})
            lxMat <- rbind(1, lxMat)  
          }
          
          # Distribution of deaths / health-loss (Point estimates)
          dx <- lx[-nage] - lx[-1]
          dx <- c(dx, lx[nage])
          
          # Distribution of deaths / health-loss (Random draws)
          if (!is.null(UI)) {
            dxMat <- lxMat[-nage, ] -  lxMat[-1, ]
            dxMat <- rbind(dxMat, lxMat[nage, ])  
          }
          
          
          #--------------------#
          # INEQUALITY MEASURE #
          #--------------------#
          
          # SD distribution of deaths / health-loss (Point estimates)
          SDdx <- CalcSD(x, ax, dx, ex, ages)
          
          # SD distribution of deaths / health-loss (Random draws)
          if (!is.null(UI)) {
            SDmat <- t(apply(array(data = c(dxMat, exMat), dim = c(dim(dxMat), 2)), 2,
                             function(z) {CalcSD(x, ax, z[, 1], z[, 2], ages)}))
            # Uncertainty intervals
            SDlow <- apply(SDmat, 2, quantile, UI[1])
            SDup <- apply(SDmat, 2, quantile, UI[2])  
          }
          
          
          #---------------#
          # STORE RESULTS #
          #---------------#
          
          # Store inequality measure
          dat <- data.frame(location_id = loc, location_name = locName,
                            ISO3 = iso3, Type = type,
                            SuperRegion = supReg, Region = reg, Year = year,
                            Sex = sex, Age = ages, Measure = measure,
                            ex = ex[which(x %in% ages)], SD = SDdx)
          
          # Uncertainty
          if (!is.null(UI)) {
            dat[, paste0('SDlow', round(100 * UI[1]))] <- SDlow
            dat[, paste0('SDup', round(100 * UI[2]))] <- SDup
          }
          
          # Output data frame
          results <- rbind(results, dat)
          rm(dat)
          
          # HLI / LI ratios
          if (exists('SDdx') & exists('LI')) {
            
            # Store ratios
            dat <- data.frame(location_id = loc, location_name = locName,
                              ISO3 = iso3, Type = type,
                              SuperRegion = supReg, Region = reg, 
                              Year = year, Sex = sex, Age = ages,
                              Ratio = SDdx / LI)
            
            # Uncertainty
            if (!is.null(UI)) {
              ratioLow <- apply(SDmat / LImat, 2, quantile, UI[1])
              ratioUp <- apply(SDmat / LImat, 2, quantile, UI[2])
              dat[, paste0('RatioLow', round(100 * UI[1]))] <- ratioLow
              dat[, paste0('RatioUp', round(100 * UI[2]))] <- ratioUp
            }
            
            # Output data frame
            ratios <- rbind(ratios, dat)
            rm(dat)
            
            # Remove objects
            rm(LI, SDdx)
            if (!is.null(UI)) rm(LImat, ratioLow, ratioUp, SDmat)
            
          }
        
          # Remove objects
          if (measure == 'Mortality') rm(qxLow, qxUp)
          if (!is.null(UI)) rm(dxMat, exMat, lxMat, qxMat, SDlow, SDup)
          rm(dx, ex, exLow, exUp, lx, qx)
      
        }
        
        # Remove objects
        rm(ax, nage, x)
        
      }
      
    }
    
    # Reduce size of input data
    lt <- lt[which(lt$location_id != loc), ]
    hale <- hale[which(hale$location_id != loc), ]
    
    # Remove objects
    rm(iso3, loc, locName, measure, reg, sex, supReg, type, year)
    
  }
  
  # Output: Inequality measures and ratios
  return(list(results = results, ratios = ratios))
  
}

