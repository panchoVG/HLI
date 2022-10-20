#############################################################################
#                                                                           #  
# HEALTHY LIFESPAN INEQUALITY: MORBIDITY COMPRESSION FROM                   #
#  A GLOBAL PERSPECTIVE                                                     #
#                                                                           #
# Inaki PERMANYER - Francisco VILLAVICENCIO - Sergi TRIAS-LLIMOS            #
#                                                                           #
# SET OF FUNCTIONS TO REPRODUCE RESULTS AND FIGURES FROM THE PAPER          #   
# Sourced by: CalcHLI.R                                                     #
#                                                                           #
#                                                                 July 2022 #
#############################################################################


############################################################
### MAIN FUNCTIONS                                       ###  
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
