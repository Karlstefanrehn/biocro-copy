## This is the Penman-Monteith hourly calculation
## This is based on the document PMhrDoc found here
## http://biomet.ucdavis.edu/evapotranspiration.html


penmon <- function(hr=12, doy=200, lat=40, long=90, stat.long=90, elev=100){

## Defining inputs
  ## hr = hour of the day (0-24)
  ## doy = day of the year (1-366)
  ## lat = latitude (degrees)
  ## long = longitude (degrees)
  
  ## Defining constants
  G.sc <- 0.082 ## solar constant MJ m-2 min-1
  sigma <- 2.04e-10 ## Steffan-Boltzman MJ m-2 hr-1 K-4

  ## Step 1
  ## preliminary calculations
  phi <- (lat * pi)/180 ## Convert latitude from degrees to radians
  J <- doy ## This is to stay consistent with the notation
  d.r <- 1 + 0.033 * cos(((2 * pi)/365) * J) ## correction for eccentricity of Earth's orbit around the sun
  delta <- 0.409 * sin(((2 * pi)/365) * J - 1.39) ## Declination of the sun above the celestial equator in radians
  b <- (2 * pi * (J - 81)) / 364
  S.c <- 0.1645 * sin(2 * b) - 0.1255 * cos(b) - 0.025 * sin(b)
  t <- hr ## local standard time
  omega <- (pi/12) * ((t - 0.5) + ((long - stat.long)/15) - 12 + S.c) ## hour angle in radians
  omega.1 <- omega - (0.5) * (pi/12) ## hour angle 1/2 hour before omega in radians
  omega.2 <- omega + (0.5) * (pi/12) ## hour angle 1/2 hour after omega in radians
  sin.theta <- (omega.2 - omega.1) * sin(phi) * sin(delta) + cos(phi) * cos(delta) * (sin(omega.2) - sin(omega.1))
  R.a <- 12/pi * (60 * G.sc) * d.r * sin.theta
  beta0 <- sin(phi) * sin(delta) + cos(phi) * cos(delta) * cos(omega) 
  beta <- 180/pi * (1/sin(beta0))

  ## Step 2
  ## R.so is the clear sky total global radiation at the Earth's surface
  R.so <- R.a * (0.75 + 2.0e-5 * elev) ## MJ m-2 hr-1

  R.so

}

## Test the function
## hrs <- 0:24
## ans <- penmon(hr=hrs, lat=0)

## convert MJ m-2 hr-1 to mu mol m-2 s-1
## ans2 <- ans * (1/(60*60)) * 2.07 * 1e6
## xyplot(ans2 ~ hrs)
