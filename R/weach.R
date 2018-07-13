##
##  BioCro/R/weach.R by Fernando Ezequiel Miguez  Copyright (C) 2007-2014
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 or 3 of the License
##  (at your option).
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  A copy of the GNU General Public License is available at
##  http://www.r-project.org/Licenses/
##
##

weach <- function(X,lat=40,ts=1,solar.units=c("MJ/m2"), temp.units=c("Fahrenheit","Celsius"),
                  rh.units=c("percent","fraction"),ws.units=c("mph","mps"),
                  pp.units=c("in","mm"),seed=1234,...){

    set.seed(seed)
    
    if(missing(lat))
        stop("latitude is missing")
  
    if((ts<1)||(24%%ts != 0))
        stop("ts should be a divisor of 24 (e.g. 1,2,3,4,6,etc.)")
  
    if(dim(X)[2] != 11)
        stop("X should have 11 columns")

    if(nrow(X) == 366 || is.leap(X[1,1])){
        warning("leap year")
        year.days <- 366
    }else{
        year.days <- 365
    }
  
  MPHTOMPERSEC <- 0.447222222222222

  solar.units <- match.arg(solar.units)
  temp.units <- match.arg(temp.units)
  rh.units <- match.arg(rh.units)
  ws.units <- match.arg(ws.units)
  pp.units <- match.arg(pp.units)

  if(solar.units != "MJ/m2") stop("MJ/m2 is the only option at the moment")
  
  year <- X[,1]
  DOYm <- X[,2]
  solar <- X[,3]
  maxTemp <- X[,4]
  minTemp <- X[,5]
  avgTemp <- X[,6]
  maxRH <- X[,7]
  minRH <- X[,8]
  avgRH <- X[,9]
  WindSpeed <- X[,10]
  precip <- X[,11]

  tint <- 24/ts
  tseq <- seq(0,23,ts)

  ## Solar radiation
  solarR <- (0.12*solar) * 2.07 * 1e6 / 3600
  ## This last line needs some explanation
  ## There is no easy straight way to convert MJ/m2 to mu mol photons / m2 / s (PAR)
  ## The above conversion is based on the following reasoning
  ## 0.12 is about how much of the total radiation is expected to ocurr during the hour of maximum insolation (it is a guesstimate)
  ## 2.07 is a coefficient which converts from MJ to mol photons (it is approximate and it is taken from ...
  ## Campbell and Norman (1998). Introduction to Environmental Biophysics. pg 151 'the energy content of solar radiation in the PAR waveband is 2.35 x 10^5 J/mol'
  ## See also the chapter radiation basics (10)
  ## Here the input is the total solar radiation so to obtain in the PAR spectrum need to multiply by 0.486
  ## This last value 0.486 is based on the approximation that PAR is 0.45-0.50 of the total radiation
  ## This means that 1e6 / (2.35e6) * 0.486 = 2.07
  ## 1e6 converts from mol to mu mol
  ## 1/3600 divides the values in hours to seconds
  solarR <- rep(solarR , each = tint)

  ltseq <- length(tseq)
  resC2 <- numeric(ltseq*year.days)
  for(i in 1:year.days)
    {
      res <- lightME(DOY = i , t.d = tseq, lat = lat, ...)
      Itot <- res$I.dir + res$I.diff
      indx <- 1:ltseq + (i-1)*ltseq 
      resC2[indx] <- (Itot-min(Itot))/max(Itot)
    }

  SolarR <- solarR * resC2

  ## Temperature
  if(temp.units == "Fahrenheit"){
    minTemp <- (minTemp - 32)*(5/9)
    minTemp <- rep(minTemp , each = tint)
    maxTemp <- (maxTemp - 32)*(5/9)
    maxTemp <- rep(maxTemp , each = tint)
    rangeTemp <- maxTemp - minTemp
  }else{
    minTemp <- rep(minTemp , each = tint)
    maxTemp <- rep(maxTemp , each = tint)
    rangeTemp <- maxTemp - minTemp
  }

  xx <- rep(tseq,year.days)
  temp1 <- sin(2 * pi * (xx - 10)/tint)
  temp1 <- (temp1 + 1)/2
  Temp <- minTemp + temp1 * rangeTemp

  ## Relative humidity
  minRH <- rep(minRH,each=tint)
  maxRH <- rep(maxRH,each=tint)

  temp2 <- cos(2 * pi * (xx - 15)/tint)
  temp2 <- (-temp2 + 1)/2
  if(rh.units == "percent"){
    RH <- (minRH + temp2 * (maxRH - minRH))/100
  }else{
    RH <- (minRH + temp2 * (maxRH - minRH))
  }

  RH <- ifelse(RH > 0.99999, 0.9999, RH) 
  
  ## Wind Speed
  temp3 <- temp1 + 0.5
  if(ws.units == "mph"){
    WS <- temp3 * rep(WindSpeed,each=tint) * MPHTOMPERSEC
  }else{
    WS <- temp3 * rep(WindSpeed,each=tint)
  }

  ## Precipitation Trying a binomial distribution to have rain
  ## unequally spaced during the day
    precip <- rep(precip, each = tint)
    resC3 <- numeric(ltseq*year.days)

    for(i in 1:year.days){
        x <- rbinom(tint, 1, 0.35)
        xp <- x/sum(x)
        indx <- 1:ltseq + (i-1)*ltseq
        resC3[indx] <- xp
    }
      
    if(pp.units == "in"){
        precip <- precip * resC3 * 2.54*10
    }else{
        precip <- precip * resC3 
    }
  
  hour <- rep(tseq,year.days)
  DOY <- 1:year.days
  doy <- rep(DOY,each=tint)


  ans <- as.data.frame(cbind(year,doy,hour,SolarR,Temp,RH,WS,precip))
  ans
}


## This function is not really needed, but here it is in case I need
## it in the future
is.leap <- function(year){
  if(year %% 4 == 0){
    if(year %% 100 == 0){
       if(year %% 400 == 0){
         ans <- TRUE
         }else{
           ans <- FALSE
         }
     }else{
       ans <- TRUE
     }
  }else{
    ans <- FALSE
  }
  ans
}
