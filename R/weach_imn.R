##
##  BioCro/R/weach_imn.R by Fernando Ezequiel Miguez  Copyright (C) 2011-2015
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

weach_imn <- function(data, ts=1,
                      rad.units=c("kilocalories"),
                      temp.units=c("Fahrenheit","Celsius"),
                      rh.units=c("percent","fraction"),
                      ws.units=c("mph","mps"),
                      pp.units=c("in","mm"),...){
  
  if((ts<1)||(24%%ts != 0))
    stop("ts should be a divisor of 24 (e.g. 1,2,3,4,6,etc.)")

  if(ts != 1) stop("This function only accepts hourly data at the moment")
  
  if(dim(data)[2] != 9)
    stop("X should have 9 columns")

  ## Need to get rid of the last row
  if(nrow(data) == 8761) data <- data[-8761,]

  MPHTOMPERSEC <- 0.447222222222222

  rad.units <- match.arg(rad.units)
  temp.units <- match.arg(temp.units)
  rh.units <- match.arg(rh.units)
  ws.units <- match.arg(ws.units)
  pp.units <- match.arg(pp.units)

  ## Collect the variables

  site <- data[,1]
  loc  <- data[,2]
  date <- data[,3]
  hour <- data[,4]
  temp <- data[,5]
  solarR <- data[,6]
  precip <- data[,7]
  RH <- data[,8]
  windS <- data[,9]

  ## Transform them in to needed input

  date <- as.Date(date)
  year <- as.numeric(format(date, "%Y"))
  doy <- as.numeric(format(date, "%j"))
  hour <- as.numeric(as.vector(sapply(as.character(hour),
                                      function(x) strsplit(x, ":")[[1]][1])))

##  temp <- ifelse(temp == -99, NA, temp)
  
  if(temp.units == "Fahrenheit"){
    temp <- (temp - 32) * (5/9)
  }

  ## the solar radiation is given in kilo calories per hour per meter squared
  ## To convert from kilocalories to joules
  ## 1 kilocalorie = 4184 joules
  ## To convert to Mega Joules
  solarR <- ifelse(solarR == -99, NA, solarR)
##  solarR <- solarR * 0.86
  solarR0 <- (solarR * 4184) * 1e-6 ## This is in MJ/hr/m2
  solarR <- solarR0 * 2.07 * 1e6 / 3600 ## Look for comments in the
                                        ## source code of the weach function
                                        ## For details
  if(pp.units == "in"){
    precip <- precip*2.54*10
  }

  if(rh.units == "percent"){
    RH <- RH / 100
  }
  
  if(ws.units == "mph"){
    windS <- windS * MPHTOMPERSEC
  }

  res <- data.frame(year = year, doy = doy, hour = hour, solarR = solarR,
                    temp = temp, RH = RH, windS = windS, precip = precip)

  res

}


weach_imn2 <- function(data,
                       rad.units = c("Watt/m2"),
                       temp.units=c("Fahrenheit","Celsius"),
                      rh.units=c("percent","fraction"),
                      ws.units=c("mph","mps"),
                      pp.units=c("in","mm"),na.chr=-99){

    ## This function is meant to be used with data extracted from  
    ## http://mesonet.agron.iastate.edu/agclimate/hist/hourly.php
 
  if(dim(data)[2] != 16)
    warning("data should have 16 columns")

##  if(dim(data)[1] != 8760)
##      warning("data should have 8760 rows \n missing data will be included")

  clnms <- c('station','valid','tmpf','relh','solar','precip','speed','drct','et',
             'soil04t','soil12t','soil24t','soil50t','soil12vwc','soil24vwc','soil50vwc')

  dtnms <- names(data)

  if(!all(!is.na(match(clnms,dtnms)))){
    stop(cat('columns should be ',clnms,'\n')) 
  }
  ## Conversion factors needed
  MPHTOMPERSEC <- 0.447222222222222
  
  ## I won't use the station name at this point
  ## Process the data for time
  date0 <- as.vector(sapply(as.character(data$valid),
                                      function(x) strsplit(x, " ")[[1]][1]))
  date0 <- as.Date(date0)
  year <- as.numeric(as.character(format(date0,"%Y")))
  doy <- as.numeric(as.character(format(date0,"%j")))

  hour0 <- as.vector(sapply(as.character(data$valid),
                                      function(x) strsplit(x, " ")[[1]][2]))
  
  hour <- as.numeric(as.vector(sapply(as.character(hour0),
                                      function(x) strsplit(x, ":")[[1]][1])))

  rad.units <- match.arg(rad.units)
  temp.units <- match.arg(temp.units)
  rh.units <- match.arg(rh.units)
  ws.units <- match.arg(ws.units)
  pp.units <- match.arg(pp.units)

  ## Collect the variables
    temp <- data[,3]
    RH <- data[,4]
    solarR <- data[,5]
    precip <- data[,6]
    windS <- data[,7]
    drct <- data[,8]
    et <- data[,9]
    soil04t <- data[,10]
    soil12t <- data[,11]
    soil24t <- data[,12]
    soil50t <- data[,13]
    soil12vwc <- data[,14]
    soil24vwc <- data[,15]
    soil50vwc <- data[,16]
  
  ## Transform them in to needed input

  temp <- ifelse(temp == na.chr, NA, temp)
    soil04t <- ifelse(soil04t == na.chr, NA, soil04t)
    soil12t <- ifelse(soil12t == na.chr, NA, soil12t)
    soil24t <- ifelse(soil24t == na.chr, NA, soil24t)
    soil50t <- ifelse(soil50t == na.chr, NA, soil50t)
  
  if(temp.units == "Fahrenheit"){
    temp <- (temp - 32) * (5/9)
    soil04t <- (soil04t - 32) * (5/9)
    soil12t <- (soil12t - 32) * (5/9)
    soil24t <- (soil24t - 32) * (5/9)
    soil50t <- (soil50t - 32) * (5/9)
  }

    soil12vwc <- ifelse(soil12vwc == na.chr, NA, soil12vwc)
    soil24vwc <- ifelse(soil24vwc == na.chr, NA, soil24vwc)
    soil50vwc <- ifelse(soil50vwc == na.chr, NA, soil50vwc)
    
  ## The solar radiation is given in W/m2, in the past it was in 
  ## kilo calories per hour per meter squared
  ## Convert from W/m2 to kcal h / m2 multiply by 0.86
  ## http://pveducation.org/pvcdrom/appendices/units-and-conversions
  ## To convert from kilocalories to joules
  ## 1 kilocalorie = 4184 joules
  ## To convert to Mega Joules
  solarR <- ifelse(solarR == na.chr, NA, solarR)
  
  solarR <- solarR * 0.86
  solarR0 <- (solarR * 4184) * 1e-6 ## This is in MJ/hr/m2
  solarR <- solarR0 * 2.07 * 1e6 / 3600 ## Look for comments in the
                                        ## source code of the weach function
                                        ## For details
##  cat("solarR length",length(solarR),"\n")
  
  if(pp.units == "in"){
    precip <- precip*2.54*10
##    cat("precip length",length(precip),"\n")
  }

  if(rh.units == "percent"){
    RH <- RH / 100
##    cat("RH length",length(RH),"\n")
  }
  
  if(ws.units == "mph"){
    windS <- windS * MPHTOMPERSEC
##    cat("windS length",length(windS),"\n")
  }

  res <- data.frame(year = year, doy = doy,
                    hour = hour, solarR = solarR,
                    temp = temp, RH = RH, windS = windS, precip = precip,
                    soil04t = soil04t, soil12t = soil12t, soil24t = soil24t, soil50t = soil50t,
                    soil12vwc = soil12vwc, soil24vwc = soil24vwc, soil50vwc = soil50vwc, et = et)
  
  ## The code below
  ## Filling in with missing data
  ## res <- data.frame(year = unique(year), doy = rep(1:365, each=24),
  ##                   hour = rep(0:23,365), solarR = NA,
  ##                   temp = NA, RH = NA, windS = NA, precip = NA)

  ## tmp <- cbind(solarR,temp,RH,windS,precip)

  ## Trying an approach where I include one row at a time
  ## k <- 1
  ## for(i in 1:length(unique(doy))){
  ##     for(j in 0:23){
  ##         doy1 <- which(res[,2] == unique(doy)[i])
  ##         doyhr1 <- doy1[1] + hour[j]
  ##         res[doyhr1,4:8] <- 
  ##     }
  ## }
  ## ## Which is the first row in the observed data
  ## doy1 <- which(res[,2] == doy[1])
  ## doyhr1 <- doy1[1] + hour[1]
  ## ## Which is the last row in the observed data
  ## doyn <- which(res[,2] == doy[length(doy)])
  ## doyhrn <- doyn[1] + hour[length(hour)]


  
  ## print(c(doyhr1,doyhrn))
  ## print(dim(tmp))
  ## res[doyhr1:doyhrn,4:8] <- tmp
  
  ## res <- as.data.frame(res)

  res

}

