##
##  BioCro/R/MaizeGro.R by Fernando Ezequiel Miguez  Copyright (C) 2012-2015
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

MaizeGro <- function(WetDat, plant.day=NULL,
                     emerge.day=NULL,
                     harvest.day=NULL,
                     plant.density = 7, 
                     timestep=1,
                     lat=40,
                     canopyControl=list(),
                     MaizeSeneControl=list(),
                     photoControl=list(),
                     MaizePhenoControl=list(),
                     MaizeCAllocControl=list(),
                     laiControl=list(),
                     soilControl=list(),
                     MaizeNitroControl=list(),
                     centuryControl=list())
  {

    ## Sanity check for weather data
    if(any(is.na(WetDat))) stop("No missing data allowed")
    nrowdata <- nrow(WetDat)
    if(nrowdata < 8760){
      warning("The weather data does not have \n hourly data for a whole year")
      warning("NA will be introduced and the last day \n with data will be the harvest day")
      if(missing(harvest.day)) harvest.day <- max(WetDat[,2])
      doys <- rep(1:365, each = 24)
      hours <- rep(0:23, 365)
      dummy.wetdat <- data.frame(year = unique(WetDat[,1]), doy = doys, hour = hours,
                                  solarR = 0, Temp = 0, RH = 0, WS = 0, precip = 0)
      dummy.wetdat[1:nrowdata,] <- WetDat
      WetDat <- dummy.wetdat
    }

    
    ## How should I guess the first day of growth?
    if(missing(plant.day))
      stop("planting date should be supplied")

    if(missing(emerge.day)) emerge.day <- -1

    ## If harvest day is missing it will be the day of the first frost
    ## in fall. The latest day allowed is 270 or Sept 27 th
    ## Completely arbitrary for now
    if(missing(harvest.day)){
        half <- as.integer(dim(WetDat)[1]/2)
        WetDat1 <- WetDat[half:dim(WetDat)[1],c(2,5)]
        if(min(WetDat1[,2]) > 0){
          harvest.day <- 270
        }else{
          WetDat1s <- WetDat1[which(WetDat1[,2]<0),]
          harvest.day <- min(WetDat1s[,1])
          if(harvest.day > 270) harvest.day <- 270
        }
      }
    
    if((plant.day<0) || (plant.day>365) || (harvest.day<0) || (harvest.day>365))
      stop("plant.day and harvest.day should be between 0 and 365")

    if(plant.day > harvest.day)
      stop("plant.day should be smaller than harvest.day")

    if( (timestep<1) || (24%%timestep != 0))
      stop("timestep should be a divisor of 24 (e.g. 1,2,3,4,6,etc.)")

    ## Getting the Parameters
    ## Canopy parameters will stay the same
    canopyP <- canopyParms()
    canopyP[names(canopyControl)] <- canopyControl
    canopyP <- c(unlist(canopyP))

    ## Getting the senescence parameters
    seneP <- MaizeSeneParms()
    seneP[names(MaizeSeneControl)] <- MaizeSeneControl
    seneP <- c(unlist(seneP))
    
    ## Soil and nitro parameters will stay the same
    soilP <- soilParms()
    soilP[names(soilControl)] <- soilControl

    maizenitroP <- MaizeNitroParms()
    maizenitroP[names(MaizeNitroControl)] <- MaizeNitroControl
    maizenitroP <- c(unlist(maizenitroP))

    ## Need a specialized phenology control
    maizephenoP <- MaizePhenoParms()
    maizephenoP[names(MaizePhenoControl)] <- MaizePhenoControl
    maizephenoP <- c(unlist(maizephenoP))
    
    ## Need special parameters for photosynthesis
    maizephotoP <- MaizePhotoParms()
    maizephotoP[names(photoControl)] <- photoControl
    maizephotoP <- c(unlist(maizephotoP))

    ## Need special parameters for photosynthesis
    laiP <- laiParms()
    laiP[names(laiControl)] <- laiControl
    laiP <- c(unlist(laiP))

    ## C allocation is distinct from phenology
    mCallocP <- MaizeCAllocParms()
    mCallocP[names(MaizeCAllocControl)] <- MaizeCAllocControl
    mCallocP <- c(unlist(mCallocP))

    ## Century model parameters
    centuryP <- centuryParms()
    centuryP[names(centuryControl)] <- centuryControl

    centCoefs <- as.vector(unlist(centuryP)[1:24])

    if(centuryP$timestep == "year"){
        stop("Not developed yet")
        centTimestep <- 365 ## This is really the growing season
    }
    if(centuryP$timestep == "week") centTimestep <- 7
    if(centuryP$timestep == "day") centTimestep <- 1
    if(centuryP$timestep == "hour") centTimestep <- 1/24
    
    
    res <- .Call("maizeGro",
                 as.integer(WetDat[,2]), ## Day of the year                         1
                 as.integer(WetDat[,3]), ## Hour of the day                         2
                 as.double(WetDat[,4]), ## Solar radiation                          3
                 as.double(WetDat[,5]), ## Temperature                              4
                 as.double(WetDat[,6]), ## Relative Humidity                        5
                 as.double(WetDat[,7]), ## Wind Speed                               6
                 as.double(WetDat[,8]), ## Precipitation                            7
                 as.integer(c(plant.day,
                              emerge.day, harvest.day)), ## Dates                   8
                 as.double(maizephenoP), ## Phenology parameters                    9
                 as.double(maizephotoP), ## Photosynthetic paramters               10
                 as.double(canopyP),     ## Canopy parameters                      11
                 as.double(maizenitroP), ## Nitro parameters                       12
                 as.double(laiP),        ## LAI parameters                         13
                 as.double(mCallocP),    ## C allocation coefficients              14
                 as.double(lat),         ## Latitude                               15
                 as.integer(timestep),   ## Time step                              16
                 as.double(plant.density),## plant density                         17
                 as.double(soilP[c(1:5,11:13,15:18)]),
                                 ## Soil parameters                                18
                 as.integer(soilP[c(7,8,10,14)]), ## Soil type and number of layers      19
                 as.double(soilP$soilDepths), ## Soil depths                       20
                 as.double(soilP$iWatCont), ## initial water status                 21
                 as.double(soilP$respcoef), ## respiration coefs                    22
                 as.double(soilP$acoef), ## soilT coef                              23
                 as.double(seneP), ## senescence parameters                         24
                 as.double(centCoefs), ## Century coefficients                      25
                 as.double(centTimestep), ## Century timestep                      26
                 as.double(centuryP$Ks) ## Century rate parameters                  27
                 )

    ## Transform the phenology

    res$PhenoStage <- phenoConv(res$PhenoStage)
    res$LAImat <- t(res$LAImat)

    res$cwsMat <- t(res$cwsMat)
    colnames(res$cwsMat) <- soilP$soilDepths[-1]
    res$rdMat <- t(res$rdMat)
    colnames(res$rdMat) <- soilP$soilDepths[-1]
    res$psimMat <- t(res$psimMat)
    colnames(res$psimMat) <- soilP$soilDepths[-1]
    
    structure(res, class = "MaizeGro")


  }


MaizePhenoParms <- function(base.temp = 10, max.leaves = 20, plant.emerg = 100,
                            phyllochron1 = 46.7, phyllochron2 = 31.1, R1 = 747,
                            R2 = 858, R3 = 969, R4 = 1080, R5 = 1136, R6 = 1452){

  list(base.temp = base.temp, max.leaves = max.leaves, plant.emerg = plant.emerg,
       phyllochron1 = phyllochron1, phyllochron2 = phyllochron2,
       R1 = R1, R2 = R2, R3 = R3, R4 = R4, R5 = R5, R6 = R6)


}

MaizeSeneParms <- function(senStem=3000,senLeaf=3500,senRoot=4000){

  list(senStem=senStem,senLeaf=senLeaf,senRoot=senRoot)

}

MaizePhotoParms <- function(vmax = 56, alpha = 0.05, kparm = 0.7,
                            theta = 0.83, beta = 0.93, Rd = 1.3,
                            Catm = 380, b0 = 0.01, b1 = 3, ws = c("gs", "vmax")){

  ws <- match.arg(ws)
  if (ws == "gs") 
    ws <- 1
  else ws <- 0

  Rd <- rep(Rd, length.out=2) ## This is a vector of length 2 first for Vegetative and second for Reproductive
  vmax <- rep(vmax, length.out=2) ## This is a vector of length 2 first for Vegetative and second for Reproductive
  
  list(vmax = vmax, alpha = alpha, kparm = kparm, theta = theta, 
        beta = beta, Rd = Rd, Catm = Catm, b0 = b0, b1 = b1, 
        ws = ws)

}

MaizeNitroParms <- function(iLeafN=5, kLN=0.25, Vmax.b1=0.569, alpha.b1=0,
                       kpLN=0.2, lnb0 = -5, lnb1 = 18, lnFun=c("none","linear")){

  lnFun <- match.arg(lnFun)
  if(lnFun == "none"){
    lnFun <- 0
  }else{
    lnFun <- 1
  }
  
  list(iLeafN=iLeafN, kLN=abs(kLN), Vmax.b1=Vmax.b1, alpha.b1=alpha.b1, kpLN=kpLN,
       lnb0 = lnb0, lnb1 = lnb1, lnFun = lnFun)

}


phenoConv <- function(x){

  if(!is.numeric(x)) stop("x should be a numeric vector")

##  if(any(x < -1)) stop("x should not be less than -1")

  if(any(x > 6)) stop("x should not be greater than 6")
  
  phenoLevels <- c(NA, 'VE', 'V0', 'V1', 'V2', 'V3', 'V4',
                   'V5', 'V6', 'V7', 'V8', 'V9',
                   'V10', 'V11', 'V12', 'V13', 'V14',
                   'V15', 'V16', 'V17', 'V18', 'V19',
                   'V20', 'V21', 'V22', 'V23', 'V24',
                   'V25', 'V26', 'VT', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6')

  pheno.char <- factor(length(x), levels = phenoLevels)

  ## If the crop has not been initiated yet
  wch1 <- which(x == -1)
  if(length(wch1) > 0) pheno.char[wch1] <- NA

  ## The crop has emerged
  wch2 <- which(x == 0)
  if(length(wch2) > 0) pheno.char[wch2] <- 'VE'

  ## The crop is at the vegetative stage
  wch3 <- which(x > 0 & x < 1)
  if(length(wch3) > 0){
    x.veg <- x[wch3] * 100
    pheno.char[wch3] <- paste('V', x.veg, sep = '')
  }

  ## The crop is at the reproductive stage
  wch4 <- which(x >= 1)
  if(length(wch4) > 0){
    x.rep <- x[wch4]
    pheno.char[wch4] <- paste('R', x.rep, sep = '')
  }

  pheno.char

}


laiParms <- function(lai.method=c("ind-leaf-Lizaso","spla","TT","Birch-Discontinuous", "Birch-Continuous"),
                     TTcoef = 0.003, max.lai = 6, Aex = 650, LT = 20, a1 = -5.61,
                     a2 = -0.59, k0 = 0.23, L0 = 133, LLx = 800, Lx = LLx-L0,
                     LNl = NULL, Amax = NULL, c1 = 8.35, c2 = 5.4, c3 = 1.58,
                     c4 = -2.62, a = NULL, b = NULL, x0 = NULL, f = 0.0048, g = NULL){
  lai.method <- match.arg(lai.method)

  if(lai.method == "TT") lai.method <- 0
  if(lai.method == "spla") lai.method <- 1
  if(lai.method == "ind-leaf-Lizaso") lai.method <- 2
  if(lai.method == "Birch-Discontinuous") lai.method <- 3
  if(lai.method == "Birch-Continuous") lai.method <- 4
  
  if(missing(LNl)){
    LNl <- -1
  }else{
    if(LNl < 5 | LNl > 20){
      warning("Impossible value for LNl")
    }
  }
    
  if(missing(Amax)){
    Amax <- -1
  }else{
    if(Amax < 400 | Amax > 1600){
      warning("Impossible value for Amax")
      Amax <- -1
    }
  }
  
  if(missing(a)){
    a <- -1
  }else{
    if(a < -1 | a > 0){
      warning("Impossible value for a")
      a <- -1
    }
  }
  
  if(missing(b)){
    b <- -1
  }else{
    if(b > 1 | b < 0){
      warning("Impossible value for b")
      b <- -1
    }
  }
  
  if(missing(x0)){
    x0 <- -1
  }else{
    if(x0 < 6 | x0 > 20){
      warning("Impossible value for x0")
      x0 <- -1
    }
  }
  
  if(missing(g)){
    g <- -1
  }else{
    if(g < 500 | g > 3000){
      warning("Impossible value for g")
      g <- -1
    }
  }
  
  list(lai.method = lai.method, TTcoef = TTcoef, max.lai = max.lai,
       Aex = Aex, LT = LT, a1 = a1, a2 = a2, k0 = k0, L0 = L0, LLx = LLx,
       Lx = Lx, LNl = LNl, Amax = Amax, c1 = c1, c2 = c2, c3 = c3, c4 = c4,
       a = a, b = b, x0 = x0, f = f, g = g)

}

MaizeCAllocParms <- function(kStem1 = 0.45, kLeaf1 = 0.45, kRoot1 = 0.10,
                             kStem2 = 0.60, kLeaf2 = 0.30, kRoot2 = 0.10,
                             kStem3 = 0.70, kLeaf3 = 0.25, kRoot3 = 0.05,
                             kStem4 = 0.000, kLeaf4 = 0.000, kRoot4 = 0.000, kGrain4 = 1){

  list(kStem1 = kStem1, kLeaf1 = kLeaf1, kRoot1 = kRoot1,
       kStem2 = kStem2, kLeaf2 = kLeaf2, kRoot2 = kRoot2,
       kStem3 = kStem3, kLeaf3 = kLeaf3, kRoot3 = kRoot3,
       kStem4 = kStem4, kLeaf4 = kLeaf4, kRoot4 = kRoot4, kGrain4 = kGrain4)

}


## Colors Stem, Leaf, Root, Rhizome, LAI

plot.MaizeGro <- function (x, obs = NULL, stem = TRUE, leaf = TRUE, root = TRUE, 
                         LAI = TRUE, grain = TRUE,
                         xlab=NULL,ylab=NULL,
                         pch=21, lty=1, lwd=1,
                         col=c("blue","green","red","black","purple"),
                         x1=0.1,y1=0.8,plot.kind=c("DB","SW","LAI","pheno","ET","cumET"),...) 
{

  if(missing(xlab)){
    xlab = expression(paste("Thermal Time (",degree,"C d)"))
  }

  if(missing(ylab)){
    ylab = expression(paste("Dry Biomass (Mg ",ha^-1,")"))
  }  

  sim <- x
  Ymax <- max(c(max(sim$Stem, na.rm=TRUE)), c(max(sim$Grain, na.rm=TRUE)), c(max(sim$Leaf, na.rm=TRUE))) + 2
  pchs <- rep(pch,length=5)
  ltys <- rep(lty,length=5)
  cols <- rep(col,length=5)
  lwds <- rep(lwd,length=5)
  plot.kind <- match.arg(plot.kind)
  if(plot.kind == "DB"){
  if (missing(obs)) {
        sim <- x
        plot1 <- xyplot(sim$Stem ~ sim$TTTc, type = "l", ...,
                        ylim = c(0, Ymax),
                        xlab = xlab,
                        ylab = ylab, 
                        panel = function(x, y, ...) {
                          if (stem == TRUE) {
                            panel.xyplot(sim$TTTc, sim$Stem, col = cols[1], 
                                         lty = ltys[1], lwd = lwds[1],...)
                          }
                          if (leaf == TRUE) {
                            panel.xyplot(sim$TTTc, sim$Leaf, col = cols[2], 
                                         lty = ltys[2], lwd = lwds[2],...)
                          }
                          if (root == TRUE) {
                            panel.xyplot(sim$TTTc, sim$Root, col = cols[3], 
                                         lty=ltys[3], lwd = lwds[3],...)
                          }
                          if (grain == TRUE) {
                            panel.xyplot(sim$TTTc, sim$Grain, col = cols[4], 
                                         lty=ltys[4], lwd = lwds[4],...)
                          }
                          if (LAI == TRUE) {
                            panel.xyplot(sim$TTTc, sim$LAI, col = cols[5], 
                                         lty = ltys[5], lwd = lwds[5],...)
                          }

                        }, key = list(text = list(c("Stem", "Leaf", "Root", 
                                        "Grain", "LAI")), col = cols, lty = ltys, lwd = lwds,
                             lines = TRUE, x = x1, y = y1))
        print(plot1)
    }
    else {
      if(ncol(obs) != 6)
        stop("obs should have 6 columns")
      sim <- x
      plot1 <- xyplot(sim$Stem ~ sim$TTTc, ..., ylim = c(0,Ymax),
                      xlab = xlab,
                      ylab = ylab, 
                      panel = function(x, y, ...) {
                        if (stem == TRUE) {
                          panel.xyplot(sim$TTTc, sim$Stem, col = cols[1], 
                                       lty = ltys[1], lwd = lwds[1], type = "l", ...)
                        }
                        if (leaf == TRUE) {
                          panel.xyplot(sim$TTTc, sim$Leaf, col = cols[2], 
                                       lty = ltys[2], lwd = lwds[2], type = "l", ...)
                        }
                        if (root == TRUE) {
                          panel.xyplot(sim$TTTc, sim$Root, col = cols[3], 
                                       lty = ltys[3], lwd = lwds[3], type = "l", ...)
                        }
                        if (grain == TRUE) {
                          panel.xyplot(sim$TTTc, sim$Grain, col = cols[4], 
                                       lty = ltys[4], lwd = lwds[4], type = "l", ...)
                        }
                        if (LAI == TRUE) {
                          panel.xyplot(sim$TTTc, sim$LAI, col = cols[5], 
                                       lty = ltys[5], lwd = lwds[5], type = "l", ...)
                        }
                        
                        panel.xyplot(obs[, 1], obs[, 2], col = cols[1], 
                                     pch=pchs[1],...)
                        panel.xyplot(obs[, 1], obs[, 3], col = cols[2], 
                                     pch=pchs[2],...)
                        panel.xyplot(obs[, 1], obs[, 4], col = cols[3], 
                                     pch=pchs[3],...)
                        panel.xyplot(obs[, 1], obs[, 5], col = cols[4], 
                                     pch=pchs[4],...)
                        panel.xyplot(obs[, 1], obs[, 6], col = cols[5], 
                                     pch=pchs[5],...)
                      }, key = list(text = list(c("Stem", "Leaf", "Root", 
                                      "Grain", "LAI")), col = cols, lines = TRUE, points=TRUE,
                           lty = ltys, pch = pchs, lwd = lwds, x = x1, y = y1))
      print(plot1)
    }
}else
  if(plot.kind == "SW"){
    matplot(x$DayofYear,as.matrix(x$cwsMat),type="l",ylab="Soil Water Content",xlab="Thermal Time")
  }else
  if(plot.kind == "LAI"){
    matplot(x$TTTc,as.matrix(x$LAImat),type="l",ylab="Individual Leaf Area",xlab="Thermal Time")
  }else
  if(plot.kind == "pheno"){
    xyplot(x$PhenoStage ~ x$TTTc,type="l",ylab="Phenological Stage",xlab="Thermal Time",...)
  }else
  if(plot.kind == "ET"){
    ## First summarize by day
    tmp <- aggregate(x$CanopyTrans, by = list(doy = x$DayofYear), FUN = sum)
    et <- tmp$x * 1e-1
    ## Soil evaporation
    tmp2 <- aggregate(x$SoilEvaporation, by = list(doy = x$DayofYear), FUN = sum)
    se <- tmp2$x * 1e-1
    xyplot(et + se ~ tmp$doy,
           ylab = "Daily EvapoTranspiration (mm)",
           xlab = "DOY",
           key = simpleKey(text = c("transp", "evapo")), ...)
  }else
  if(plot.kind == "cumET"){
    ## First summarize by day
    tmp <- aggregate(x$CanopyTrans, by = list(doy = x$DayofYear), FUN = sum)
    et <- tmp$x * 1e-1
    cumet <- cumsum(et)
    ## Soil evaporation
    tmp2 <- aggregate(x$SoilEvaporation, by = list(doy = x$DayofYear), FUN = sum)
    se <- tmp2$x * 1e-1
    cumse <- cumsum(se)
    xyplot(cumet + cumse + I(cumet + cumse) ~ tmp$doy, type="l",
           ylab = "Cummulative EvapoTranspiration (mm)",
           xlab = "DOY",
           key = simpleKey(text = c("transp", "evapo", "ET")), ...)
  }
  
}
      
print.MaizeGro <- function(x,level=1,...){

  if(level == 0){
    print(summary(as.data.frame(unclass(x)[1:11])))
  }else
  if(level == 1){
    print(summary(as.data.frame(unclass(x)[c(1,3,4,7:11)])))
  }else
  if(level == 2){
    print(summary(as.data.frame(unclass(x)[c(1,3,4,7:11,13,14)])))
  }else
  if(level == 3){
    print(summary(as.data.frame(unclass(x)[c(1,3,4,7:11,13,14,20,21,22,23)])))
  }

}

## idea for plotting canopy assimilation
## levelplot(CanopyAssim ~ DayofYear*Hour, data = res, col.regions=rainbow(10, s=1, v=0.7)[10:1], cuts=9,xlim=c(150,230))

## levelplot(CanopyAssim ~ DayofYear*Hour, data = res, col.regions=rainbow(10, s=1, v=1)[10:1], cuts=9,xlim=c(150,230))
