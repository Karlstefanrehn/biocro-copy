## This functions illustrate the effect of soil water content in the soil

wtrstr <- function(precipt,evapo,cws,soildepth,fieldc,wiltp,phi1=0.01,phi2 =10, smthresh=0.3, wsFun = c("linear","logistic","exp","none","thresh")){

  wsFun <- match.arg(wsFun)
  ## precipitation data are entered in mm but need to convert to m
  precipM = precipt * 1e-3;
  ## precipitation in meters can now be added to avaliable water
  ## to indicate the new avaliable water
  aw = precipM + cws;

  ## If water exceeds saturation it is considered to runoff
  if(aw > 0.5){ 
       runoff = aw - 0.5;
       aw = 0.5;
     }

  ## available water in 1 ha
  awha = aw * 1e4 * soildepth;
  ## Since evaporation is entered in Mg H2O ha-1
  ## This substraction can be made
  Newawha = awha - evapo;

  ## Need to convert back to new available water
  naw = Newawha * 1e-4 * (1/soildepth);
  ## I impose a limit of water which is the wilting point
  if(naw < 0) naw = 0;

  ## This is an empirical way of simulating water stress
  if(wsFun == "linear"){
    slp = 1/(fieldc - wiltp);
    intcpt = 1 - fieldc / (fieldc - wiltp);
    wsPhoto = slp * naw + intcpt ;
  }else
  if(wsFun == "logistic"){
    phi10 = (fieldc + wiltp)/2;
    wsPhoto = 1/(1 + exp((phi10 - naw)/phi1));
  }else
  if(wsFun == "exp"){
    slp = (1 - wiltp)/(fieldc - wiltp);
    intcpt = 1 - fieldc * slp;
    theta = slp * naw + intcpt ;
    wsPhoto = (1 - exp(-1 * (theta - wiltp)/(1 - wiltp))) / (1 - exp(-1));
  }else
  if(wsFun == "thresh"){
    rawc = (naw - wiltp)/(fieldc - wiltp)
    if(rawc > smthresh){
      wsPhoto = 1
    }else{
      wsPhoto = rawc / smthresh
    }
  }else{
    wsPhoto = 1;
  }
    
  if(wsPhoto < 0) wsPhoto = 0;
  if(wsPhoto > 1) wsPhoto = 1
  wsSpleaf = naw^phi2 * 1/fieldc^phi2;
  if(wsSpleaf > 1) wsSpleaf = 1;

  list(wsPhoto=wsPhoto,wsSpleaf=wsSpleaf,awc=naw)

}

wsRcoef <- function(aw,fieldc,wiltp,phi1,phi2, smthresh, wsFun = c("linear","logistic","exp","none","thresh") ){

  wsFun <- match.arg(wsFun)

  if(aw > fieldc){
    aw = fieldc
    warning("In this function available water should be lower than field capacity")
  }
  
  if(wsFun == "linear"){
    slp = 1/(fieldc - wiltp);
    intcpt = 1 - fieldc / (fieldc - wiltp);
    wsPhoto = slp * aw + intcpt ;
  }else
  if(wsFun == "logistic"){
    phi10 = (fieldc + wiltp)/2;
    wsPhoto = 1/(1 + exp((phi10 - aw)/phi1));
  }else
  if(wsFun == "exp"){
    slp = (1 - wiltp)/(fieldc - wiltp);
    intcpt = 1 - fieldc * slp;
    theta = slp * aw + intcpt ;
    wsPhoto = (1 - exp(-2.5 * (theta - wiltp)/(1 - wiltp))) / (1 - exp(-2.5));
  }else
  if(wsFun == "thresh"){
    raw = (aw - wiltp)/(fieldc - wiltp)
    if(raw > smthresh){
      wsPhoto = 1
    }else{
      wsPhoto = raw / smthresh
    }
  }else{
    wsPhoto = 1;
  }

  wsPhoto[wsPhoto < 0] <- 0
  wsSpleaf <- aw^phi2 * 1/fieldc^phi2
  wsSpleaf[wsSpleaf > 1] <- 1
  list(wsPhoto=wsPhoto,wsSpleaf=wsSpleaf)
}

