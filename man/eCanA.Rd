\name{eCanA}
\alias{eCanA}
\title{Simulates canopy assimilation (von Caemmerer model)}

\description{
  It represents an integration of the photosynthesis function
  \code{\link{eC4photo}}, canopy evapo/transpiration and the
  multilayer canopy model \code{\link{sunML}}.
}
\usage{
eCanA(LAI, doy, hour, solarR, AirTemp, RH, WindS, Vcmax, 
    Vpmax, Vpr, Jmax, Ca = 380, Oa = 210, StomataWS = 1)
}

\arguments{

  
\item{LAI}{leaf area index.}
\item{doy}{day of the year, (1--365).}
\item{hour}{hour of the day, (0--23).}
\item{solarR}{solar radiation (\eqn{\mu mol \; m^{-2} \; s^{-1}}{micro mol /m2/s}).}
\item{AirTemp}{temperature (Celsius).}
\item{RH}{relative humidity (0--1).}
\item{WindS}{wind speed (\eqn{m \; s^{-1}}{m/s}).}
\item{Vcmax}{Maximum rubisco activity (\eqn{\mu mol \; m^{-2} \; s^{-1}}{micro mol /m2/s}).}
\item{Vpmax}{Maximum PEP carboxylase activity (\eqn{\mu mol \; m^{-2} \; s^{-1}}{micro mol /m2/s}).}
\item{Vpr}{PEP regeneration rate  (\eqn{\mu mol \; m^{-2} \; s^{-1}}{micro mol /m2/s}).}
\item{Jmax}{Maximal electron transport rate (\eqn{\mu}{micro}mol electrons
  \eqn{m^{-2}}{m-2} \eqn{s^{-1}}{s-1}).}
\item{Ca}{atmospheric carbon dioxide concentration (ppm or \eqn{\mu}{micro}bar) (e.g. 380).}
\item{Oa}{atmospheric oxygen concentration (mbar) (e.g. 210).}
\item{StomataWS}{ Effect of water stress on assimilation. }
}

\value{

  \code{\link{numeric}} 

  returns a single value which is hourly
  canopy assimilation (mol \eqn{m^{-2}}{m-2} ground \eqn{hr^{-1}}{s-1})
}    

\examples{

## No examples yet

}

\keyword{models}