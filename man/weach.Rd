\name{weach}
\alias{weach}
\alias{is.leap}
\title{Simulates the hourly conditions from daily}

\description{
  Manipulates weather data in the format obtained from WARM
  (see link below)
  and returns the format and units needed for most functions
  in this package. This function should be used for one year at a time.
  It returns hourly (or sub-daily) weather information. 

}
\usage{
weach(X,lat=40,ts=1, solar.units=c("MJ/m2"),
     temp.units=c("Fahrenheit","Celsius"),
     rh.units=c("percent","fraction"),
     ws.units=c("mph","mps"),
     pp.units=c("in","mm"),seed=1234,...)

is.leap(year)
}

\arguments{

  \item{X}{a matrix (or data frame) containing weather information.
    The input format is strict but it is meant to be used with the data
    usually obtained from weather stations in Illinois. The data frame
    should  have 11 columns (see details).}
  \item{lat}{latitude at the specific location}
  \item{ts}{timestep for the simulation of sub-daily data from
    daily. For example a value of 3 would return data every 3
    hours. Only divisors of 24 work (i.e. 1,2,3,4, etc.).}
  \item{solar.units}{ Assumed units for solar radiation input.}
  \item{temp.units}{ Option to specify the units in which the
    temperature is entered. Default is Farenheit.}
  \item{rh.units}{ Option to specify the units in which the
    relative humidity is entered. Default is percent.}
  \item{ws.units}{ Option to specify the units in which the
    wind speed is entered. Default is miles per hour.}
  \item{pp.units}{ Option to specify the units in which the
    precipitation is entered. Default is inches.}
  \item{\ldots}{additional arguments to be passed to
    \code{\link{lightME}}}
  \item{year}{ input for the \code{is.leap} function. It checks
    if a year is leap or not.}
  \item{seed}{seed used for the random generation of the distribution for
    downscaling rainfall from daily to hourly.}
}

\details{

  This function was originally used to transform daily data to hourly
  data. Some flexibility has been added so that other units can be
  used. The input data used originally looked as follows.
\enumerate{

\item{col 1}{ year}
\item{col 2}{ day of the year (1--365). Does not consider leap years.}
\item{col 3}{ total daily solar radiation (MJ/m^2).}
\item{col 4}{ maximum temperature (Fahrenheit).}
\item{col 5}{ minimum temperature (Fahrenheit).}
\item{col 6}{ average temperature (Fahrenheit).}
\item{col 7}{ maximum relative humidity (\%).}
\item{col 8}{ minimum relative humidity (\%).}
\item{col 9}{ average relative humidity (\%).}
\item{col 10}{ average wind speed (miles per hour).}
\item{col 11}{ precipitation (inches).}
}

All the units above are the defaults but they can be changed as part of
the arguments.
}




\value{
  a \code{\link{data.frame}} 
  returning hourly (or sub-daily) weather data. Dimensions 8760 (if hourly) by 8.
If it is a leap year it will be 8784 rows.
  
\item{year}{Year.}
\item{doy}{Day of the year.}
\item{hour}{Hour of the day (0--23, depending on the timestep).}
\item{SolarR}{Direct solar radiation (\eqn{\mu mol \; m^{-2} \;
    s^{-1}}{micro mol per meter squared per second}).}
\item{Temp}{Air temperature (Celsius).}
\item{RH}{Relative humidity (0--1).}
\item{WS}{Average wind speed (\eqn{m \;s^{-1}}{meter per second}).}
\item{precip}{Precipitation (\eqn{mm}{mm})}

}    

\examples{
\dontrun{
## The data typically used can be obtained from
## http://www.isws.illinois.edu/warm/weather/
## I have downloaded and cleaned up the data
## For Champaign and it is available here
url <- "http://www.agron.iastate.edu/miguezlab/research/biocro/data/cmi-weather.csv"
cmi <- read.csv(url)

cmi05 <- weach(subset(cmi, year == 2005), lat=40)
}
}
\keyword{datagen}
