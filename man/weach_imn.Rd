\name{weach_imn}
\alias{weach_imn}
\alias{weach_imn2}
\alias{boo14.200}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Weather change Iowa Mesonet

}
\description{
  Manipulates weather data in the format obtained from Iowa Mesonet
  (see link below)
  and returns the format and units needed for most functions
  in this package. This function should be used for one year at a time.
  It takes and returns hourly weather information only.
}
\usage{
weach_imn(data, ts = 1, temp.units = c("Fahrenheit", "Celsius"),
rh.units = c("percent", "fraction"), ws.units = c("mph", "mps"),
pp.units = c("in", "mm"), ...)

weach_imn2(data, temp.units = c("Fahrenheit", "Celsius"), rh.units = c("percent", 
    "fraction"), ws.units = c("mph", "mps"), pp.units = c("in", 
    "mm"), na.chr = -99) 
NULL
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    data as obtained from the Iowa Mesonet (see details)
}

\item{ts}{
  Time step, at the moment it is assumed to be 1.
}
\item{temp.units}{
  Temperature units.
%%     ~~Describe \code{temp.units} here~~
}
\item{rh.units}{
  Relative humidity units
%%     ~~Describe \code{rh.units} here~~
}
\item{ws.units}{
  wind speed units
%%     ~~Describe \code{ws.units} here~~
}
\item{pp.units}{
  precipitation units
%%     ~~Describe \code{pp.units} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
\item{na.chr}{missing data character}
}
\details{
    The first function (weach_imn) should be used to transform data from the Iowa Mesonet
    at hourly intervals from here:

    http://mesonet.agron.iastate.edu/agclimate/hist/hourlyRequest.php

    When selecting to download variables:
    Air Temperature (Fahrenheit)
    Solar Radiation (kilocalories per meter squared)
    Precipitation (inches)
    Relative humidity (percent)
    Wind Speed (mph)

    You can read the data directly as it is downloaded making sure you
    skip the first 6 lines (This includes the title row).
    
    The imported data frame should have 9 columns with:
    
\enumerate{

\item{col 1}{ site ID}
\item{col 2}{ site name }
\item{col 3}{ date in format "year-month-day", e.g. '2010-3-25' }
\item{col 4}{ hour in format "hour:minute", e.g. '15:00'}
\item{col 5}{ temperature (Fahrenheit) }
\item{col 6}{ solar radiation (kilocalories per meter squared }
\item{col 7}{ precipitation (inches) }
\item{col 8}{ relative humidity (\%).}
\item{col 9}{ wind speed (mph) }

}

The second function is in development and it is meant to be used with
data downloaded from:

http://mesonet.agron.iastate.edu/agclimate/hist/hourly.php

It might not work well at the moment. It still needs work.
}
\value{

  It will return a data frame in the same format as the
  \code{\link{weach}} function.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Iowa Mesonet
  http://mesonet.agron.iastate.edu/index.phtml
}
\author{
Fernando E. Miguez
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{weach}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

\dontrun{
## Read an example data set from my website
## This can take a while as it processes data for the period 1986-2013
url <- "http://www.agron.iastate.edu/miguezlab/research/biocro/data/ames-weather.csv"
ams <- read.csv(url)
ams2 <- weach_imn(ams)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }

