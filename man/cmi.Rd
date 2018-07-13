\name{cmi-weather}
\alias{weather06}
\alias{weather05}
\alias{weather04}
\alias{cmi05}
\alias{cmi04}
\alias{cmi06}
\docType{data}
\title{
Weather data for Champaign Illinois
}
\description{
    Weather data as produced by the \code{\link{weach}} function.
  These are for 2004, 2005 and 2006. The original source is daily
  weather data from the WARM network IL.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(cmi06)}
\format{
  A data frame with 8760 (8784 if leap year) observations on the following 8 variables.
  \describe{
    \item{\code{year}}{year}
    \item{\code{doy}}{day of the year}
    \item{\code{hour}}{hour}
    \item{\code{SolarR}}{solar radiation (micro mol/m2/s)}
    \item{\code{Temp}}{temperature (Celsius)}
    \item{\code{RH}}{relative humidity (0-1)}
    \item{\code{WS}}{wind speed (m s-1)}
    \item{\code{precip}}{precipitation (mm)}
  }
}
\source{
  These data are generated using the \code{weach} function.
  The original data is from the Illinois WARM network.
  http://www.isws.illinois.edu/warm/data/cdfs/cmiday.txt

  Datasets \code{weather04}, \code{weather05}, \code{weather06} were
  generated with an old version of weach which presumably did not
  produce the most realistic sub-daily profiles of realtive humidity and
  wind speed (they were constant through out the day before).
}
\examples{
data(cmi06)
summary(cmi06)
## maybe str(weather06) ; plot(weather06) ...
}
\keyword{datasets}
