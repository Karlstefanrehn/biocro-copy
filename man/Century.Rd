\name{Century}
\alias{Century}
\alias{CenturyC}
\alias{somc}
\alias{isom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ This function implements the Century model from Parton. }
\description{
Calculates flows of soil organic carbon and nitrogen based on the
Century model. There are two versions one written in R (Century) and one
in C (CenturyC). The C version only
runs at weekly timesteps.
}
\usage{
Century(LeafL, StemL, RootL, RhizL, smoist, stemp, precip, leachWater,
centuryControl = list(), soilType=5, verbose = FALSE)

CenturyC(LeafL, StemL, RootL, RhizL, smoist, stemp, precip, 
    leachWater, centuryControl = list(), soilType = 5) 

somc(om = 2, cc = 0.5, depth = 0.3, bd=1.3, pp=c(0.01,0.19,0.8))
isom(SC6=0.39, SC7=7.41, SC8=31.2, cc=0.5, depth=0.3, bd=1.3)

}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LeafL}{ Leaf litter.}
  \item{StemL}{ Stem litter. }
  \item{RootL}{ Root litter. }
  \item{RhizL}{ Rhizome litter. }
  \item{smoist}{ Soil moisture (0-1). }
  \item{stemp}{ Soil temperature (Celsius). } 
  \item{precip}{ Precipitation (mm). }
  \item{leachWater}{ Leached water. }
  \item{centuryControl}{ See \code{\link{centuryParms}}.}
  \item{verbose}{ Only used in the R version for debugging. }
  \item{soilType}{ See \code{\link{showSoilType}}.}

  function somc
  
  \item{om}{organic matter (percent)}
  \item{cc}{carbon content (proportion 0-1) of the om}
  \item{depth}{soil depth (meters)}
  \item{bd}{bulk density (g/cm3)}
  \item{pp}{partition proportion, the firts one is the proportion of the
    soil organic carbon in microbes, second in slow and third in passive.}

  function isom

  Some arguments are the same as somc
  
  \item{SC6}{microbial biomass (Mg/ha)}
  \item{SC7}{slow pool (Mg/ha)}
  \item{SC8}{passive pool (Mg/ha)}

}
\details{
 Most of the details can be found in the papers by Parton et al. (1987,
 1988, 1993)

 Function somc can used to set the values of the main soil organic
 carbon pools. For an application the model would need to be spined (run
 for a long time).

 Function isom is for doing the inverse conversion that somc does.
 
 At this point the R version and the C version will differ
 slightly. TODO is to fix this. At least part of the reason is that the
 texture is different. 
}
\value{
  A list with,
  \item{SCs }{Soil carbon pools 1-9.}
  \item{SNs }{Soil nitrogen pools 1-9.}
  \item{MinN}{Mineralized nitrogen.}
  \item{Resp}{Soil respiration.}
}
\references{ ~put references to the literature/web site here ~ }
\author{ Fernando E. Miguez }
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
% \seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
% Add one or more standard keywords, see file 'KEYWORDS' in the
Century(0,0,0,0,0.3,5,2,2)$Resp
Century(0,0,0,0,0.3,5,2,2)$MinN


}
\keyword{ models }

