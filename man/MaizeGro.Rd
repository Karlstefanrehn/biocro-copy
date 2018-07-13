\name{MaizeGro}
\alias{MaizeGro}
\alias{MaizeCAllocParms}
\alias{MaizeNitroParms}
\alias{MaizePhenoParms}
\alias{MaizePhotoParms}
\alias{MaizeSeneParms}
\alias{phenoConv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Simulation of Maize, Growth, LAI, Photosynthesis and phenology
}
\description{

  It takes weather data as input (hourly timesteps) and several
  parameters
  and it produces phenology, photosynthesis, LAI, etc.

}
\usage{
MaizeGro(WetDat, plant.day = NULL, emerge.day = NULL, harvest.day = NULL, 
    plant.density = 7, timestep = 1, lat = 40, canopyControl = list(), 
    MaizeSeneControl = list(), photoControl = list(), MaizePhenoControl = list(), 
    MaizeCAllocControl = list(), laiControl = list(), soilControl = list(), 
    MaizeNitroControl = list(), centuryControl = list()) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{WetDat}{weather data as produced by the
    \code{\link{weach}} function.}

  \item{plant.day}{
    Planting date (format 0-365)
}
  \item{emerge.day}{
    Emergence date (format 0-365)
}
  \item{harvest.day}{
    Harvest date (format 0-365)
}
  \item{plant.density}{
    Planting density (plants per meter squared, default = 7)
  }

  \item{timestep}{Simulation timestep, the default of 1 requires houlry
  weather data. A value of 3 would require weather data every 3 hours.
  This number should be a divisor of 24.}

  \item{lat}{
    latitude, default 40.
  }

  \item{canopyControl}{List that controls aspects of the canopy simulation. It should be supplied through the \code{canopyParms} function. 

  
\code{Sp} (specific leaf area) here the units are ha \eqn{Mg^{-1}}.
  If you have data in \eqn{m^2} of leaf per kg of dry matter (e.g. 15)
  then divide by 10 before inputting this coefficient.

  \code{SpD} decrease of specific leaf area. Empirical
  parameter. Default 0. example value (1.7e-3).

\code{nlayers} (number of layers of the canopy) 
  Maximum 50. To increase the number of layers (more than 50) the \code{C}
  source code needs to be changed slightly.

\code{kd} (extinction coefficient for diffuse light) between 0 and 1.

\code{mResp} (maintenance respiration) a vector of length 2 with the first component for leaf and stem and the second component for rhizome and root.

 For other terms see \code{\link{BioGro}}.
}


\item{MaizeSeneControl}{ List that controls aspects of senescence simulation. It should be supplied through the \code{MaizeSeneParms} function.

\code{senLeaf} Thermal time at which leaf senescence will start.

\code{senStem} Thermal time at which stem senescence will start.

\code{senRoot} Thermal time at which root senescence will start.

}

\item{photoControl}{  List that controls aspects of photosynthesis simulation. It should be supplied through the \code{MaizePhotoParms} function.

\code{vmax} Vmax passed to the \code{\link{c4photo}} function.

\code{alpha} alpha parameter passed to the \code{\link{c4photo}}
  function.

\code{kparm} kparm parameter passed to the \code{\link{c4photo}}
function.

\code{theta} theta parameter passed to the \code{\link{c4photo}}
function.

\code{beta} beta parameter passed to the \code{\link{c4photo}}
function.

\code{Rd} Rd parameter passed to the \code{\link{c4photo}}
  function.

\code{Catm} Catm parameter passed to the \code{\link{c4photo}}
  function.

\code{b0} b0 parameter passed to the \code{\link{c4photo}}
  function.

\code{b1} b1 parameter passed to the \code{\link{c4photo}}
  function.

}


\item{MaizePhenoControl}{
      argument used to pass parameters related to phenology characteristics
%%     ~~Describe \code{MaizePhenoControl} here~~
}

\item{MaizeCAllocControl}{
Carbon allocation partitioning for maize. Similar to \code{\link{BioGro}}.
  }

  \item{soilControl}{
 See function \code{\link{BioGro}} for details.
}
\item{MaizeNitroControl}{
  See \code{MaizeNitroParms}
}
  \item{centuryControl}{
  See function \code{\link{BioGro}} for details. 
}

\item{laiControl}{
  argument used to pass parameters related to leaf area index. Use the
  function laiParms for this. The additional parameters depend on the
  method used.

  For the "TT" method the leaf area index only depends on the thermal
  time and the relevant parameters are TTcoef and max.lai.

  The method "ind-leaf-Lizaso" has been tested more and the reference
  documentation is the Lizaso et al. (2003) paper. It is currently the
  default. At the moment it is not sensitive to stress. This need to be
  fixed.
  
  }

}
\details{

  The phenology follows the 'Corn Growth and Development' Iowa State Publication.

}
\value{

  It currently returns a list with the following components

  \item{DayofYear}{Day of the year (1-365)}

  \item{Hour}{Hour of the day (0-23)}

  \item{TTTc}{Accumulated thermal time}

  \item{PhenoStage}{Phenological stage of the crop}

  \item{CanopyAssim}{Hourly canopy assimilation,
    (Mg \eqn{ha^-1} ground \eqn{hr^-1}).}

  \item{CanopyTrans}{Hourly canopy transpiration,
    (Mg \eqn{ha^-1} ground \eqn{hr^-1}).}

  \item{LAI}{Leaf Area Index}

  \item{Leaf}{Leaf biomass (Mg ha-1)}

  \item{Stem}{Stem biomass (Mg ha-1)}

  \item{Root}{Root biomass (Mg ha-1)}

  \item{Grain}{Grain biomass - yield - (Mg ha-1)}

  \item{LAImat}{Leaf Area Index matrix - leaf area index by leaf - (m^2
    m^-2)}

  \item{VmaxVec}{Vmax vector - it changes if a relationship with
    Nitrogen is defined}

  \item{LeafNVec}{Leaf N concentration (g m^-2)}

  \item{cwsMat}{current water status matrix - this will have more than
    one column if the soil model is multilayered (m^3 m^-3)}

  \item{psimMat}{soil water potential (MPa)}

  \item{rdMat}{root density matrix (Mg ha^-1)}

  \item{SCpools}{Soil carbon pools }

  \item{SNpools}{soil nitrogen pools }

  \item{LeafPsimVec}{Leaf water potential (MPa)}

  \item{SoilEvaporation}{soil evaporation (Mg ha-1). plotting with the
    plot(x,plot.kind="ET") will use mm/day units.}

  \item{SoilWatCont}{Average soil water content for all layers
    (volumetric)}

  \item{SomatalCondCoefs}{Stomatal conductance coefficient (0-1)}

}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Fernando E Miguez
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{BioGro}}
  Use plotting methods such as
  plot(x, plot.kind="DB") for dry biomass
  plot(x, plot.kind="SW") for soil water
  plot(x, plot.kind="ET") for EvapoTranspiration
  plot(x, plot.kind="cumET") for cummulative EvapoTranspiration
  plot(x, plot.kind="LAI") for Leaf Area Index
  plot(x, plot.kind="pheno") for phenology

}
\examples{
\dontrun{
data(cmi05)
res <- MaizeGro(cmi05, plant.day = 110,
                emerge.day = 120, harvest.day=300)
                

## Explore data with plots
plot(res, plot.kind="DB") ## for dry biomass
plot(res, plot.kind="SW") ## for soil water
plot(res, plot.kind="ET") ## for EvapoTranspiration
plot(res, plot.kind="cumET") ## for cummulative EvapoTranspiration
plot(res, plot.kind="LAI") ## for Leaf Area Index
plot(res, plot.kind="pheno") ## for phenology

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
