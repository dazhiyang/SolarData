\name{OSMG.loc}
\alias{OSMG.loc}
\docType{data}
\title{
Location data for the OSMG
}
\description{
Location metadata for the Oahu Solar Measurement Gird data
}
\usage{data("OSMG.loc")}
\format{
  A data frame with 20 observations on the following 4 variables.
  \describe{
    \item{\code{Location}}{a factor with levels \code{AP1} \code{AP2} \code{AP3} \code{AP4} \code{AP5} \code{AP6} \code{AP6T} \code{AP7} \code{DH1} \code{DH10} \code{DH11} \code{DH1T} \code{DH2} \code{DH3} \code{DH4} \code{DH5} \code{DH6} \code{DH7} \code{DH8} \code{DH9}}
    \item{\code{Pakbus}}{a factor with levels \code{201} \code{202} \code{203} \code{204} \code{205} \code{206} \code{207} \code{208} \code{208T} \code{209} \code{209T} \code{210} \code{211} \code{212} \code{213} \code{214} \code{215} \code{216} \code{217} \code{230}}
    \item{\code{Latitude}}{a numeric vector}
    \item{\code{Longitude}}{a numeric vector}
  }
}
\details{
These metadata include 17 LI-200 horizontal GHI pyranometers, 2 tilted paranometers, as well as a rotating shadowband radiometer that measures all three component, namely, GHI, DNI, and DHI (or DIF, diffuse horizontal irradiance)
}
\source{
\url{https://midcdmz.nrel.gov/oahu_archive/}
}
\examples{
data(OSMG.loc)
}
\keyword{OSMG}
