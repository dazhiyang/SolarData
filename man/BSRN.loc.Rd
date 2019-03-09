\name{BSRN.loc}
\alias{BSRN.loc}
\docType{data}
\title{
Location metadata for BSRN
}
\description{
BSRN location metadata
}
\usage{data("BSRN.loc")}
\format{
  A data frame with 76 rows (stations) on the following 19 variables.
  \describe{
    \item{\code{stn}}{a character vector for station abbreviation}
    \item{\code{full}}{a character vector for station full name}
    \item{\code{lat}}{a numeric vector, in degress}
    \item{\code{lon}}{a numeric vector, in degress}
    \item{\code{elev}}{a numeric vector, in meters}
    \item{\code{status}}{a character vector showing the status of the stations, which can be either \code{Active}, \code{Closed}, or \code{Candidate}.}
    \item{\code{tz}}{a character vector}
    \item{\code{LTF.Jan}}{a numeric vector}
    \item{\code{LTF.Feb}}{a numeric vector}
    \item{\code{LTF.Mar}}{a numeric vector}
    \item{\code{LTF.Apr}}{a numeric vector}
    \item{\code{LTF.May}}{a numeric vector}
    \item{\code{LTF.Jun}}{a numeric vector}
    \item{\code{LTF.Jul}}{a numeric vector}
    \item{\code{LTF.Aug}}{a numeric vector}
    \item{\code{LTF.Sep}}{a numeric vector}
    \item{\code{LTF.Oct}}{a numeric vector}
    \item{\code{LTF.Nov}}{a numeric vector}
    \item{\code{LTF.Dec}}{a numeric vector}
  }
}
\details{
This dataset contains the metadata of the 76 BSRN stations. To facilitating fast retrieval of Linke turbidity values at these locations, they are preloaded.
}

\keyword{datasets}
\keyword{BSRN}
