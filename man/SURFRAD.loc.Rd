\name{SURFRAD.loc}
\alias{SURFRAD.loc}
\docType{data}
\title{
Location metadata for SURFRAD
}
\description{
SURFRAD location metadata
}
\usage{data("SURFRAD.loc")}
\format{
  A data frame with 7 observations on the following 19 variables.
  \describe{
    \item{\code{stn}}{a factor with levels \code{bon} \code{dra} \code{fpk} \code{gcm} \code{psu} \code{sxf} \code{tbl}}
    \item{\code{full}}{a factor with levels \code{Bondville, Illinois} \code{Desert Rock, Nevada} \code{Fort Peck, Montana} \code{Goodwin Creek, Mississippi} \code{Penn. State Univ., Pennsylvania} \code{Sioux Falls, South Dakota} \code{Table Mountain, Boulder, Colorado}}
    \item{\code{lat}}{a numeric vector}
    \item{\code{lon}}{a numeric vector}
    \item{\code{elev}}{a numeric vector}
    \item{\code{tz}}{a numeric vector}
    \item{\code{install}}{a factor with levels \code{1994-04} \code{1994-11} \code{1994-12} \code{1995-07} \code{1998-03} \code{1998-06} \code{2003-06}}
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
This dataset contains the metadata of the 7 SURFRAD stations. To facilitating fast retrieval of Linke turbidity values at these locations, they are preloaded.
}

\keyword{datasets}
\keyword{SURFRAD}
