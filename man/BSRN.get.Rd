\name{BSRN.get}
\alias{BSRN.get}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get BSRN data
}
\description{
Download Baseline Solar Radiation Network(BSRN) data in .dat.gz format.
}
\usage{
BSRN.get(station, mmyy, directory, user, pwd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{station}{a three letter character string indicating the station name, e.g., "\code{ale}", "\code{asp}", "\code{bar}", "\code{ber}", etc.}
  \item{mmyy}{a vector of character strings indicating the month, e.g., for 2000 June, "\code{0600}"}
  \item{directory}{the directory for the downloads, default to current working directory}
  \item{user}{BSRN ftp user name}
  \item{pwd}{BSRN ftp password}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
One needs to make sure that the station-to-archive files are available. The user name and password can be obtained from Amelie Driemel \url{Amelie.Driemel@awi.de}.
}
\value{
A .dat.gz file, or .dat.gz files if \code{length(mmyy) > 1}, saved into your intended directory
}
\author{
D. Yang
}
\seealso{
\code{\link{BSRN.list}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{get}% use one of  RShowDoc("KEYWORDS")
\keyword{BSRN}% __ONLY ONE__ keyword per line
