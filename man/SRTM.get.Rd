\name{SRTM.get}
\alias{SRTM.get}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get NASA SRTM data
}
\description{
Download NASA Shuttle Radar Topography Mission (SRTM) digital elevation model (DEM) data in .hgt format.
}
\usage{
SRTM.get(resolution, files, directory = "data-raw")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resolution}{a number indicating the resolution of the DEM, two options are available: 1 or 3, see \code{\link{SRTM.list}}}
  \item{files}{character strings indicating the file names for downloading, see \code{\link{SRTM.list}} for available files}
  \item{directory}{the directory for the downloads, default to "\code{~/data-raw}"}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
place holder for now
}
\value{
A .dat file, or .dat files if \code{length(day_of_year) > 1}, saved into your intended directory
}
\references{
Bernhard Rabus, Michael Eineder, Achim Roth, Richard Bamler, The shuttle radar topography mission—a new class of digital elevation models acquired by spaceborne radar, \emph{ISPRS Journal of Photogrammetry and Remote Sensing}, Volume 57, Issue 4, 2003, Pages 241-262, \url{https://doi.org/10.1016/S0924-2716(02)00124-7}.
}
\author{
D. Yang
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{get}% use one of  RShowDoc("KEYWORDS")
\keyword{SRTM}% __ONLY ONE__ keyword per line
