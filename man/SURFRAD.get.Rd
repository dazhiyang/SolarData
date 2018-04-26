\name{SURFRAD.get}
\alias{SURFRAD.get}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get NOAA SURFRAD data
}
\description{
Download NOAA Surface Radiation (SURFRAD) data in .dat format.
}
\usage{
SURFRAD.get(station, year, day.of.year, directory = "data-raw")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{station}{a character string indicating the station name, options are "\code{Bondville_IL}", "\code{Boulder_CO}", "\code{Desert_Rock_NV}", "\code{Fort_Peck_MT}", "\code{Goodwin_Creek_MS}", "\code{Penn_State_PA}", "\code{Sioux_Falls_SD}". Alternatively, abbreviations of the stations can be used, i.e., "\code{bon}", "\code{tbl}", "\code{dra}", "\code{fpk}", "\code{gwn}", "\code{psu}", "\code{sxf}", respectively}
  \item{year}{a character string indicating the four-digit year}
  \item{day.of.year}{days of year to be downloaded, as a numeric array}
  \item{directory}{the directory for the downloads, default to "\code{~/data-raw}"}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
NOAA high-resolution, long-term, ground-based irradiance measurements at 7 locations.
}
\value{
A .dat file, or .dat files if \code{length(day_of_year) > 1}, saved into your intended directory
}
\references{
J.A. Augustine, J.J. DeLuisi, and C.N. Long, SURFRAD–A National Surface Radiation Budget Network for Atmospheric Research. \emph{Bull. Amer. Meteor. Soc.}, Volume 81, Pages 2341–2358, \url{https://doi.org/10.1175/1520-0477(2000)081<2341:SANSRB>2.3.CO;2}
}
\author{
D. Yang
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{get}% use one of  RShowDoc("KEYWORDS")
\keyword{SURFRAD}% __ONLY ONE__ keyword per line
