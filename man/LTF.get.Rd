\name{LTF.get}
\alias{LTF.get}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Read Linke turbidity factor from tiff maps
}
\description{
Read monthyl Linke turbidity factor from 12 tiff maps for any location(s) in the world.
}
\usage{
LTF.get(lon, lat, directory)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{lon}{a number or a numeric array, longitude of the location(s) for downloading}
  \item{lat}{a number or a numeric array, latitude of the location(s) for downloading}
  \item{directory}{the directory of the \emph{downloaded} .tiff images}
}
\details{
Twelve monthly maps of Linke turbidity factor values were created given by latitude and longitude. Latitude is positive North, longitude is positive eastwards of longitude 0. The data are in gridded, raw format, no header (tiff images), 1 byte per value (unsigned int encoding), 2160 rows and 4320 columns. Cell size is 5' (approx. 10 km at mid-latitude). Upper left corner is 90 N, 180 W. Then, point 90 N, 179.5 W etc. Lower right is 90 S, 180 E.
}
\value{
A \eqn{n\times 12} matrix, where \eqn{n} is the number of lat-lon pairs.
}
\author{
D. Yang
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{get}% use one of  RShowDoc("KEYWORDS")
\keyword{LTF}% __ONLY ONE__ keyword per line
