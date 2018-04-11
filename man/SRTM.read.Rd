\name{SRTM.read}
\alias{SRTM.read}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read SRTM .hgt files
}
\description{
This function reads SRTM .hgt files and outputs a \code{RasterLayer} obejct or a \code{data.frame} (not recommended).
}
\usage{
SRTM.read(files, as.data.frame = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{a vector of character strings indicating the file names to be read}
  \item{as.data.frame}{boolean, if \code{TRUE} output a \code{data.frame}}
}
\details{
The SRTM1 has \eqn{3601\times 3601} cells, and SRTM3 has \eqn{1201\times 1201} cells. It is thus not recommened to convert the output into a \code{data.frame}. Moreover, the \code{RasterLayer} object is easier to work with using the \href{https://cran.r-project.org/web/packages/raster/index.html}{raster} package.
}
\value{
A \code{RasterLayer} object or a \code{data.frame}
}
\author{
D. Yang
}
\seealso{
\code{\link{SRTM.list}}, \code{\link{SRTM.get}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{read}% use one of  RShowDoc("KEYWORDS")
\keyword{SRTM}% __ONLY ONE__ keyword per line
