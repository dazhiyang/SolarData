\name{OSMG.read}
\alias{OSMG.read}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read and process the OSMG dataset
}
\description{
This function reads and concatenates the daily OSMG files.
}
\usage{
OSMG.read(files, directory.LI200, directory.RSR = NULL,
          clear.sky = FALSE, AP2 = FALSE, agg = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{character strings indicating the file names for reading}
  \item{directory.LI200}{the directory for the LI-200 files}
  \item{directory.RSR}{the directory for the LI-200 files}
  \item{clear.sky}{boolean, indicating whether the Ineichen--Perez clera-sky irradiance should be calculated}
  \item{AP2}{boolean, indicating whether 3-second data from AP2 should be joint to the 1-second data from other stations}
  \item{agg}{numberic, aggregation interval in seconds}
}
\details{
Since the RSR files and the LI-200 files have the same file names, it is advised to store these two group of files in two separate folders, and thus I define two path variables for this function.
}
\value{
a tibble of the read and aggregated data.
}
\author{
D. Yang
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{OSMG.loc}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{OSMG}% use one of  RShowDoc("KEYWORDS")
\keyword{read}% __ONLY ONE__ keyword per line
