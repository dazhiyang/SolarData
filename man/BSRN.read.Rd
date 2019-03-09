\name{BSRN.read}
\alias{BSRN.read}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read and process the BSRN dataset
}
\description{
This function reads and concatenates the daily BSRN files.
}
\usage{
BSRN.read(file, directory, use.qc = TRUE, test = NULL, use.agg = FALSE, agg = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{character strings indicating the file names for reading}
  \item{directory}{the directory for files}
  \item{use.qc}{boolean, indicating whether the built-in QC should be used}
  \item{test}{vector of character strings, indicating which QC tests are used. Options include "phy", "ext", "closr", "df", "clim", and "all".}
  \item{use.agg}{boolean, indicating whether aggregation needs to be performed. If \code{TRUE}, \code{agg} must be stated.}
  \item{agg}{numeric, aggregation interval in minute}
}
\details{
BSRN files are monthly station-to-archive files. This function processes one file at a time.
}
\value{
a tibble of the read and aggregated data.
}
\author{
D. Yang
}
\seealso{
\code{\link{BSRN.loc}}
}
\keyword{BSRN}% use one of  RShowDoc("KEYWORDS")
\keyword{read}% __ONLY ONE__ keyword per line

