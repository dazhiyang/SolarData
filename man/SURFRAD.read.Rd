\name{SURFRAD.read}
\alias{SURFRAD.read}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read and process the SURFRAD dataset
}
\description{
This function reads and concatenates the daily SURFRAD files.
}
\usage{
SURFRAD.read(files, directory, use.original.qc = FALSE, use.qc = TRUE,
             test = NULL, progress.bar = TRUE, agg = 1, additional.variables = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{character strings indicating the file names for reading}
  \item{directory}{the directory for files}
  \item{use.original.qc}{boolean, indicating whether the SURFRAD QC flags should be used}
  \item{use.qc}{boolean, indicating whether the built-in QC should be used}
  \item{test}{vector of character strings, indicating which QC tests are used. Options include "phy", "ext", "closr", "df", "clim", and "all".}
  \item{progress.bar}{boolean, indicating whether the progress bar should be used}
  \item{agg}{numeric, aggregation interval in minute}
  \item{additional.variables}{four compulsory variables are included, namely, "dw_solar", "direct_n", "diffuse", and "pressure". This function parameter allows the user to set additional variables for extraction, such as "temp" or "windspd". See the code for the available variable list.}
}
\details{
It is not recommended to put all files in a same folder. The SURFRAD FTP folder directory is good, keep that, and use a loop (outside of this function) to access all folders.
}
\value{
a tibble of the read and aggregated data.
}
\author{
D. Yang
}
\seealso{
\code{\link{SURFRAD.loc}}, \code{\link{SURFRAD.get}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{SURFRAD}% use one of  RShowDoc("KEYWORDS")
\keyword{read}% __ONLY ONE__ keyword per line
