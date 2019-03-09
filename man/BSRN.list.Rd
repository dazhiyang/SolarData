\name{BSRN.list}
\alias{BSRN.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
List the available BSRN files
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function retrieves the directory listing, i.e., available files from \url{ftp://ftp.bsrn.awi.de/}.
}
\usage{
BSRN.list(station, user, pwd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{station}{a vector of text string indicating the station abbreviations}
  \item{user}{BSRN ftp user name}
  \item{pwd}{BSRN ftp password}
}
\details{
BSRN ftp listing is dynamic. This function retrieves the most updated listing of BSRN station-to-archive files. You do not need to use this function very often. It is recommended to use it yearly to check on file updates.
}
\value{
A list of tibbles, indicating the available file names.
}

\author{
D. Yang
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{BSRN.read}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{list}% use one of  RShowDoc("KEYWORDS")
\keyword{BSRN}% __ONLY ONE__ keyword per line
