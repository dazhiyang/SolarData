\name{PSM.get}
\alias{PSM.get}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Get NREL PSM version 3 data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Use API to download NREL Physical Solar Model (PSM) version 3 data in .csv format.
}
\usage{
PSM.get(lon, lat, api.key, attributes, name, affiliation, year, leap.year,
interval, utc, reason.for.use, email, mailing.list, directory = "data-raw")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lon}{a number or a numeric array, longitude of the location(s) for downloading}
  \item{lat}{a number or a numeric array, latitude of the location(s) for downloading}
  \item{api.key}{the API key as a character string, can be obtained at \url{https://developer.nrel.gov/signup/}}
  \item{attributes}{the parameters to be downloaded, options are: "\code{air_temperature}", "\code{clearsky_dhi}", "\code{clearsky_dni}", "\code{clearsky_ghi}", "\code{cloud_type}", "\code{dew_point}", "\code{dhi}", "\code{dni}", "\code{fill_flag}", "\code{ghi}", "\code{relative_humidity}", "\code{solar_zenith_angle}", "\code{surface_albedo}", "\code{surface_pressure}", "\code{total_precipitable_water}", "\code{wind_direction}", "\code{wind_speed}". This argument should be passed in as a character string. If more than one, separate them using comma without spaces.}
  \item{name}{your name as a character string, e.g., "\code{John+Smith}". Spaces need to be replace with the + sign}
  \item{affiliation}{your affiliation as a character string, e.g., "\code{National+Renewable+Energy+Lab}". Spaces need to be replace with the + sign.}
  \item{year}{the year to be downloaded, options are: "\code{1998}", "\code{1999}", \eqn{\dots}, "\code{2016}", "\code{tmy}". This argument should be passed in as a character string. The API only allows downloading one year at a time}
  \item{leap.year}{a character string ("\code{true}" or "\code{false}") indicating whether you want the data on 29 Feb if it is a leap year}
  \item{interval}{a character string ("\code{30}" or "\code{60}") indicating whether you want 30 min or 60 min data}
  \item{utc}{a character string ("\code{true}" or "\code{false}") indicating whether you want the time to be UTC}
  \item{reason.for.use}{a character string indicating your purpose, e.g., "\code{research+development}"}
  \item{email}{your email as a character string, e.g., "\code{john.smith@gmail.com}"}
  \item{mailing.list}{a character string ("\code{true}" or "\code{false}") indicating whether you want to be on the NREL mainling list}
  \item{directory}{the directory for the downloads, default to "\code{~/data-raw}"}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
NREL PSM v3 contains half-hourly, regularly-gridded, satellite-derived irradiance and other meteorological parameters.
}
\value{
A .csv file, or .csv files if \code{length(lat) > 1}, saved into your intended directory
}
\references{
Manajit Sengupta, Yu Xie, Anthony Lopez, Aron Habte, Galen Maclaurin, and James Shelby, The National Solar Radiation Data Base (NSRDB), \emph{Renewable and Sustainable Energy Reviews}, Volume 89,
2018, Pages 51-60, \url{https://doi.org/10.1016/j.rser.2018.03.003}.
}
\author{
D. Yang
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{get}% use one of  RShowDoc("KEYWORDS")
\keyword{PSM}% __ONLY ONE__ keyword per line
