#' @importFrom httr write_disk GET

#' @export
PSM.get <- function(lon, lat, api.key, attributes, name, affiliation, year, leap.year, interval, utc, reason.for.use, email, mailing.list, directory = "data-raw")
{
  if(length(lat)!=length(lon))
    stop("'lon' and 'lat' must have the same length")

  if(nchar(year)>4)
    stop("NSRDB API only allows downloading one year at a time")

  if(!file.exists(directory))
    dir.create(directory)

  for(i in 1:length(lat))
  {
    # Declare url string
    URL <- paste0('http://developer.nrel.gov/api/solar/nsrdb_psm3_download.csv?wkt=POINT(', lon[i], '+', lat[i], ')&names=', year, '&leap_day=', leap.year, '&interval=', interval, '&utc=', utc, '&full_name=', name, '&email=', email, '&affiliation=', affiliation, '&mailing_list=', mailing.list, '&reason=', reason.for.use, '&api_key=', api.key, '&attributes=', attributes)
    # name the output file
    output_file <- paste0(lat[i], "_", lon[i], "_", year, ".csv")
    # API request and saving
    httr::GET(url = URL, httr::write_disk(paste(directory, output_file, sep = "/")))
    cat("PSM3 file", output_file, "is written in folder '", directory, "'\n")
  }

}



# pack <- "SolarData"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")), "CMD", "Rd2pdf", shQuote(path)))
