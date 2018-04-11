#' @importFrom httr write_disk GET

#' @export
PSM.get <- function(lon, lat, api_key, attributes, name, affiliation, year, leap_year, interval, utc, reason_for_use, email, mailing_list, directory = "data-raw")
{
  if(length(lat)!=length(lon))
    stop("'lat' and 'lon' must have the same length")

  if(nchar(year)>4)
    stop("NSRDB API only allows downloading one year at a time")

  if(!file.exists(directory))
    dir.create(directory)

  for(i in 1:length(lat))
  {
    # Declare url string
    URL <- paste0('http://developer.nrel.gov/api/solar/nsrdb_psm3_download.csv?wkt=POINT(', lon[i], '+', lat[i], ')&names=', year, '&leap_day=', leap_year, '&interval=', interval, '&utc=', utc, '&full_name=', name, '&email=', email, '&affiliation=', affiliation, '&mailing_list=', mailing_list, '&reason=', reason_for_use, '&api_key=', api_key, '&attributes=', attributes)
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
