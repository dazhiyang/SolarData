#' @importFrom httr write_disk GET

#' @export
PSM.get <- function(lon, lat, api.key, attributes, name, affiliation, year, leap.year, interval, utc, reason.for.use, email, mailing.list, directory = "data-raw")
{
  if(length(lat)!=length(lon))
    stop("'lon' and 'lat' must have the same length")

  year <- as.character(year)
  if(nchar(year)>4)
    stop("NSRDB API only allows downloading one year at a time")

  interval <- as.character(interval)
  if(interval %in% c("5", "15") & !(year %in% c("2018", "2019")))
    stop("For 5-min data, only years 2018 and 2019 are available at the moment")

  if(!file.exists(directory))
    dir.create(directory)

  for(i in 1:length(lat))
  {
    # Declare url string
    if(interval %in% c("5", "15"))
    {
      URL <- paste0('https://developer.nrel.gov/api/nsrdb/v2/solar/psm3-5min-download.csv?wkt=POINT(', lon[i], '+', lat[i], ')&names=', year, '&leap_day=', leap.year, '&interval=', interval, '&utc=', utc, '&full_name=', name, '&email=', email, '&affiliation=', affiliation, '&mailing_list=', mailing.list, '&reason=', reason.for.use, '&api_key=', api.key, '&attributes=', attributes)
    }else{
      URL <- paste0('https://developer.nrel.gov/api/nsrdb/v2/solar/psm3-download.csv?wkt=POINT(', lon[i], '+', lat[i], ')&names=', year, '&leap_day=', leap.year, '&interval=', interval, '&utc=', utc, '&full_name=', name, '&email=', email, '&affiliation=', affiliation, '&mailing_list=', mailing.list, '&reason=', reason.for.use, '&api_key=', api.key, '&attributes=', attributes)
    }

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
