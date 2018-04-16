#' @importFrom utils download.file

#' @export
SURFRAD.get <- function(station, year, day_of_year, directory = "data-raw")
{
  if(nchar(station) == 3)
    station <- switch(tolower(station), bon={"Bondville_IL"}, tbl={"Boulder_CO"}, dra ={"Desert_Rock_NV"}, fpk = {"Fort_Peck_MT"}, gwn = {"Goodwin_Creek_MS"}, psu = {"Penn_State_PA"}, sxf = {"Sioux_Falls_SD"})

  if(nchar(year) != 4)
    stop("Enter a four digit character string for year, e.g., '2016'")

  station.abb <- switch(station, Bondville_IL={"bon"}, Boulder_CO={"tbl"}, Desert_Rock_NV ={"dra"}, Fort_Peck_MT = {"fpk"}, Goodwin_Creek_MS = {"gwn"}, Penn_State_PA = {"psu"}, Sioux_Falls_SD = {"sxf"})

  if(!file.exists(directory))
    dir.create(directory)

  for(i in 1:length(day_of_year))
  {
    file <- paste(station.abb, substr(year, 3, 4), sprintf("%03.f", day_of_year[i]), ".dat", sep = "")
    dest <- paste(directory, file, sep = "/")
    URL <- paste("ftp://aftp.cmdl.noaa.gov/data/radiation/surfrad", station, year, file, sep = "/")
    utils::download.file(URL, destfile = dest)

    cat("SURFRAD file", file, "is written in folder '", directory, "'\n")
  }

}


