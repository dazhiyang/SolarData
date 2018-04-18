#' @importFrom utils download.file
#' @importFrom lubridate ymd_hm
#' @importFrom tibble as_tibble


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


# SURFRAD.read <- function(files, use.original.qc = FALSE, directory)
# {
#   setwd(directory)
#
#   #get the Linke turbidity for the station
#   if(length(unique(substr(files, 1, 3))) != 1)
#     stop("Please process one location at a time")
#   data("SURFRAD.loc")
#   stn <- substr(files[1], 1, 3)
#   LT <- as.numeric(SURFRAD.loc[match(stn, SURFRAD.loc$stn), 8:19])
#   header = c("year","jday","month","day","hour","min","dt","zen",
#              "dw_solar","qc_dwsolar","uw_solar","qc_uwsolar","direct_n","qc_direct_n","diffuse","qc_diffuse",
#              "dw_ir","qc_dwir","dw_casetemp","qc_dwcasetemp","dw_dometemp","qc_dwdometemp",
#              "uw_ir","qc_uwir","uw_casetemp","qc_uwcasetemp","uw_dometemp","qc_uwdometemp",
#              "uvb","qc_uvb","par","qc_par","netsolar","qc_netsolar","netir","qc_netir","totalnet","qc_totalnet",
#              "temp","qc_temp","rh","qc_rh","windspd","qc_windspd","winddir","qc_winddir","pressure","qc_pressure")
#
#   for(x in files)
#   {
#     date <- strptime(x = substr(x, 4, nchar(x)-4), format = "%y%j", tz = "UTC") #convert file name to date
#     yr <- date$year+1900 #year
#     res <- ifelse(yr >=2009, 1, 3) # data resolution, 1 min if yr>2009, 3 min otherwise
#     #read data
#     tmp <- read.table(x, header = FALSE, skip = 2)
#     names(tmp) <- header
#     #convert date time
#     tmp <- tmp %>%
#       mutate(., Time = lubridate::ymd_hm(paste(paste(tmp$year, tmp$month, tmp$day, sep = "-"), paste(tmp$hour, tmp$min, sep = ":"), sep = " "), tz = "UTC")) %>%
#       dplyr::select(-(1:7)) %>%
#       tibble::as_tibble(.)
#
#     if(use.original.qc)
#     {
#
#     }
#
#     #complete time stamps. Even if no missing, still left_join
#     ct <- tibble::as_tibble(data.frame(Time = seq(date, date+(60*24-1)*60, by  = res*60)))
#     tmp <- tmp %>% left_join(ct, .)
#
#
#
#   }
#
#
#
#
# }
#
#
#
#
#
#
#
#
# directory = "/Volumes/Macintosh Research/Data/surfrad/raw/bon/2015"
# files <- dir()


# dir <- "/Volumes/Macintosh Research/Data/surfrad/Linke Turbidity"
# data("SURFRAD.loc")
# TF <- t(LTF.get(SURFRAD.loc$lon, SURFRAD.loc$lat, directory = dir))
# SURFRAD.loc <- cbind(SURFRAD.loc, TF)
# names(SURFRAD.loc)[8:19] <- paste("LTF", names(SURFRAD.loc)[8:19], sep = ".")
# setwd("/Volumes/Macintosh Research/Data/surfrad")
# save(SURFRAD.loc, file = "SURFRAD.loc.RData")


