#' @import insol
#' @import dplyr
#' @importFrom lubridate tz
#' @importFrom geosphere destPoint alongTrackDistance
#' @importFrom fields rdist.earth
#' @importFrom utils read.table
#' @importFrom grDevices hcl

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "Time", "SURFRAD.loc"))

#################################################################################
# function to calculate sun position and extraterrestial irradiance
#################################################################################
calZen <- function(Tm, lat, lon, tz = 0, LT, alt = 0)
{
  jd = JD(Tm)
  sunv = sunvector(jd, lat, lon, tz)
  azi = round(sunpos(sunv)[,1],3)#azimuth of the sun
  zen = round(sunpos(sunv)[,2],3)#zenith angle
  #surface.norm = normalvector(tilt, orientation)
  #inc = round(as.numeric(degrees(acos(sunv%*% as.vector(surface.norm)))),3)
  dec = declination(jd)*pi/180
  re = 1.000110+0.034221*cos(dec)+0.001280*sin(dec)+0.00719*cos(2*dec)+0.000077*sin(2*dec)
  Io = round(1362*re,3)#extraterrestrial direct normal irradiance
  Ioh = round(1362*re*cos(radians(zen)))#horizontal extraterrestrial irradiance
  Ioh <- ifelse(zen>=90, 0, Ioh)

  # Equation of time (L. O. Lamm, 1981, Solar Energy 26, p465)
  dn <- round(as.numeric(format(Tm, "%j"))-1 + as.numeric(format(Tm, "%H"))/24, 3)
  coef <- matrix(c(0, 0.00020870, 0,
                   1, 0.0092869, -0.12229,
                   2, -0.052258, -0.15698,
                   3, -0.0013077, -0.0051602,
                   4, -0.0021867, -0.0029823,
                   5, -0.00015100, -0.00023463), ncol = 3, byrow = TRUE)
  EOT <- rowSums(sapply(1:6, function(i) coef[i,2]*cos(2*pi*coef[i,1]*dn/365.25) + coef[i,3]*sin(2*pi*coef[i,1]*dn/365.25)))*60 #EOT in minutes
  Tsolar <- Tm - 4*60*(tz*15-lon) + EOT*60

  #find pressure from altitude
  pressure = 100 * ((44331.514 - alt)/11880.516)^(1/0.1902632)

  #Perez-Ineichen clear sky model (Ineichen and Perez, 2002), monthly Linke turbidity can be obtained from SoDa service (following: Gueymard and Ruiz-Aries, 2015)
  fh1 <- exp(-alt/8000)
  fh2 <- exp(-alt/1250)
  cg1 <- (0.0000509*alt + 0.868)
  cg2 <- (0.0000392*alt + 0.0387)
  z <- pmin(zen, 90)
  #M <- 1/(cos(radians(z))+0.50572*(96.07995 - z)^(-1.6364))
  AM <- 1/(cos(radians(z))+0.00176759*(z)*((94.37515 - z)^(-1.21563)))
  AM <- AM/101325*pressure #elevation corrected AM
  Ics <- cg1*Io*cos(radians(z))*exp(-cg2*AM*(fh1 + fh2*(LT-1)))*exp(0.01*pmin(AM,12)^1.8)
  Ics <- ifelse(zen>=90, 0, Ics)
  Icsd <- (0.664+0.163/fh1)*Io*exp(-0.09*(LT-1)*AM)*cos(z)
  Icsd <- ifelse(zen>=90, 0, Icsd)

  out = list(zen, Io, Ioh, Ics, Icsd, Tsolar)
  names(out) = c("zenith", "Io", "Ioh", "Ics", "Icsd", "Tsolar")
  out
}

################################################################################
#function to calculate the geographical distance amongst the stations
################################################################################
GetGeoDist <- function(lon, lat)
{
  dist = matrix(NA, length(lat), length(lat))
  for(i in 1:length(lat)){
    for(j in 1:length(lat)){
      dist[i,j] = fields::rdist.earth(matrix(c(lon[i],lat[i]),ncol =2),matrix(c(lon[j],lat[j]),ncol =2), miles = FALSE)
    }
  }
  dist = as.matrix(dist)
}

GetAlongWindDist <- function(lon, lat, wind.dir)
{
  x = cbind(lon, lat)
  dist = matrix(NA, nrow(x), nrow(x))
  for(i in 1:nrow(x))
  {
    end.point = geosphere::destPoint(x[i,], wind.dir, 2000)
    temp = geosphere::alongTrackDistance(x[i,], end.point, x)/1000
    dist[i, ] = ifelse(is.na(temp), 0, temp)
  }
  dist
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

################################################################################
#function to calculate the ceiling of time, lubridate doesn't work for 50 sec!!
################################################################################
ceiling.time <- function(Tm, agg.interval)
{
  if(class(Tm)[1] != "POSIXct")
    stop("'Tm' must be a POSIXct object")

  TZ <- lubridate::tz(Tm[1])
  #convert to UTC
  attributes(Tm)$tzone <- "UTC"
  tmp <- ceiling(unclass(Tm)/(agg.interval))*(agg.interval)
  origin <- as.POSIXct('1970-01-01 00:00:00', tz = "UTC")
  Tm <- as.POSIXct(tmp, origin = origin, tz = "UTC")
  attributes(Tm)$tzone <- TZ
  Tm
}


################################################################################
#function to read and aggregate OSMG
################################################################################
#' @export
OSMG.read <- function(files, directory.LI200, directory.RSR = NULL, clear.sky = FALSE, AP2 = FALSE, agg = 1)
{
  #LI-200 station names
  stations = c("DH3", "DH4", "DH5", "DH10", "DH11", "DH9", "DH2", "DH1", "DH1T", "AP6", "AP6T", "AP1", "AP3", "AP5", "AP4", "AP7", "DH6", "DH7", "DH8")

  setwd(directory.LI200) #read LI-200 files first
  data_all <- NULL
  for(x in files)
  {
    cat("Reading and processing", x, "...\n")
    setwd(directory.LI200)
    #read data
    data <- read.table(x, header = FALSE, sep = ",", colClasses = c(rep("character", 4), rep("numeric", 19))) #read data
    names(data) <- c("SS", "year", "doy", "HM", stations)

    #arrange the date and time
    data <- data %>%
      mutate(., SS = ifelse(nchar(data$SS)==2, data$SS, paste0("0", data$SS))) %>%
      mutate(., HM = ifelse(nchar(data$HM)==4, data$HM, paste0("0", data$HM)))
    Tm <- as.POSIXct(paste(substr(x, 1, 8), data$HM, data$SS, sep = "-"), format = "%Y%m%d-%H%M-%S", tz = "HST") #Time format

    #select the horizontal stations and round to 3 decimal
    data <- data %>%
      dplyr::select(., c(5:23)) %>%
      mutate_all(., funs(round(., 3))) %>%
      mutate_all(., funs(replace(., .<0, 0)))

    if(clear.sky)
    {
      #since loading the tiff images is slow, the values are listed here
      LT <- c(3.50, 3.80, 4.05, 4.55, 4.45, 4.55, 4.50, 4.35, 4.45, 4.45, 3.85, 3.75)
      month <- as.numeric(substr(x, 5, 6))
      solpos <- calZen(Tm = Tm, lat = 21.31234, lon = -158.0841, LT = LT[month], alt = 3.9)
      data <- data %>%
        mutate(., zen = solpos$zenith, Ics = solpos$Ics, Ioh = solpos$Ioh) %>%
        dplyr::select(., c(20:22, 9, 11, 1:8, 10, 12:19))
    }

    if(AP2 & is.null(directory.RSR))
    {
      stop("Please specify the directory for RSR (AP2) dataset")
    }

    if(agg < 3 & AP2){
      stop("AP2 is 3-sec data, please choose a higher 'agg' value")
    }else if(agg >= 3 & AP2)
    {
      setwd(directory.RSR)
      data_AP2 = utils::read.table(x, header = FALSE, sep = ",", colClasses = c(rep("character", 4), rep("numeric", 3))) #read data
      names(data_AP2) <- c("SS", "year", "doy", "HM", "AP2", "AP2.dif", "AP2.dir")
      #arrange the date and time
      data_AP2 <- data_AP2 %>%
        mutate(., SS = ifelse(nchar(data_AP2$SS)==2, data_AP2$SS, paste0("0", data_AP2$SS))) %>%
        mutate(., HM = ifelse(nchar(data_AP2$HM)==4, data_AP2$HM, paste0("0", data_AP2$HM)))
      Tm_AP2 <- as.POSIXct(paste(substr(x, 1, 8), data_AP2$HM, data_AP2$SS, sep = "-"), format = "%Y%m%d-%H%M-%S", tz = "HST") #Time format
      #select GHI, DIF, and DIR
      data_AP2 <- data_AP2 %>%
        dplyr::select(., c(5:7)) %>%
        mutate_all(., funs(round(., 3))) %>%
        mutate_all(., funs(replace(., .<0, 0))) %>%
        mutate(., Time = ceiling.time(Tm_AP2, agg)) %>%
        group_by(Time) %>%
        summarise_all(funs(mean), na.rm = TRUE)

      #aggregate 1-sec data
      data <- data %>%
        mutate(., Time = ceiling.time(Tm, agg)) %>%
        group_by(Time) %>%
        summarise_all(funs(mean), na.rm = TRUE)

      data_day <- left_join(data, data_AP2)
    }else if(!AP2){
      #aggregate 1-sec data
      data <- data %>%
        mutate(., Time = ceiling.time(Tm, agg)) %>%
        group_by(Time) %>%
        summarise_all(funs(mean), na.rm = TRUE)

      data_day <- data
    } #end joining if's

    data_all <- rbind(data_all, data_day)

  } #end for loop for files

  data_all
}


