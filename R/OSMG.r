#' @import insol
#' @import dplyr
#' @importFrom lubridate ceiling_date
#' @importFrom geosphere destPoint alongTrackDistance
#' @importFrom fields rdist.earth
#' @importFrom utils read.table

. <- "Shut up"
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
  Ioh = round(1362*re*cos(d2r(zen)))#horizontal extraterrestrial irradiance
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

  #Perez-Ineichen clear sky model (Ineichen and Perez, 2002), monthly Linke turbidity can be obtained from SoDa service (following: Gueymard and Ruiz-Aries, 2015)
  fh1 <- exp(-alt/8000)
  fh2 <- exp(-alt/1250)
  cg1 <- (0.0000509*alt + 0.868)
  cg2 <- (0.0000392*alt + 0.0387)
  AM <- 1/(cos(zen*pi/180)+0.50572*(96.07995 - zen)^(-1.6364))
  Ics <- cg1*Io*cos(zen*pi/180)*exp(-cg2*AM*(fh1 + fh2*(LT-1)))*exp(0.01*AM^1.8)
  Ics <- ifelse(zen>=90, 0, Ics)
  Icsd <- (0.664+0.163/fh1)*Io*exp(-0.09*(LT-1)*AM)*cos(zen*pi/180)
  Icsd <- ifelse(zen>=90, 0, Icsd)

  out = list(zen, Io, Ioh, Ics, Icsd, Tsolar)
  names(out) = c("zenith", "Io", "Ioh", "Ics", "Icsd", "Tsolar")
  out
}

#degree to radian
d2r <-function(x)
{
  x*pi/180
}
#radian to degree
r2d <-function(x)
{
  x*180/pi
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

################################################################################
#function to read and aggregate OSMG
################################################################################
#' @export
OSMG.read <- function(files, directory_LI200, directory_RSR = NULL, clear_sky = FALSE, AP2 = FALSE, agg = 1)
{
  #LI-200 station names
  stations = c("DH3", "DH4", "DH5", "DH10", "DH11", "DH9", "DH2", "DH1", "DH1T", "AP6", "AP6T", "AP1", "AP3", "AP5", "AP4", "AP7", "DH6", "DH7", "DH8")

  setwd(directory_LI200) #read LI-200 files first
  data_all <- NULL
  for(x in files)
  {
    cat("Reading and processing", x, "...\n")
    setwd(directory_LI200)
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

    if(clear_sky)
    {
      #since loading the tiff images is slow, the values are listed here
      LT <- c(3.50, 3.80, 4.05, 4.55, 4.45, 4.55, 4.50, 4.35, 4.45, 4.45, 3.85, 3.75)
      month <- as.numeric(substr(x, 5, 6))
      solpos <- calZen(Tm = Tm, lat = 21.31234, lon = -158.0841, LT = LT[month], alt = 3.9)
      data <- data %>%
        mutate(., zen = solpos$zenith, Ics = solpos$Ics, Ioh = solpos$Ioh) %>%
        dplyr::select(., c(20:22, 9, 11, 1:8, 10, 12:19))
    }

    if(AP2 & is.null(directory_RSR))
    {
      stop("Please specify the directory for RSR (AP2) dataset")
    }

    if(agg < 3 & AP2){
      stop("AP2 is 3-sec data, please choose a higher 'agg' value")
    }else if(agg >= 3 & AP2)
    {
      setwd(directory_RSR)
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
        mutate(., Time = lubridate::ceiling_date(Tm_AP2, paste0(agg, " seconds"))) %>%
        group_by(Time) %>%
        summarise_all(funs(mean), na.rm = TRUE)

      #aggregate 1-sec data
      data <- data %>%
        mutate(., Time = lubridate::ceiling_date(Tm, paste0(agg, " seconds"))) %>%
        group_by(Time) %>%
        summarise_all(funs(mean), na.rm = TRUE)

      data_day <- left_join(data, data_AP2)
    }else if(!AP2){
      #aggregate 1-sec data
      data <- data %>%
        mutate(., Time = lubridate::ceiling_date(Tm, paste0(agg, " seconds"))) %>%
        group_by(Time) %>%
        summarise_all(funs(mean), na.rm = TRUE)

      data_day <- data
    } #end joining if's

    data_all <- rbind(data_all, data_day)

  } #end for loop for files

  data_all
}







# loc = read.csv("/Users/DYang/Dropbox/Working papers/HEM/Data/HawaiiStations.csv", header = TRUE)
# loc = loc[-c(9, 11),] #remove tilted stations
# lat = loc$Latitude; lon = loc$Longitude; loc.name = loc$Name
# dist.aw = GetAlongWindDist(lon, lat, wind.dir = 60)
# aw.order = order(dist.aw[14,]) #14 is AP7, order in increasing distance
# loc = loc[aw.order,]#order the locations in along-wind direction
# lat = loc$Latitude; lon = loc$Longitude;
# dist <- GetGeoDist(lon, lat)
# dist.cw = GetAlongWindDist(lon, lat, wind.dir = -30) #get cross wind distance
# dist.aw = GetAlongWindDist(lon, lat, wind.dir = 60)
# #get distance in meters
# dist.aw = dist.aw*1000; dist.cw = dist.cw*1000; dist = dist*1000;
#
# #Define along-wind and cross-wind pairs
# #AP7 = 1, AP4 = 2, AP3 = 3, AP2 = 4, AP6 = 5, DH5 = 6, AP1 = 7, DH2 = 8, AP5 = 9, DH3 = 10, DH4 = 11, DH1 = 12, DH7 = 13, DH10 = 14, DH11 = 15, DH9 = 16, DH6 = 17, DH8 = 18
# pair.aw = matrix(c( 7,  6, 10, 4,  7, 11,  4,  6, 2, 1,  3,  1,  1,
#                     10, 11, 14, 7, 14, 17, 10, 17, 9, 3, 15, 15, 18), ncol = 2, byrow = FALSE)
# pair.cw = matrix(c(11, 6, 13, 7, 6, 8, 4,
#                    10, 7, 14, 9, 9, 9, 5), ncol = 2, byrow = FALSE)
#
#
#
# directory_LI200 <- "/Volumes/Macintosh Research/Data/Oahu/raw"
# directory_RSR <- "/Volumes/Macintosh Research/Data/Oahu/raw AP2"
# setwd(directory_LI200)
# files <- dir()
#
# data <- OSMG.read(files, directory_LI200 = "/Volumes/Macintosh Research/Data/Oahu/raw", directory_RSR = "/Volumes/Macintosh Research/Data/Oahu/raw AP2", clear_sky = TRUE, AP2 = TRUE, agg = 60)
# data <- data.frame(data)
#
# Tm = as.POSIXlt(data[,1], format = "%Y-%m-%d %H:%M:%S", tz = "HST")
# select = which(Tm$hour>8&Tm$hour<15) #select 09:00 to 15:00
# data = data[select, 7:24] #remove time column from the data frame, as well as those "non-noon" time
# data = apply(data, 2, function(x) diff(x))
# data = data[,aw.order] #order the data in along-wind direction
# C = cor(data, use = "pairwise.complete.obs")
#
# mat.aw = mat.cw = mat.other.cor <- matrix(FALSE, 18, 18)
# mat.other.cor[lower.tri(mat.other.cor)] <- TRUE
# for(i in 1:nrow(pair.aw)){
#   mat.aw[pair.aw[i, 2], pair.aw[i, 1]] <- TRUE
#   mat.other.cor[pair.aw[i, 2], pair.aw[i, 1]] <- FALSE
# }
# C.aw <- C[mat.aw]; d.aw <- dist[mat.aw]; a.aw <- dist.aw[mat.aw]; b.aw <- dist.cw[mat.aw];
#
# for(i in 1:nrow(pair.cw)){
#   mat.cw[pair.cw[i, 2], pair.cw[i, 1]] <- TRUE
#   mat.other.cor[pair.aw[i, 2], pair.aw[i, 1]] <- FALSE
# }
# C.cw <- C[mat.cw]; d.cw <- dist[mat.cw]; a.cw <- dist.cw[mat.cw]; b.cw <- dist.cw[mat.cw];
#
# data.plot.point.aw <- data.frame(x = d.aw, y = C.aw)
#
#
#
# plot(data.plot.point.aw[,1], data.plot.point.aw[,2])
#
#
#
# data.plot.point.aw = rbind(data.plot.point.aw, data.frame(x = d.aw, y = C.aw, timescale = paste("Avg-", time.scale[j], "s", sep = "")))
# data.plot.point.cw = rbind(data.plot.point.cw, data.frame(x = d.cw, y = C.cw, timescale = paste("Avg-", time.scale[j], "s", sep = "")))
# data.plot.point.neither = rbind(data.plot.point.neither, data.frame(x = dist[mat.other.cor], y = C[mat.other.cor], timescale = paste("Avg-", time.scale[j], "s", sep = "")))








