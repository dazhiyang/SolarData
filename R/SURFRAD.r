#' @importFrom utils download.file txtProgressBar setTxtProgressBar data
#' @importFrom lubridate ymd_hm month
#' @importFrom tibble as_tibble
#' @importFrom stats complete.cases


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

#' @export
SURFRAD.read <- function(files, use.original.qc = FALSE, progress.bar = TRUE, directory)
{
  setwd(directory)

  # match the desired variables
  header = c("year","jday","month","day","hour","min","dt","zen",
             "dw_solar","qc_dwsolar","uw_solar","qc_uwsolar","direct_n","qc_direct_n","diffuse","qc_diffuse",
             "dw_ir","qc_dwir","dw_casetemp","qc_dwcasetemp","dw_dometemp","qc_dwdometemp",
             "uw_ir","qc_uwir","uw_casetemp","qc_uwcasetemp","uw_dometemp","qc_uwdometemp",
             "uvb","qc_uvb","par","qc_par","netsolar","qc_netsolar","netir","qc_netir","totalnet","qc_totalnet",
             "temp","qc_temp","rh","qc_rh","windspd","qc_windspd","winddir","qc_winddir","pressure","qc_pressure")
  variables <- c("dw_solar","direct_n","diffuse","pressure")
  choice <- match(variables, header)
  choice <- sort(c(1:8, choice, choice+1))
  colClasses <- rep("NULL", length(header))
  colClasses[choice] <- "numeric"

  #get the Linke turbidity for the station
  if(length(unique(substr(files, 1, 3))) != 1)
    stop("Please process one location at a time")

  utils::data("SURFRAD.loc")
  stn <- substr(files[1], 1, 3)
  stn <- match(stn, SURFRAD.loc$stn)
  LT <- as.numeric(SURFRAD.loc[stn, 8:19])

  data_all <- vector("list", length(files))
  if(progress.bar)
    pb <- utils::txtProgressBar(min = 0, max = length(files), style = 3)
  for(i in 1:length(files))
  {
    date <- strptime(x = substr(files[i], 4, nchar(files[i])-4), format = "%y%j", tz = "UTC") # convert file name to date
    yr <- date$year+1900 # year
    mon <- lubridate::month(date) # month
    res <- ifelse(yr >=2009, 1, 3) # data resolution, 1 min if yr>2009, 3 min otherwise
    # read data
    tmp <- read.table(files[i], header = FALSE, skip = 2, colClasses = colClasses)
    names(tmp) <- header[choice]
    # convert date time
    tmp <- tmp %>%
      mutate(., Time = lubridate::ymd_hm(paste(paste(tmp$year, tmp$month, tmp$day, sep = "-"), paste(tmp$hour, tmp$min, sep = ":"), sep = " "), tz = "UTC")) %>%
      dplyr::select(-(1:7)) %>%
      tibble::as_tibble(.)
    #
    if(use.original.qc)
    {
      # new variable: sum of all irradiance qc values
      tmp <- tmp %>%
        mutate(qc_all = tmp %>% dplyr::select(starts_with("qc")) %>% rowSums(.))
      # rm the non-zero qc rows, and rm all qc columns
      tmp <- tmp %>%
        filter(tmp$qc_all == 0) %>%
        dplyr::select(-starts_with("qc"))
    }else{
      # rm all qc columns
      tmp <- tmp %>%
        dplyr::select(-starts_with("qc"))
    }

    # compute (SolarData) QC parameters, note the time shift in solpos
    solpos <- calZen(Tm = tmp$Time - res*30, lat = SURFRAD.loc$lat[stn], lon = SURFRAD.loc$lon[stn], tz = 0, LT[mon], alt = SURFRAD.loc$elev[stn])
    tmp <- tmp %>%
      filter(complete.cases(.)) %>%
      mutate(Ics = solpos$Ics) %>%
      mutate(Ioh = solpos$Ioh) %>%
      mutate(Mu0 = ifelse(tmp$zen > 90, 0, cos(radians(tmp$zen)))) %>%
      mutate(Sa = solpos$Io) #%>%
      #mutate(dw_solar = ifelse(tmp$dw_solar<0 | tmp$zen>90, 0, tmp$dw_solar)) %>%
      #mutate(direct_n = ifelse(tmp$direct_n<0 | tmp$zen>90, 0, tmp$direct_n)) %>%
      #mutate(diffuse = ifelse(tmp$diffuse<0 | tmp$zen>90, 0, tmp$diffuse))
    # calculate closure
    tmp <- tmp %>%
      mutate(sum = tmp$diffuse + tmp$Mu0*tmp$direct_n)
    # calculate Rayleigh limit
    tmp <- tmp %>%
      mutate(RL = 209.3*tmp$Mu0 - 708.38*(tmp$Mu0)^2 + 1128.7*tmp$Mu0^3 -911.2*tmp$Mu0^4 + 287.85*tmp$Mu0^5+0.046725*tmp$Mu0*tmp$pressure)

    #perform basic QC on 1-min or 3-min data
    tmp <- QC.Basic(tmp)

    # complete time stamps. Even if no missing, still left_join
    ct <- tibble::as_tibble(data.frame(Time = seq(date, date+(60*24-1)*60, by  = res*60)))
    tmp <- tmp %>% left_join(ct, ., by = "Time")

    data_all[[i]] <- tmp
    if(progress.bar)
      utils::setTxtProgressBar(pb, i)
  }#end for x in files
  if(progress.bar)
    close(pb)
  data_all <- do.call("rbind", data_all)
  data_all
}#end read function


QC.Basic <- function(df)
{
  #physical and extremely-rare limit
  phy.lim.Gh <- df$Sa*1.5*(df$Mu0)^1.2 + 100
  ext.lim.Gh <- df$Sa*1.2*(df$Mu0)^1.2 + 50
  phy.lim.Dh <- df$Sa*0.95*(df$Mu0)^1.2 + 50
  ext.lim.Dh <- df$Sa*0.75*(df$Mu0)^1.2 + 30
  phy.lim.BI <- df$Sa*df$Mu0
  ext.lim.BI <- df$Sa*0.95*(df$Mu0)^1.2 + 10

  #check physical limit and flag
  df <- df %>%
    mutate(phy.lim.G = ifelse(df$dw_solar > phy.lim.Gh | df$dw_solar < -4, 1, 0)) %>%
    mutate(phy.lim.D = ifelse(df$diffuse > phy.lim.Dh | df$diffuse < -4, 1, 0)) %>%
    mutate(phy.lim.I = ifelse(df$direct_n*df$Mu0 > phy.lim.BI | df$direct_n*df$Mu0 < -4, 1, 0))
  #check extreme-rare limit and flag
  df <- df %>%
    mutate(ext.lim.G = ifelse(df$dw_solar > ext.lim.Gh | df$dw_solar < -2, 1, 0)) %>%
    mutate(ext.lim.D = ifelse(df$diffuse > ext.lim.Dh | df$diffuse < -2, 1, 0)) %>%
    mutate(ext.lim.I = ifelse(df$direct_n*df$Mu0 > ext.lim.BI | df$direct_n*df$Mu0 < -2, 1, 0))

  #check closure
  df <- df %>%
    mutate(closr = ifelse(df$dw_solar > 50 & df$zen < 75 & abs((df$sum-df$dw_solar)/df$dw_solar)*100 > 8 , 1, 0))
  df <- df %>%
    mutate(closr = ifelse(df$dw_solar > 50 & df$zen > 75 & df$zen < 93 & abs((df$sum-df$dw_solar)/df$dw_solar)*100 > 15, 1, df$closr))

  #check diffuse ratio
  df <- df %>%
    mutate(d.ratio = ifelse(df$diffuse/df$dw_solar > 1.05 & df$dw_solar>50 & df$zen < 75, 1, 0))
  df <- df %>%
    mutate(d.ratio = ifelse(df$diffuse/df$dw_solar > 1.10 & df$dw_solar>50 & df$zen > 75, 1, df$d.ratio))

  #check climatology
  df <- df %>%
    mutate(clim1 = ifelse(df$diffuse/df$dw_solar > 0.85 & df$dw_solar/df$Ics > 0.85 & df$diffuse>50, 1, 0)) %>%
    mutate(clim2 = ifelse(df$diffuse<df$RL-1 & df$diffuse/df$dw_solar < 0.80 & df$dw_solar>50, 1, 0))

  df
}




# directory = "/Volumes/Macintosh Research/Data/surfrad/raw/bon/2015"
# setwd(directory)
# files <- dir()
#
# dat <- SURFRAD.read(files, use.original.qc = F, directory = directory)

#plot(dat$dw_solar[dat$phy.lim.G==1])
#plot(dat$diffuse[dat$phy.lim.D==1])
#plot(as.numeric(dat$direct_n[dat$phy.lim.I==1]))

# dir <- "/Volumes/Macintosh Research/Data/surfrad/Linke Turbidity"
# data("SURFRAD.loc")
# TF <- t(LTF.get(SURFRAD.loc$lon, SURFRAD.loc$lat, directory = dir))
# SURFRAD.loc <- cbind(SURFRAD.loc, TF)
# names(SURFRAD.loc)[8:19] <- paste("LTF", names(SURFRAD.loc)[8:19], sep = ".")
# setwd("/Volumes/Macintosh Research/Data/surfrad")
# save(SURFRAD.loc, file = "SURFRAD.loc.RData")


