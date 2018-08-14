#' @importFrom utils download.file txtProgressBar setTxtProgressBar data
#' @importFrom lubridate ymd_hm month
#' @importFrom tibble as_tibble
#' @importFrom stats complete.cases


#' @export
SURFRAD.get <- function(station, year, day.of.year, directory = "data-raw")
{
  if(nchar(station) == 3)
    station <- switch(tolower(station), bon={"Bondville_IL"}, tbl={"Boulder_CO"}, dra ={"Desert_Rock_NV"}, fpk = {"Fort_Peck_MT"}, gwn = {"Goodwin_Creek_MS"}, psu = {"Penn_State_PA"}, sxf = {"Sioux_Falls_SD"})

  if(nchar(year) != 4)
    stop("Enter a four digit character string for year, e.g., '2016'")

  station.abb <- switch(station, Bondville_IL={"bon"}, Boulder_CO={"tbl"}, Desert_Rock_NV ={"dra"}, Fort_Peck_MT = {"fpk"}, Goodwin_Creek_MS = {"gwn"}, Penn_State_PA = {"psu"}, Sioux_Falls_SD = {"sxf"})

  if(!file.exists(directory))
    dir.create(directory)

  for(i in 1:length(day.of.year))
  {
    file <- paste(station.abb, substr(year, 3, 4), sprintf("%03.f", day.of.year[i]), ".dat", sep = "")
    dest <- paste(directory, file, sep = "/")
    URL <- paste("ftp://aftp.cmdl.noaa.gov/data/radiation/surfrad", station, year, file, sep = "/")
    utils::download.file(URL, destfile = dest)

    cat("SURFRAD file", file, "is written in folder '", directory, "'\n")
  }

}

# directory <- "/Volumes/Macintosh Research/Data/surfrad/raw/bon/2005"
# setwd(directory)
# files <- dir()
# use.original.qc = FALSE
# use.qc = TRUE
# test = c("ext", "closr")
# progress.bar = TRUE


#' @export
SURFRAD.read <- function(files, directory, use.original.qc = FALSE, use.qc = TRUE, test = NULL, progress.bar = TRUE, agg = 1, additional.variables = NULL)
{
  setwd(directory)

  # match the desired variables
  header = c("year","jday","month","day","hour","min","dt","zen",
             "dw_solar","qc_dwsolar","uw_solar","qc_uwsolar","direct_n","qc_direct_n","diffuse","qc_diffuse",
             "dw_ir","qc_dwir","dw_casetemp","qc_dwcasetemp","dw_dometemp","qc_dwdometemp",
             "uw_ir","qc_uwir","uw_casetemp","qc_uwcasetemp","uw_dometemp","qc_uwdometemp",
             "uvb","qc_uvb","par","qc_par","netsolar","qc_netsolar","netir","qc_netir","totalnet","qc_totalnet",
             "temp","qc_temp","rh","qc_rh","windspd","qc_windspd","winddir","qc_winddir","pressure","qc_pressure")

  variables <- c("dw_solar","direct_n","diffuse","pressure", additional.variables)
  choice <- match(variables, header)
   if(length(which(is.na(choice)))!=0)
    stop("Names in 'additional.variable' are defined incorrectly")
  choice <- sort(c(1:8, choice, choice+1))
  colClasses <- rep("NULL", length(header))
  colClasses[choice] <- "numeric"

  if(length(unique(substr(files, 1, 3))) != 1)
    stop("Please process one location at a time")

  if(length(unique(substr(files, 4, 5))) != 1)
    stop("Please process one year at a time")

  #get the Linke turbidity for the station
  utils::data("SURFRAD.loc")
  stn <- substr(files[1], 1, 3)
  stn <- match(stn, SURFRAD.loc$stn)
  LT <- as.numeric(SURFRAD.loc[stn, 8:19])

  data_all <- vector("list", length(files))
  if(progress.bar)
    pb <- utils::txtProgressBar(min = 0, max = length(files), style = 3)
  for(i in 1:length(files))
  {
    #get year and raw data resolution
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

    # complete time stamps. Even if no missing, still left_join
    ct <- tibble::as_tibble(data.frame(Time = seq(date, date+(60*24-1)*60, by  = res*60)))
    tmp <- tmp %>% left_join(ct, ., by = "Time")

    # compute (SolarData) QC parameters, note the time shift in solpos
    solpos <- calZen(Tm = tmp$Time - res*30, lat = SURFRAD.loc$lat[stn], lon = SURFRAD.loc$lon[stn], tz = 0, LT[mon], alt = SURFRAD.loc$elev[stn])
    # fill zenith angle for missing data points
    tmp <- tmp %>%
      mutate(zen = ifelse(is.na(tmp$zen), solpos$zenith, tmp$zen))
    # fill other solar positioning parameters
    tmp <- tmp %>%
      mutate(Ics = solpos$Ics) %>%
      mutate(Ioh = solpos$Ioh) %>%
      mutate(Mu0 = ifelse(tmp$zen > 90, 0, cos(radians(tmp$zen)))) %>%
      mutate(Sa = solpos$Io)
    # calculate closure
    tmp <- tmp %>%
      mutate(sum = tmp$diffuse + tmp$Mu0*tmp$direct_n)
    # calculate Rayleigh limit
    tmp <- tmp %>%
      mutate(RL = 209.3*tmp$Mu0 - 708.38*(tmp$Mu0)^2 + 1128.7*tmp$Mu0^3 -911.2*tmp$Mu0^4 + 287.85*tmp$Mu0^5+0.046725*tmp$Mu0*ifelse(tmp$pressure>0 & !is.na(tmp$pressure), tmp$pressure, mean(tmp$pressure, na.rm = TRUE))) #if pressure is missing or negative, replace with the mean value

    #perform basic QC on 1-min or 3-min data
    if(use.qc)
    {
      if(is.null(test))
        stop("You must specify the test(s), avialable ones include 'phy', 'ext', 'closr', 'dr', 'clim', 'all'.")
      tmp <- QC.Basic(tmp, test)
      # new variable: sum of all irradiance qc values
      tmp <- tmp %>%
        mutate(qc_all = tmp %>% dplyr::select(starts_with("qc")) %>% rowSums(.))
      # set NA at the non-zero qc rows, and rm all qc columns
      tmp <- tmp %>%
        mutate_at(.vars = c(3:5), funs(ifelse(tmp$qc_all>0, NA, .))) %>%
        dplyr::select(-starts_with("qc"))
    }

    data_all[[i]] <- tmp
    if(progress.bar)
      utils::setTxtProgressBar(pb, i)
  }#end for x in files
  if(progress.bar)
    close(pb)

  cat("Concatenating files ...\n")
  data_all <- do.call("rbind", data_all)
  data_all <- data_all %>%
    mutate(dw_solar = ifelse(data_all$dw_solar<0 | data_all$zen>90, 0, data_all$dw_solar)) %>%
    mutate(direct_n = ifelse(data_all$direct_n<0 | data_all$zen>90, 0, data_all$direct_n)) %>%
    mutate(diffuse = ifelse(data_all$diffuse<0 | data_all$zen>90, 0, data_all$diffuse)) %>%
    dplyr::select(-matches("Mu0|Sa|sum|RL"))

  #time difference
  diffs <- as.numeric(data_all$Time[2:nrow(data_all)]-data_all$Time[1:(nrow(data_all)-1)])

  #get serially complete data
  start <- data_all$Time[1]
  end <- data_all$Time[nrow(data_all)]
  Time_all <- data_frame(Time = seq(start, end, by = 60*res))
  if(nrow(Time_all) != nrow(data_all))
  {
    data_all <- data_all %>%
      right_join(., Time_all, by = "Time")
    missing <- which(is.na(data_all$zen))
    solpos <- calZen(Tm = data_all$Time[missing] - res*30, lat = SURFRAD.loc$lat[stn], lon = SURFRAD.loc$lon[stn], tz = 0, LT[mon], alt = SURFRAD.loc$elev[stn])
    data_all$Ics[missing] <- solpos$Ics
    data_all$Ioh[missing] <- solpos$Ioh
    data_all$zen[missing] <- solpos$zen
  }

  #aggregate
  if(min(diffs) == 3 & agg < 3)
  {
    stop("you have 3-min data, 'agg' must be at least 3\n")
  }else if(min(diffs) == 3 & agg == 3){
    cat("you have 3-min data, no aggregation required\n")
  }else if(min(diffs) == 1 & agg == 1){
    cat("you have 1-min data, no aggregation required\n")
  }else if(agg>res){
    cat("Aggregating files ...\n")
    data_all <- data_all %>%
      mutate(., Time = ceiling.time(data_all$Time, agg*60))
    #check number of points in each aggregation interval
    n.point.in.each.interval <- array(0, length(unique(data_all$Time)))
    if(length(which(is.na(data_all$dw_solar)))!=0)
    {
      non.empty.interval <- match(unique(data_all$Time[-which(is.na(data_all$dw_solar))]), unique(data_all$Time))
      n.point.in.each.interval[non.empty.interval] <- rle(as.numeric(data_all$Time[-which(is.na(data_all$dw_solar))]))$length
      #assign NAs, so that when aggregate, these intervals will be NA-valued.
      bad.interval <- unique(data_all$Time)[which(n.point.in.each.interval < floor(agg/res/2))][-1]
      remove <- which(data_all$Time %in% bad.interval)
      data_all[remove, c(3:5)] <- NA
    }
    data_all <- data_all %>%
      group_by(Time) %>%
      summarise_all(funs(mean), args = list(na.rm = TRUE))
  }

  # output
  data_all
}#end read function


QC.Basic <- function(df, test)
{
  #check physical limit and flag
  if("phy" %in% test | "all" %in% test)
  {
    df <- df %>%
      mutate(qc_phy_G = ifelse(df$dw_solar > df$Sa*1.5*(df$Mu0)^1.2 + 100 | df$dw_solar < -4, 1, 0)) %>%
      mutate(qc_phy_D = ifelse(df$diffuse > df$Sa*0.95*(df$Mu0)^1.2 + 50 | df$diffuse < -4, 1, 0)) %>%
      mutate(qc_phy_I = ifelse(df$direct_n*df$Mu0 > df$Sa*df$Mu0 | df$direct_n*df$Mu0 < -4, 1, 0))
  }

  #check extreme-rare limit and flag
  if("ext" %in% test | "all" %in% test)
  {
    df <- df %>%
      mutate(qc_ext_G = ifelse(df$dw_solar > df$Sa*1.2*(df$Mu0)^1.2 + 50 | df$dw_solar < -2, 1, 0)) %>%
      mutate(qc_ext_D = ifelse(df$diffuse > df$Sa*0.75*(df$Mu0)^1.2 + 30 | df$diffuse < -2, 1, 0)) %>%
      mutate(qc_ext_I = ifelse(df$direct_n*df$Mu0 > df$Sa*0.95*(df$Mu0)^1.2 + 10 | df$direct_n*df$Mu0 < -2, 1, 0))
  }

  #check closure
  if("closr" %in% test | "all" %in% test)
  {
    df <- df %>%
      mutate(qc_closr = ifelse(df$dw_solar > 50 & df$zen < 75 & abs((df$sum-df$dw_solar)/df$dw_solar)*100 > 8 , 1, 0))
    df <- df %>%
      mutate(qc_closr = ifelse(df$dw_solar > 50 & df$zen > 75 & df$zen < 93 & abs((df$sum-df$dw_solar)/df$dw_solar)*100 > 15, 1, df$qc_closr))
  }

  #check diffuse ratio
  if("dr" %in% test | "all" %in% test)
  {
    df <- df %>%
      mutate(qc_d_ratio = ifelse(df$diffuse/df$dw_solar > 1.05 & df$dw_solar>50 & df$zen < 75, 1, 0))
    df <- df %>%
      mutate(qc_d_ratio = ifelse(df$diffuse/df$dw_solar > 1.10 & df$dw_solar>50 & df$zen > 75, 1, df$qc_d_ratio))
  }

  #check climatology
  if("clim" %in% test | "all" %in% test)
  {
    df <- df %>%
      mutate(qc_clim1 = ifelse(df$diffuse/df$dw_solar > 0.85 & df$dw_solar/df$Ics > 0.85 & df$diffuse>50, 1, 0)) %>%
      mutate(qc_clim2 = ifelse(df$diffuse<df$RL-1 & df$diffuse/df$dw_solar < 0.80 & df$dw_solar>50, 1, 0))
  }

  df
}



# dir <- "/Volumes/Macintosh Research/Data/surfrad/Linke Turbidity"
# data("SURFRAD.loc")
# TF <- t(LTF.get(SURFRAD.loc$lon, SURFRAD.loc$lat, directory = dir))
# SURFRAD.loc <- cbind(SURFRAD.loc, TF)
# names(SURFRAD.loc)[8:19] <- paste("LTF", names(SURFRAD.loc)[8:19], sep = ".")
# setwd("/Volumes/Macintosh Research/Data/surfrad")
# save(SURFRAD.loc, file = "SURFRAD.loc.RData")


