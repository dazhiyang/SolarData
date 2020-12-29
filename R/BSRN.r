#' @importFrom utils data
#' @importFrom lubridate ymd_hm month
#' @importFrom tibble as_tibble tibble
#' @importFrom stats complete.cases
#' @importFrom RCurl getURL getCurlHandle getBinaryURL


#' @export
BSRN.list <- function(station, user, pwd)
{
  if(is.null(user) | is.null(pwd))
    stop("You must have BSRN account and password, ask Amelie for it.")
  userpwd <- paste0(user, ":", pwd)

  if(length(station) == 0)
    stop("Please enter the station name(s) as string or vector of strings")
  stn<- tolower(station)

  BSRNdir <- list()
  for(i in 1:length(stn))
  {
    url <- paste0("ftp://ftp.bsrn.awi.de/", stn[i], "/")
    cat("Retrieving station", stn[i], "...\n")
    filenames <- RCurl::getURL(url, userpwd=userpwd, ftp.use.epsv = FALSE, dirlistonly = FALSE)
    destnames <- strsplit(filenames, "\n")[[1]]
    destnames <- read.table(text=destnames)
    BSRNdir[[i]] <- tibble(filename = destnames$V9, filesize = destnames$V5, last.update = as.POSIXct(paste(destnames$V8, destnames$V6, destnames$V7, sep = "-"), format = "%Y-%b-%d", tz = "UTC") )
  }
  names(BSRNdir) <- stn

  BSRNdir
}

#' @export
BSRN.get <- function(station, mmyy, directory, user, pwd)
{
  if(is.null(directory))
    directory <- getwd()

  # build a list of URLs for requested files
  files <- sapply(mmyy, function(x) paste0("ftp://ftp.bsrn.awi.de/", station, "/", station, x, ".dat.gz"), simplify = F)

  # build a list of paths to save gotten files.
  wdfiles <- sapply(mmyy, function(x) paste0(directory, "/", station, x, ".dat.gz"), simplify = F)

  # request files using RCurl
  con <-  RCurl::getCurlHandle(ftp.use.epsv = FALSE, userpwd= paste0(user, ":", pwd))
  invisible(mapply(function(x,y) writeBin(RCurl::getBinaryURL(x, curl = con), y), x = files, y = wdfiles))
}



#' @export
BSRN.read <- function(file, directory, use.qc = TRUE, test = NULL, use.agg = FALSE, agg = 1)
{
  setwd(directory)

  if(length(file) != 1)
    stop("Please process one file at a time")

  if(!grepl("*.dat|*.gz", file))
    stop("File extension is not recognized, use .dat or. gz")

  #read data
  if(grepl("*.gz", file))
  {
    zz <- gzfile(file, open = "rt")
    lines <- readLines(zz, encoding = "ACSII")
    close(zz)
  }else{
    lines <- readLines(file)
  }

  #define start and end of *0100 records
  LR <- which(grepl("\\*", lines)) #all lines indicating the start of logical record
  LR.pos <- which(grepl("\\*U0100|\\*C0100", lines[LR])) #position of '0100 record' in LR
  start <-  LR[LR.pos] + 1 #first row of *0100 data
  end <- ifelse(LR.pos == length(LR), length(lines), LR[LR.pos+1]-1)

  #read *0100 data
  line1 <- read.table(text=lines[seq(start, end-1, by=2)], colClasses = rep("numeric", 10))
  line2 <- read.table(text=lines[seq(start+1, end, by=2)], colClasses = rep("numeric", 11))
  tmp <- cbind(line1, line2)
  tmp <- tmp[,c(1:3, 7, 11, 15, 19:21)]
  colnames(tmp) <- c("day_number", "minute_number", "dw_solar", "direct_n", "diffuse", "dw_longwave", "temp", "rh", "pressure") # "SWD", "DIR", "DIF", "LWD", "T2", "RH", "PoPoPoPo"
  tmp <- dplyr::as_tibble(tmp)
  tmp <- tmp %>%
    mutate(., Time = ymd_hm(paste0(substr(file, 6, 7), "-", substr(file, 4, 5), "-", tmp$day_number, " ", tmp$minute_number%/%60, ":", tmp$minute_number%% 60))) %>%
    select(., Time, everything()) %>%
    dplyr::select(-matches("day_number|minute_number"))

  #generate complete time sequence and bind data
  month <- as.Date(paste0("01",substr(file, 4, 7)), "%d%m%y") #which month?
  n_days <- as.numeric(lubridate::days_in_month(month)) #number of days in that month?
  res <- as.numeric(tmp$Time[2]- tmp$Time[1]) #check data resolution, 1 or 3 min
  Time.complete <- tibble(Time = seq(from = ymd_hm(paste0(month, ' 00:00')), to = ymd_hm(paste0(month, ' 00:00')) + 1440*60*(n_days)-60*res, by = paste0(res, ' min'))) #get complete time stamp
  tmp <- Time.complete %>% left_join(., tmp, by = "Time")

  #get the Linke turbidity for the station
  utils::data("BSRN.loc")
  stn <- substr(file, 1, 3)
  stn <- match(stn, BSRN.loc$stn)
  LT <- as.numeric(BSRN.loc[stn, which(substr(colnames(BSRN.loc), 1, 3)=="LTF")])

  # compute (SolarData) QC parameters, note the time shift in solpos
  mon <- lubridate::month(tmp$Time) # month
  solpos <- calZen(Tm = tmp$Time - res*30, lat = BSRN.loc$lat[stn], lon = BSRN.loc$lon[stn], tz = 0, LT[mon], alt = BSRN.loc$elev[stn])
  # append the solpos parameters
  tmp <- tmp %>%
    mutate(zen = solpos$zenith) %>%
    mutate(Ics = solpos$Ics) %>%
    mutate(Icsb = solpos$Icsb) %>%
    mutate(Icsd = solpos$Icsd) %>%
    mutate(Ioh = solpos$Ioh)
  tmp <- tmp %>%
    mutate(Mu0 = ifelse(tmp$zen > 90, 0, cos(radians(tmp$zen)))) %>%
    mutate(Sa = solpos$Io)
  # calculate closure
  tmp <- tmp %>%
    mutate(sum = tmp$diffuse + tmp$Mu0*tmp$direct_n)
  # calculate Rayleigh limit
  tmp <- tmp %>%
    mutate(RL = 209.3*tmp$Mu0 - 708.38*(tmp$Mu0)^2 + 1128.7*tmp$Mu0^3 -911.2*tmp$Mu0^4 + 287.85*tmp$Mu0^5+0.046725*tmp$Mu0*ifelse(tmp$pressure>0 & !is.na(tmp$pressure), tmp$pressure, as.numeric(NA))) #if pressure is missing or neg

  #assign NA to missing values, the missing code is -999 or -99.9
  tmp <- tmp %>%
    mutate_if(is.numeric, function(x) ifelse(x == -999 | x==-99.9, as.numeric(NA), x))

  #perform basic QC on 1-min or 3-min data
  if(use.qc)
  {
    if(is.null(test))
      stop("You must specify the test(s), avialable ones include 'phy', 'ext', 'closr', 'dr', 'clim', 'all'.")
    tmp <- QC.Basic.v2(tmp, test)
  }

  data_all <- tmp

  #set nighttime to zero and drop the unused columns
  data_all <- data_all %>%
    mutate(dw_solar = ifelse(data_all$dw_solar<0 & data_all$zen>90, 0, data_all$dw_solar)) %>%
    mutate(direct_n = ifelse(data_all$direct_n<0 & data_all$zen>90, 0, data_all$direct_n)) %>%
    mutate(diffuse = ifelse(data_all$diffuse<0 & data_all$zen>90, 0, data_all$diffuse)) %>%
    dplyr::select(-matches("Mu0|Sa|sum|RL"))

  #aggregate
  if(use.agg){
    if(res == 3 & agg < 3)
    {
      stop("you have 3-min data, 'agg' must be at least 3\n")
    }else if(res == 3 & agg == 3){
      cat("you have 3-min data, no aggregation required\n")
    }else if(res == 1 & agg == 1){
      cat("you have 1-min data, no aggregation required\n")
    }else if(agg>res){
      #cat("Aggregating files ...\n")
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
        data_all[remove, c(2:4)] <- NA
      }
      data_all <- data_all %>%
        group_by(Time) %>%
        summarise_all(.funs = mean, na.rm = TRUE)
    }
  }

  data_all
}


QC.Basic.v2 <- function(df, test)
{
  #check physical limit and flag
  if("phy" %in% test)
  {
    df <- df %>%
      mutate(qc_phy_G = ifelse(df$dw_solar > df$Sa*1.5*(df$Mu0)^1.2 + 100 | df$dw_solar < -4, 1, 0)) %>%
      mutate(qc_phy_D = ifelse(df$diffuse > df$Sa*0.95*(df$Mu0)^1.2 + 50 | df$diffuse < -4, 1, 0)) %>%
      mutate(qc_phy_I = ifelse(df$direct_n*df$Mu0 > df$Sa*df$Mu0 | df$direct_n*df$Mu0 < -4, 1, 0))

    df <- df %>%
      mutate(dw_solar = ifelse(df$qc_phy_G == 1 | is.na(df$qc_phy_G), as.numeric(NA), df$dw_solar)) %>%
      mutate(diffuse = ifelse(df$qc_phy_D == 1 | is.na(df$qc_phy_D), as.numeric(NA), df$diffuse)) %>%
      mutate(direct_n = ifelse(df$qc_phy_I == 1 | is.na(df$qc_phy_I), as.numeric(NA), df$direct_n)) %>%
      dplyr::select(-matches("qc_phy_G|qc_phy_D|qc_phy_I"))
  }

  #check extreme-rare limit and flag
  if("ext" %in% test | "all" %in% test)
  {
    df <- df %>%
      mutate(qc_ext_G = ifelse(df$dw_solar > df$Sa*1.2*(df$Mu0)^1.2 + 50 | df$dw_solar < -2, 1, 0)) %>%
      mutate(qc_ext_D = ifelse(df$diffuse > df$Sa*0.75*(df$Mu0)^1.2 + 30 | df$diffuse < -2, 1, 0)) %>%
      mutate(qc_ext_I = ifelse(df$direct_n*df$Mu0 > df$Sa*0.95*(df$Mu0)^1.2 + 10 | df$direct_n*df$Mu0 < -2, 1, 0))

    df <- df %>%
      mutate(dw_solar = ifelse(df$qc_ext_G == 1 | is.na(df$qc_ext_G), as.numeric(NA), df$dw_solar)) %>%
      mutate(diffuse = ifelse(df$qc_ext_D == 1 | is.na(df$qc_ext_D), as.numeric(NA), df$diffuse)) %>%
      mutate(direct_n = ifelse(df$qc_ext_I == 1 | is.na(df$qc_ext_I), as.numeric(NA), df$direct_n)) %>%
      dplyr::select(-matches("qc_ext_G|qc_ext_D|qc_ext_I"))
  }

  #check closure
  if("closr" %in% test | "all" %in% test)
  {
    df <- df %>%
      mutate(qc_closr = ifelse(df$dw_solar > 50 & df$zen < 75 & abs((df$sum-df$dw_solar)/df$dw_solar)*100 > 8 , 1, 0))
    df <- df %>%
      mutate(qc_closr = ifelse(df$dw_solar > 50 & df$zen > 75 & df$zen < 93 & abs((df$sum-df$dw_solar)/df$dw_solar)*100 > 15, 1, df$qc_closr))
    # the question is: if the test is indeterministic, should we remove the data point
    # If the answer is "yes", then uncomment the line below, otherwise, comment the line below
    # df <- df %>%
    #   mutate(qc_closr = ifelse(is.na(df$dw_solar) | is.na(df$sum), as.numeric(NA), df$qc_closr))

    df <- df %>%
      mutate_at(.vars = c("dw_solar", "direct_n", "diffuse"), function(x) ifelse(df$qc_closr == 1 | is.na(df$qc_closr), as.numeric(NA), x)) %>%
      dplyr::select(-matches("qc_closr"))
  }

  #check diffuse ratio
  if("dr" %in% test | "all" %in% test)
  {
    df <- df %>%
      mutate(qc_d_ratio = ifelse(df$diffuse/df$dw_solar > 1.05 & df$dw_solar>50 & df$zen < 75, 1, 0))
    df <- df %>%
      mutate(qc_d_ratio = ifelse(df$diffuse/df$dw_solar > 1.10 & df$dw_solar>50 & df$zen > 75, 1, df$qc_d_ratio))
    df <- df %>%
      mutate(qc_d_ratio = ifelse(is.na(df$dw_solar) | is.na(df$diffuse), as.numeric(NA), df$qc_d_ratio))

    df <- df %>%
      mutate_at(.vars = c("dw_solar", "diffuse"), function(x) ifelse(df$qc_d_ratio == 1 | is.na(df$qc_d_ratio), as.numeric(NA), x)) %>%
      dplyr::select(-matches("qc_d_ratio"))
  }

  #check climatology
  if("clim" %in% test | "all" %in% test)
  {
    df <- df %>%
      mutate(qc_clim1 = ifelse(df$diffuse/df$dw_solar > 0.85 & df$dw_solar/df$Ics > 0.85 & df$diffuse>50, 1, 0)) %>%
      mutate(qc_clim2 = ifelse(df$diffuse<df$RL-1 & df$diffuse/df$dw_solar < 0.80 & df$dw_solar>50, 1, 0))

    df <- df %>%
      mutate(qc_clim1 = ifelse(is.na(df$dw_solar) | is.na(df$diffuse), as.numeric(NA), df$qc_clim1)) %>%
      mutate(qc_clim2 = ifelse(is.na(df$dw_solar) | is.na(df$diffuse), as.numeric(NA), df$qc_clim2))

    df <- df %>%
      mutate_at(.vars = c("dw_solar", "diffuse"), function(x) ifelse(df$qc_clim1 == 1 | is.na(df$qc_clim1) | df$qc_clim2 == 1 | is.na(df$qc_clim2), as.numeric(NA), x)) %>%
      dplyr::select(-matches("qc_clim1|qc_clim2"))
  }

  df
}

