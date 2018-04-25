# In this example, I show how to validate the Engerer2 (Engerer, SE, 2015) separation model using SURFRAD data

library(SolarData) #load the package
data("SURFRAD.loc")

#################################################################################
# read raw SURFRAD data, 2015 for fitting and 2016 for testing
#################################################################################
# NOTE: This part assumes that you have two years of daily files from dra downloaded
dir_SURFRAD <- "/Volumes/Macintosh Research/Data/surfrad/raw/dra"
setwd(file.path(dir_SURFRAD, "2015"))
files_2015 <- dir() # get all the files in the directory
setwd(file.path(dir_SURFRAD, "2016"))
files_2016 <- dir() # get all the files in the directory
stn <- substr(files_2016[1], 1, 3)
# read in the data as 1 min
# WHY we don't use the agg parameter?
# Since the PSM is satellite-derived data, i.e., snapshot, the high-res SURFRAD data needs to be rounded, instead of ceiling
SURFRAD_2015 <- SURFRAD.read(files_2015, use.original.qc = FALSE, use.qc = TRUE, test = c("all"), directory = file.path(dir_SURFRAD, "2015"), agg = 1)
SURFRAD_2015 <- SURFRAD_2015 %>%
  mutate(Time = lubridate::round_date(Time, "30 min")) %>%
  group_by(Time) %>%
  summarise_all(funs(mean), args = list(na.rm = TRUE))

SURFRAD_2016 <- SURFRAD.read(files_2016, use.original.qc = FALSE, use.qc = TRUE, test = c("all"), directory = file.path(dir_SURFRAD, "2016"), agg = 1)
SURFRAD_2016 <- SURFRAD_2016 %>%
  mutate(Time = lubridate::round_date(Time, "30 min")) %>%
  group_by(Time) %>%
  summarise_all(funs(mean), args = list(na.rm = TRUE))

loc <- SURFRAD.loc[match(stn, SURFRAD.loc$stn),]

#################################################################################
# download PSM data for the same two years
#################################################################################
dir_PSM <- "/Users/DYang/Dropbox/Working papers/R"
setwd(dir_PSM)
i <- 1; data <- list(); #initialize for the loop
for(yr in 2015:2016)
{
  # !!! REMEMBER to adjust the funcition parameters
  PSM.get(lon = loc$lon, lat = loc$lat, api_key = "YourAPIKey", attributes = "ghi,clearsky_ghi", name = "John+Smith", affiliation = "Some+Institute", year = as.character(yr), leap_year = "true", interval = "30", utc = "true", reason_for_use = "research", email = "email@gmail.com", mailing_list = "false", directory = dir_PSM)

  # read PSM data
  tmp <- read.csv(paste0("36.62373_-116.01947_", yr, ".csv"), header = TRUE, skip = 2)
  # get date.time for the PSM data
  tmp <- tmp %>%
    mutate(., Time = lubridate::ymd_hm(paste(paste(Year, Month, Day, sep = "-"), paste(Hour, Minute,  sep = ":"), sep = " "), tz = "UTC")) %>%
    dplyr::select(., -(1:5))
  # join the PSM and SURFRAD data
  tmp <- tmp %>%
    left_join(get(paste0("SURFRAD_", yr)), ., by = "Time")

  data[[i]] <- tmp %>%
    filter(complete.cases(.)) %>% # NA filter
    filter(zen < 80) %>% # zenith angle filter
    filter(dw_solar > 30 & GHI > 30) #filter out small GHI values
  i <- i+1
}

#################################################################################
# site adaptation
#################################################################################
p1 <- ggplot(data = data[[2]], aes(x = dw_solar, y = GHI)) +
  stat_binhex(bins = 30) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.2) +
  geom_abline(intercept = 0, slope = 1.1, linetype = "dashed", color = "black", size = 0.2) +
  geom_abline(intercept = 0, slope = 0.9, linetype = "dashed", color = "black", size = 0.2) +
  xlab(expression(paste("SURFRAD data (-88.37309", degree, ", 40.05192", degree, ")"))) +
  ylab(expression(paste("PSM data (-88.37309", degree, ", 40.05192", degree, ")"))) +
  scale_fill_gradientn(colours = c("lightblue","yellow","red"), name = "Count") +
  theme_gray() +
  theme(plot.margin = unit(c(0,-0.4,0,0.1), "lines"), legend.position = "none", text = element_text(family = "Times", size = 7))

# now we have some nice data to be site adapted
# and a linear model seems sufficient
p1

# build a linear model
PMS_GHI_2015 <- data[[1]]$GHI
SURFRAD_GHI_2015 <- data[[1]]$dw_solar
reg <- lm(SURFRAD_GHI_2015~PMS_GHI_2015)

# adaptation
PMS_GHI_2016 <- data[[2]]$GHI
SURFRAD_GHI_2016 <- data[[2]]$dw_solar
PMS_GHI_2016_pred <- PMS_GHI_2016*coef(reg)[2] + coef(reg)[1]

#error
nMBE <- function(meas, pred)
{
  mean(meas-pred)/mean(meas)*100
}
nMBE(meas = SURFRAD_GHI_2016, pred = PMS_GHI_2016)
nMBE(meas = SURFRAD_GHI_2016, pred = PMS_GHI_2016_pred) #somewhat smaller

nRMSE <- function(meas, pred)
{
  sqrt(mean((meas-pred)^2))/sqrt(mean(meas^2))*100
}
nRMSE(meas = SURFRAD_GHI_2016, pred = PMS_GHI_2016)
nRMSE(meas = SURFRAD_GHI_2016, pred = PMS_GHI_2016_pred) # lol, ok, maybe linear doesn't work

