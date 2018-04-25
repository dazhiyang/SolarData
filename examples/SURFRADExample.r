# In this example, I show how to validate the Engerer2 (Engerer, SE, 2015) separation model using SURFRAD data

library(SolarData) #load the package
data("SURFRAD.loc")

#################################################################################
# two separation models, and normalize RMSE function
#################################################################################
Erbs <- function(Kt)
{
  if(Kt <= 0.22)
  {
    Kd <- 1- 0.09*Kt
  }else if(Kt > 0.22 & Kt <= 0.80){
    Kd <- 0.9511 - 0.1604*Kt + 4.388*(Kt)^2 - 16.638*(Kt)^3 + 12.336*(Kt)^4
  }else{
    Kd <- 0.165
  }
  Kd
}

Engerer2 <- function(Kt, AST, theta_z, delta_Ktc, Kde)
{
  C <- 4.2336e-2
  beta0 <- -3.7912; beta1 <- 7.5479; beta2 <- -1.0036e-2;
  beta3 <- 3.1480e-3; beta4 <- -5.3146; beta5 <- 1.7073;

  C + (1-C)/(1+exp(beta0+beta1*Kt + beta2*AST + beta3*theta_z + beta4*delta_Ktc)) + beta5*Kde
}

nRMSE <- function(meas, pred)
{
  sqrt(mean((meas-pred)^2))/sqrt(mean(meas^2))*100
}

#################################################################################
# read raw data, filter data, and calculate some required parameters
#################################################################################
# NOTE: This part assumes that you have the 365 daily files from bon downloaded
directory = "/Volumes/Macintosh Research/Data/surfrad/raw/bon/2016"
setwd(directory)
files <- dir() # get all the files in the directory
stn <- substr(files[1], 1, 3)
data <- SURFRAD.read(files, use.original.qc = FALSE, use.qc = TRUE, test = c("all"), directory = directory, agg = 1) # use all the tests, since we want the most accurate data points
tz <- 0 #since the time is in UTC
loc <- SURFRAD.loc[match(stn, SURFRAD.loc$stn),]
lon <- loc$lon
elev <- loc$elev

data <- data %>%
  filter(complete.cases(.)) %>% # rm NA rows
  filter(zen < 85)  # see section 3.3 in Gueymard and Ruiz-Arias (2016)

theta_z <- data$zen
Kt <- data$dw_solar/data$Ioh # clearness index, Kt, be careful, this is NOT clear-sky index
Kde <- pmax(0, 1-data$Ics/data$dw_solar) # see Eq. (32) of Engerer 2015
Ktc <- data$Ics/data$Ioh # clear-sky value of the clearness index
delta_Ktc <- Ktc-Kt # Eq. (30) of Engerer 2015
#calculate apparent solar time
jd <- JD(data$Time) # Julian Day from POSIXct
ha <- hourangle(jd, lon, tz) # hour angle
AST <- degrees(ha)/15 + 12

#################################################################################
# separation models
#################################################################################
#predicted diffuse fraction
Kd_Erbs <- sapply(Kt, Erbs)
Kd_Engerer2 <- Engerer2(Kt, AST, theta_z, delta_Ktc, Kde)
#predicted DNI
DNI_Erbs <- (data$dw_solar - Kd_Erbs * data$dw_solar)/cos(radians(theta_z))
DNI_Engerer2 <- (data$dw_solar - Kd_Engerer2 * data$dw_solar)/cos(radians(theta_z))
#error in percentage, comparable to Gueymard's report at bon, see his supplementary material mmc2 (spreadsheet)
nRMSE(meas = data$direct_n, pred = DNI_Erbs)
nRMSE(meas = data$direct_n, pred = DNI_Engerer2)
#plot
plot(Kt, data$diffuse/ data$dw_solar, pch = ".", ylim = c(0,1))
points(Kt, Kd_Erbs, col = 2, pch = ".")
points(Kt, Kd_Engerer2, col = 4, pch = ".")
