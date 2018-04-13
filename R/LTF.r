
#' @importFrom tiff readTIFF

#################################################################################
# function to retrieve Linke turbidity from the monthly .tif map
# http://www.soda-is.com/eng/services/climat_free_eng.php#c5
# Remund et al. (2003)
#################################################################################
#' @export
LTF.get <- function(lon, lat, directory) #location in (lat, lon), path of the .tif files
{
  if(length(lat)!=length(lon))
    stop("'lon' and 'lat' must have the same length")

  setwd(directory)
  files <- dir(pattern = ".tif")
  lat_map <- seq(90, -90, by = -1/12)[-1] #1st (top) to last (bottom) pixel in the step of 1/12 degree
  lon_map <- seq(-180, 180, by = 1/12)[-1] #1st (left) to last (right) pixel in the step of 1/12 degree
  LT_map <- lapply(1:12, function(i) tiff::readTIFF(files[i], as.is = TRUE)/20) #12 monthly maps of size 2160x4320 pixels

  LT <- list()
  for(i in 1:length(lat))
  {
    LT[[i]] <- unlist(lapply(LT_map, function(x) x[which.min(abs(lat[i]-lat_map)), which.min(abs(lon[i]-lon_map))]))
  }
  LT_mat <- matrix(unlist(LT), nrow = 12)

  row.names(LT_mat) <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  colnames(LT_mat) <- paste("(", lon, ", ", lat, ")", sep = "")

  LT_mat
}

#lat <-c(10, 20,50)
#lon <- c(100, -100, 20)
#t <- LTF.get(lon, lat, directory = "/Volumes/Macintosh Research/Completed works/Solar_Energy_2017c_Entropy/Data/Linke Turbidity")


