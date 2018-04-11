#' @import ggplot2
#' @importFrom dplyr %>%


#' @export
SRTM.list <- function(resolution, want.plot = TRUE)
{
  cat("retrieving directory listing from 'https://dds.cr.usgs.gov/'...")
  res <- paste0("SRTM", resolution)
  URL <- paste("https://dds.cr.usgs.gov/srtm/version2_1", res, sep = "/")

  page <- textreadr::read_html(URL)
  dir <- page[-c(1,2)]

  sub_dir <- list()
  j = 1
  for(i in 1:length(dir))
  {
    if(substr(dir[i], nchar(dir[i]), nchar(dir[i])) == "/")
    {
      sub_dir[[j]] <- textreadr::read_html(paste(URL, dir[i], sep = "/"))[-c(1,2)]
      names(sub_dir)[[j]] <- gsub('.{1}$', '', dir[i])
      j <- j+1
    }
  }

  files <- sapply(1:length(sub_dir), function(x) paste(names(sub_dir)[x], sub_dir[[x]], sep = "/"))
  cat("Done~")

  if(want.plot)
  {
    libs <- c("ggplot2", "ggmap")
    invisible(lapply(libs, library, character.only = TRUE))
    data.plot<- NULL
    for(i in 1:length(files))
    {
      loc <- strsplit(files[[i]], split = "/")
      group <- loc[[1]][1]
      lat <- sapply(loc, function(x) substr(x[[2]], 1, 3))
      lat <- ifelse(substr(lat, 1, 1) == "N", as.numeric(substr(lat, 2, 3)), -as.numeric(substr(lat, 2, 3)))
      lon <- sapply(loc, function(x) substr(x[[2]], 4, 8))
      lon <- ifelse(substr(lon, 1, 1) == "E", as.numeric(substr(lon, 2, 4)), -as.numeric(substr(lon, 2, 4)))
      data.plot <- rbind(data.plot, data.frame(lat = lat, lon = lon, group = group))
    }

    WorldData <- ggplot2::map_data('world')
    WorldData %>% dplyr::filter(WorldData$region != "Antarctica") -> WorldData
    WorldData <- ggplot2::fortify(WorldData)

    p <- ggplot()
    p <- p + geom_map(map=WorldData, fill=NA, colour="gray30", size=0.5) +
      geom_tile(data=data.plot, aes(lon+0.5,lat+0.5, color = group), alpha=0, size = 0.1) +
      coord_fixed() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      xlab(expression(paste("Longitude [", degree, "]", sep = ""))) +
      ylab(expression(paste("Latitude [", degree, "]", sep = ""))) +
      theme_bw() +
      theme(plot.margin = unit(c(0,0,0,0), "lines"), legend.position = "bottom", legend.box.spacing = unit(c(-0.5,0,0.5,0), "lines"), text = element_text(size = 7), legend.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.x = element_text(margin = margin(0.2,0,0.2,0, "lines"))    )
    print(p)
  }

  files
}

#listSRTM(1)

#' @export
SRTM.get <- function(resolution, files, directory = "data-raw")
{
  res <- paste0("SRTM", resolution)

  #if(missing(zone))
    #stop("zone must be defined! Use function 'listSRTM' to check for available zones")

  if(missing(files))
    stop("file names must be defined! Use function 'listSRTM' to check for available files")

  if(!file.exists(directory))
    dir.create(directory)

  for(i in 1:length(files))
  {
    URL <- paste("https://dds.cr.usgs.gov/srtm/version2_1", res, files[i], sep = "/")
    dest <- paste(directory, strsplit(files[i], split = "/")[[1]][2], sep = "/")
    utils::download.file(URL, destfile = dest)
    utils::unzip(dest, exdir = directory)
    unlink(dest)
  }

}

#' @export
SRTM.read <- function(files, as.data.frame = FALSE)
{
  cat("Warning: the files are big, you may want to read one file at a time\n")
  if(length(files)>1)
  {
    tmp <- list()
    for(i in 1:length(files))
    {
      cat("Reading", files[i], "...\n")
      tmp[[i]] <- raster::raster(files[i])
    }
    elevation <- do.call(raster::merge, tmp)
  }else{
    elevation <- raster::raster(files[i])
  }

  #set voids to NA
  raster::values(elevation)[raster::values(elevation) < -10994 | raster::values(elevation) > 8848] = NA

  if(as.data.frame)
  {
    elevation <- raster::as.data.frame(elevation, xy = TRUE)
  }
  elevation
}



#setwd("/Users/DYang/Git/SolarData/data-raw")
#files <- dir()
#raster::image(elevation)

#raster::image(df2)

#zone <- "North_America"

#files <- c("North_America/N16W090.hgt.zip", "North_America/N15W090.hgt.zip", "North_America/N15W091.hgt.zip", "North_America/N15W092.hgt.zip")
#SRTM.get(resolution = 3, files = files)

