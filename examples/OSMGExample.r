# In this example, I show how to obtain the correlation functions mentioed in Arias-Castro et al. 2014

library(SolarData) #load the package

#################################################################################
# correlation functions used in Arias-Castro et al. 2014
#################################################################################
ACM <- function(a, b, cs, t, r)
{
  A <- function(D) {ifelse(D>1, 0, 2*(acos(D)-sqrt(D^2-D^4))/pi)}
  D.aw = a/2/r    #Dalong
  D.cw = b/2/r    #Dcross
  D.cl = cs*t/2/r #Dcloud

  denom = 2*(0.5-0.5^(2-A(D.cl)))
  numer = 2*0.5^(2-A(sqrt(D.aw^2+D.cw^2))) - 0.5^(2-A(sqrt((D.aw-D.cl)^2+D.cw^2))) - 0.5^(2-A(sqrt((D.aw+D.cl)^2+D.cw^2)))
  numer/denom
}

Lonij <- function(a, b, cs, t)
{
  #these parameter follows Lonij et al. 2013, except that sigma.x follow Arias-Castro et al. 2014
  q <- 0.65; A <- 0.024; sigma.t <- 6*3600; sigma.x <- 1000
  cs.aw <- cs; cs.cw <- 0; #assume cloud speed in the direction of distance
  f <- function(dx, dy, dt){
    A*exp(-sqrt((dx-cs.aw*dt)^2 + (dy-cs.cw*dt)^2)/sigma.x)*exp(-(abs(dt)/sigma.t)^q)
  }
  f(a, b, t) + f(a, b, -t) -2*f(a, b, 0)
}

Lave <- function(d, cs, t)
{
  #exp(-d/0.5/cs/t) #error here! should be exp(-d/cs/t), then can obtain figure 4
  exp(-d/cs/t)
}

Perez <- function(d, cs, t)
{
  exp(log(0.2)*d/1.5/cs/t)
}


#################################################################################
# read OSMG location data, and get the inter-station distances
#################################################################################
data("OSMG.loc") #load data from package
loc <- OSMG.loc[-c(9, 11),] #rm tilted loc, assign a shorter name
#along-wind dist mat
dist_aw <- SolarData:::GetAlongWindDist(loc$Longitude, loc$Latitude, wind.dir=60)
aw_order <- order(dist_aw[14,]) #14 is AP7, order in increasing dist
loc = loc[aw_order,] #reorder the locations in along-wind direction
dist <- SolarData:::GetGeoDist(loc$Longitude, loc$Latitude)
dist_aw <- dist_aw[aw_order, aw_order]
#cross-wind dist mat
dist_cw <- SolarData:::GetAlongWindDist(loc$Longitude, loc$Latitude, wind.dir=-30)
dist_aw = dist_aw*1000; dist_cw = dist_cw*1000; dist = dist*1000;

# Define along-wind and cross-wind pairs
# AP7 = 1, AP4 = 2, AP3 = 3, AP2 = 4, AP6 = 5, DH5 = 6, AP1 = 7, DH2 = 8, AP5 = 9, DH3 = 10, DH4 = 11, DH1 = 12, DH7 = 13, DH10 = 14, DH11 = 15, DH9 = 16, DH6 = 17, DH8 = 18
pair_aw = matrix(c( 7,  6, 10, 4,  7, 11,  4,  6, 2, 1,  3,  1,  1,
                    10, 11, 14, 7, 14, 17, 10, 17, 9, 3, 15, 15, 18), ncol = 2, byrow = FALSE)
pair_cw = matrix(c(11, 6, 13, 7, 6, 8, 4,
                   10, 7, 14, 9, 9, 9, 5), ncol = 2, byrow = FALSE)

#read raw data and aggregate to 10 s
#NOTE: This part assumes that you have the 26 files (13 days mentioned in the paper) downloaded in the two folders
dir_LI200 <- "/Volumes/Macintosh Research/Data/Oahu/raw"
dir_RSR <- "/Volumes/Macintosh Research/Data/Oahu/raw AP2"
setwd(dir_LI200)
files <- dir()
data <- OSMG.read(files, directory_LI200 = dir_LI200, directory_RSR = dir_RSR, clear_sky = TRUE, AP2 = TRUE, agg = 10)

#data filters
data <- data %>%
  filter(lubridate::hour(Time) > 8) %>%
  filter(lubridate::hour(Time) < 15) %>%
  dplyr::select(., c(1, 7:24)) #only horizontal GHI data

#aggregate data into 10, 20, 30, 60 , 120, 180s
data_all <- list()
time_scale <- c(10, 20, 30, 60, 90, 120, 180, 300)
for(i in 1:length(time_scale))
{
  ts <- time_scale[i]
  data_all[[i]] <- data %>%
    mutate(Time = SolarData:::ceiling.time(Time, ts)) %>%
    group_by(Time) %>%
    summarise_all(funs(mean), na.rm = TRUE)
}

#################################################################################
# compute correlation (points), and fit the correlation functions (lines)
#################################################################################
data_point_aw = data_point_cw = data_line <- NULL
for(j in 1:length(time_scale))
{
  ts <- time_scale[j]
  data <- apply(data_all[[j]][,-1], 2, function(x) diff(x))
  data <- data[,aw_order] #order the data in along-wind direction
  C = cor(data, use = "pairwise.complete.obs") #correlation

  #select the along-wind pairs, and extract their dist and corr
  mat_aw <- matrix(FALSE, 18, 18)
  for(i in 1:nrow(pair_aw))
    mat_aw[pair_aw[i, 2], pair_aw[i, 1]] <- TRUE
  C_aw <- C[mat_aw]; d_aw <- dist[mat_aw]; a_aw <- dist_aw[mat_aw]; b_aw <- dist_cw[mat_aw];

  #select the cross-wind pairs, and extract their dist and corr
  mat_cw <- matrix(FALSE, 18, 18)
  for(i in 1:nrow(pair_cw))
    mat_cw[pair_cw[i, 2], pair_cw[i, 1]] <- TRUE
  C_cw <- C[mat_cw]; d_cw <- dist[mat_cw]; a_cw <- dist_cw[mat_cw]; b_cw <- dist_cw[mat_cw];

  data_point_aw = rbind(data_point_aw, data.frame(x = d_aw, y = C_aw, timescale = paste("Avg-", ts, "s", sep = "")))
  data_point_cw = rbind(data_point_cw, data.frame(x = d_cw, y = C_cw, timescale = paste("Avg-", ts, "s", sep = "")))

  d_plot <- seq(0, max(dist), length.out = 100)
  line_Lave <- Lave(d = d_plot, cs = 10, t = ts)
  line_Perez <- Perez(d = d_plot, cs = 10, t = ts)
  line_Lonij_aw <- Lonij(a = d_plot, b = 0, cs = 10, t = ts)/Lonij(a = 0, b = 0, cs = 10, t = ts)
  line_Lonij_cw <- Lonij(a = 0, b = d_plot, cs = 10, t= ts)/Lonij(a = 0, b = 0, cs = 10, t= ts)
  line_ACM_aw <- ACM(a = d_plot, b = 0, cs = 10, t = ts, r = 1000)
  line_ACM_cw <- ACM(a = 0, b = d_plot, cs = 10, t = ts, r = 1000)

  tmp <- data.frame(x = rep(d_plot, 6), y = c(line_Lave, line_Perez, line_Lonij_aw, line_ACM_aw, line_Lonij_cw, line_ACM_cw),
                    model = rep(c("Lave", "Perez", "Lonij_aw", "ACM_aw", "Lonij_cw", "ACM_cw"), each = 100), timescale = paste("Avg-", ts, "s", sep = ""))
  data_line <- rbind(data_line, tmp)
}


#################################################################################
# plot, replicate Fig.4 in Arias-Castro et al. 2014
#################################################################################
data_line$model <- factor(data_line$model, levels = c("ACM_aw", "ACM_cw", "Lonij_aw", "Lonij_cw", "Lave", "Perez"))
p <- ggplot(data_line) +
     geom_line(aes(x = x, y = y, linetype = model, colour = model)) +
     scale_x_continuous(limits = c(0, 1000), expand = c(0, 0)) +
     scale_y_continuous(limits = c(-0.4,1.2),expand = c(0, 0)) +
     xlab("Distance [m]") +
     ylab("Ramp correlation [dimensionless]") +
     facet_wrap(~timescale, ncol = 2) +
     scale_linetype_manual(values = c("solid", "solid", "dashed", "dashed", "solid", "dashed")) +
     scale_color_manual(values = c("blue", "red", "blue", "red", "black", "black")) +
     scale_size_manual(values = c(1,1, 0.5, 0.5, 0.5, 0.5)) +
     geom_point(data = data_point_aw, aes(x = x, y = y), size = 1, color = "blue") +
     geom_point(data = data_point_cw, aes(x = x, y = y), size = 1, color = "red") +
     theme_gray() +
     theme(plot.margin = unit(c(0.5,0.5,0,0), "lines"),
        panel.spacing = unit(0.1, "lines"),
        text = element_text(size = 9),
        legend.position = "right",
        strip.text.x = element_text(margin = margin(0.2,0,0.2,0, "lines")))
p #exactly like Fig.4 in Arias-Castro et al. 2014

