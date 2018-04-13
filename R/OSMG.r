#' @import insol

#################################################################################
# function to calculate sun position and extraterrestial irradiance
#################################################################################
calZen <- function(Tm, lat, lon, tz, LT, alt = 0)
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

