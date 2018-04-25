# Access and manipulate some publicly available solar data

There are many publicly available solar datasets. This package contains functions to download and manipulate these datasets. Currently available ones include: 
- NREL [Physical Solar Model (PSM)](https://nsrdb.nrel.gov/current-version) version 3, gridded satellite-derived irradiance data
- NREL [Oahu Solar Measurement Grid (OSMG)](https://midcdmz.nrel.gov/oahu_archive/), dense sensor network in Oahu, Hawaii
- NOAA [Surface Radiation (SURFRAD)](https://www.esrl.noaa.gov/gmd/grad/surfrad/), long-term high-resolution ground-based irradiance data
- NASA [Shuttle Radar Topography Mission (SRTM)](https://www2.jpl.nasa.gov/srtm/cbanddataproducts.html), digital elevation model data
- SoDa [Linke Turbidity Factor (LTF)](http://www.soda-pro.com/help/general-knowledge/linke-turbidity-factor), Linke turbidity data


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

This is an R package, so you need to install [R](https://www.r-project.org/) on your computer first. In addition, [RStudio](https://www.rstudio.com/) is an integrated development environment (IDE) for R; it is highly recommended.

### Installing

Once R and RStudio are installed. Open R or RStudio and install the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package, which allows you to install R package from GitHub

```
install.packages("devtools")
```

Load the package that you just installed

```
library("devtools")
```

Now, you can install the SolMod package, using

```
install_github("dazhiyang/SolarData")
```

## Running the tests

This code segment gives an example on how to run transposition modeling (horizontal to tilt) using a variety of models. (This is not up to date, I will update this section, as well as the package reference manual once the paper is accepted for publication)

```
library("SolarData")

#get SURFRAD data from Goodwin_Creek_MS (gwn) station, for the first three days in 2004
SURFRAD.get(station = 'Goodwin_Creek_MS', year = '2004', day_of_year = c(1:3))

#get PSM data for two locations
PSM.get(lat = c(42.05, 44), lon = c(-124.02, -110), api_key <- 'FVltdchrxzBCHiSNF6M7R4ua6BFe4j81fbPp8dDP', attributes <- 'ghi,dhi,dni,clearsky_dhi,clearsky_dni,clearsky_ghi,solar_zenith_angle', name = 'John+Smith', affiliation = 'Some+Institute', year = '2016', leap_year = 'true', interval = '30', utc = 'false', reason_for_use = 'research', email = 'yangdazhi.nus@gmail.com', mailing_list = 'false')

#get SRTM, i.e., digital elevation model, data for two boxes with resolution 3 arcsec
SRTM.list(3, want.plot = TRUE) #check available files
files <- c("Eurasia/N00E072.hgt.zip", "Eurasia/N00E073.hgt.zip")
SRTM.get(resolution = 3, files = files)
```

## License

This package is under the MIT license.
