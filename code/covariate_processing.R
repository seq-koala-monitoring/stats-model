# -------------------------------------------------------------------
#                 UPDATE AND PROCESS COVARIATES
# -------------------------------------------------------------------

# This R code updates and processes covariates required for the Bayesian state-space model 
# used to estimate koala densities across Southeast Queensland.
#
# The basic workflow involves checking if the input folder already has the latest covariate file. 
# If it doesn't, the code will automatically download, process, and save the most recent file
# the correct location.
#
# Before anything, start fresh by clicking Session > Restart R... > Yes in the pop-up (do not save anything if asked)
#
# To get started, please, select all lines by pressing Ctrl + A on a Windows PC or Command + A on a MAc. 
# Then, run these lines by pressing Ctrl + Enter on a Windows PC or Command + Return on a Mac.
# Alternatively, click on the Source button near the top right corner of this script
#
# When you see "THIS CODE HAS FINISHED" in the Console panel (usually at the bottom left),
# you're ready to the modelling stage. Go to File > Open file... > locate the file named data_processing.R > Open.
# This file prepares all the data for the Bayesian state-space model
#
# PS: 1) Disregard any warnings on the task bar about packages that are not installed 
#     2) This code may take anywhere from a few minutes to days to run, 
#        depending on how many files need updating and your computer's specifications.


# -------------------------------------------------------------------
# Clean global environment
rm(list=ls())
try(dev.off(dev.list()["RStudioGD"]), silent=TRUE)
gc()

# read utility functions
source("code/functions.R")

# install packages
packages <- c("terra","sf","rvest","pbapply","httr","stringi","foreach","doParallel","tidyverse")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, quiet = T)
invisible(lapply(packages, library, character.only = TRUE))

# get tokens needed
source("keys/apis.R")

# load study area
StudyArea <- st_read("input/covariates/raw_data/study_area/study_area.gdb", layer = "study_area")

# get bounding box with a 20km buffer around the study area
Bounding <- st_bbox(st_buffer(StudyArea, 20000)) %>% st_as_sfc() %>% st_sf(attribute = 1) %>% vect()

# process covariate layers

# get digital elevation model data - located in "raw_data" folder
if (!file.exists("input/covariates/output/htele.tif") | !file.exists("input/covariates/output/htslo.tif")
          | !file.exists("input/covariates/output/htrug.tif")) {
  DEM <- rast(paste(getwd(), "input/covariates/raw_data/dem/srtm-1sec-dem-v1-COG.tif",sep = "")) %>%
        project(rast(Bounding, resolution = 30))

  # Elevation
  Elev <- DEM
  writeRaster(Elev, "input/covariates/output/htele.tif", overwrite = TRUE)

  # Slope
  Slope <- terrain(Elev, "slope")
  writeRaster(Slope, "input/covariates/output/htslo.tif", overwrite = TRUE)

  # Terrain Ruggedness Index
  TRI <- terrain(Elev, "TRI")
  writeRaster(TRI, "input/covariates/output/htrug.tif", overwrite = TRUE)
}

# get soil data - accessed directly from TERN based on COGs (need an API)

# get API key here
apikey <- paste0("apikey:", TERN_API)

# cation exchange capacity - mean values to 2 m
if (!file.exists("input/covariates/output/hscec.tif")) {
  CationExch <- rast(c(
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/CEC/CEC_000_005_EV_N_P_AU_TRN_N_20220826.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/CEC/CEC_005_015_EV_N_P_AU_TRN_N_20220826.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/CEC/CEC_015_030_EV_N_P_AU_TRN_N_20220826.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/CEC/CEC_030_060_EV_N_P_AU_TRN_N_20220826.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/CEC/CEC_060_100_EV_N_P_AU_TRN_N_20220826.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/CEC/CEC_100_200_EV_N_P_AU_TRN_N_20220826.tif")
                )) %>% project(rast(Bounding, resolution = 90)) %>% app(mean)
  writeRaster(CationExch, "input/covariates/output/hscec.tif", overwrite = TRUE)
}

# total phosphorous - mean values to 2 m
if (!file.exists("input/covariates/output/hspho.tif")) {
  Phos <- rast(c(
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/PTO/PTO_000_005_EV_N_P_AU_NAT_C_20140801.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/PTO/PTO_005_015_EV_N_P_AU_NAT_C_20140801.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/PTO/PTO_015_030_EV_N_P_AU_NAT_C_20140801.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/PTO/PTO_030_060_EV_N_P_AU_NAT_C_20140801.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/PTO/PTO_060_100_EV_N_P_AU_NAT_C_20140801.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/PTO/PTO_100_200_EV_N_P_AU_NAT_C_20140801.tif")
                )) %>% project(rast(Bounding, resolution = 90)) %>% app(mean)
  writeRaster(Phos, "input/covariates/output/hspho.tif", overwrite = TRUE)
}

# total nitrogen - mean values to 2 m
if (!file.exists("input/covariates/output/hsnit.tif")) {
  Nit <- rast(c(
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/NTO/NTO_000_005_EV_N_P_AU_NAT_C_20140801.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/NTO/NTO_005_015_EV_N_P_AU_NAT_C_20140801.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/NTO/NTO_015_030_EV_N_P_AU_NAT_C_20140801.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/NTO/NTO_030_060_EV_N_P_AU_NAT_C_20140801.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/NTO/NTO_060_100_EV_N_P_AU_NAT_C_20140801.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/NTO/NTO_100_200_EV_N_P_AU_NAT_C_20140801.tif")
                )) %>% project(rast(Bounding, resolution = 90)) %>% app(mean)
  writeRaster(Nit, "input/covariates/output/hsnit.tif", overwrite = TRUE)
}
# available water capacity - mean values to 2 m
if (!file.exists("input/covariates/output/hswat.tif")) {
  AWat <- rast(c(
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/AWC/AWC_000_005_EV_N_P_AU_TRN_N_20210614.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/AWC/AWC_005_015_EV_N_P_AU_TRN_N_20210614.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/AWC/AWC_015_030_EV_N_P_AU_TRN_N_20210614.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/AWC/AWC_030_060_EV_N_P_AU_TRN_N_20210614.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/AWC/AWC_060_100_EV_N_P_AU_TRN_N_20210614.tif"),
                  paste0("/vsicurl/https://", apikey, "@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/AWC/AWC_100_200_EV_N_P_AU_TRN_N_20210614.tif")
                )) %>% project(rast(Bounding, resolution = 90)) %>% app(mean)
  writeRaster(AWat, "input/covariates/output/hswat.tif", overwrite = TRUE)
}

# persistent green

# download all persistent green data from TERN if needed
# get URL
URL <- "https://data.tern.org.au/rs/public/data/landsat/seasonal_fractional_cover_v3/persistent_green/qld/"
pg <- read_html(URL)
# find all files that contain "tif" and start with "lztmre"
tifs <- grep("tif", html_attr(html_nodes(pg, "a[href^='lztmre']"), "href"), value = TRUE)
# save files
# save the current directory path for later
wd <- getwd()
# change working directory for the download
setwd(paste0(wd, "/input/covariates/raw_data/persistent_green"))
# download them all - do not download if file aready exists
invisible(pbsapply(tifs, function(tif_file) {
  if(!file.exists(tif_file)) {GET(URL %s+% tif_file, write_disk(tif_file))}
}))
setwd(wd)

# process persistent green data - taking the mid month of each 3 month epoch as the date of processed data
# get all .tif files in the relevant directory
PGFiles <- list.files("input/covariates/raw_data/persistent_green", pattern = "\\.tif$")

# set up parallel processing cluster
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)

# create raster of boundary box with 30m resolution and wrap() so it can be passed to foreach()
BoundingWrap <- wrap(rast(Bounding, resolution = 30, vals = 1))

# get each tif file, project, clip, and save
foreach (i = 1:length(PGFiles), .packages = c("tidyverse", "terra")) %dopar% {
  if (substr(sub(".*qld_m", "", PGFiles[i]), 5, 6) == "12") {
    Date <- paste0(as.character(as.numeric(substr(sub(".*qld_m", "", PGFiles[i]), 1, 4)) + 1), "01")
  } else {
    if (as.numeric(substr(sub(".*qld_m", "", PGFiles[i]), 5, 6)) < 9) {
      Date <- paste0(substr(sub(".*qld_m", "", PGFiles[i]), 1, 4), "0", as.character(as.numeric(substr(sub(".*qld_m", "", PGFiles[i]), 5, 6)) + 1))
    } else {
      Date <- paste0(substr(sub(".*qld_m", "", PGFiles[i]), 1, 4), as.character(as.numeric(substr(sub(".*qld_m", "", PGFiles[i]), 5, 6)) + 1))
    }
  }
  if (!file.exists(paste0("input/covariates/output/hhpgr", Date, ".tif"))) {
    Rast <- rast(paste0(getwd(), "/input/covariates/raw_data/persistent_green/", PGFiles[i])) %>% project(unwrap(BoundingWrap))
    NAflag(Rast) <- 255

    writeRaster(Rast, paste0("input/covariates/output/hhpgr", Date, ".tif"), overwrite = TRUE)

    # now get buffered average values round each cell for 1km, 2km, and 3km buffers
    Buff1 <- focalMat(Rast, 1000, "circle", fillNA = TRUE)
    Buff2 <- focalMat(Rast, 2000, "circle", fillNA = TRUE)
    Buff3 <- focalMat(Rast, 3000, "circle", fillNA = TRUE)
    focal(Rast, Buff1, filename = paste0("input/covariates/output/hhpgr1km", Date, ".tif"),
                    overwrite = TRUE, na.rm = TRUE)
    focal(Rast, Buff2, filename = paste0("input/covariates/output/hhpgr2km", Date, ".tif"),
                    overwrite = TRUE, na.rm = TRUE)
    focal(Rast, Buff3, filename = paste0("input/covariates/output/hhpgr3km", Date, ".tif"),
                    overwrite = TRUE, na.rm = TRUE)
  }
}

stopCluster(cl)

# land-use

# load land-use for 1999 and 2017 - this uses the secondary ALUM land-use category, reclassifies to 7 categories, and then converts to raster at a 30m resolution
# note that the original land-use data was provided in a .gdb but was corrupted in some way so was first exported to a .shp file using ArcGIS Pro
if (!file.exists("input/covariates/output/htlus201706.tif") | !file.exists("input/covariates/output/htlus199906.tif")) {
  LandUse2017 <- vect("input/covariates/raw_data/land_use/QLD_LANDUSE_June_2019/QLD_LANDUSE_CURRENT_X.shp") %>% project(Bounding) %>% crop(Bounding) %>% rasterize(rast(Bounding, resolution = 30), field = "Secondary")
  # reclasify land use into: 1 = natural, 2  = production from natural, 3 = plantation, 4 = agriculture, 5 = intensive, 6 = wetland, 7 = water
  LULookup2017 <- read_csv("input/covariates/raw_data/land_use/land_use_2017_lookup_raw.csv")
  LandUse2017 <- classify(LandUse2017, LULookup2017)

  LandUse1999 <- vect("input/covariates/raw_data/land_use/QLD_LANDUSE_June_2019/QLD_LANDUSE_1999_X.shp") %>% project(Bounding) %>% crop(Bounding) %>% rasterize(rast(Bounding, resolution = 30), field = "Secondary")
  # reclasify land use into: 1 = natural, 2  = production from natural, 3 = plantation, 4 = agriculture, 5 = intensive, 6 = wetland, 7 = water
  LULookup1999 <- read_csv("input/covariates/raw_data/land_use/land_use_1999_lookup_raw.csv")
  LandUse1999 <- classify(LandUse1999, LULookup1999)

  # save rasters
  writeRaster(LandUse2017, "input/covariates/output/htlus201706.tif", overwrite = TRUE)
  writeRaster(LandUse1999, "input/covariates/output/htlus199906.tif", overwrite = TRUE)

  # intensive land-use buffers

  # reclassify so that intensive land-uses are 1 and everything else is 0
  Intensive2017 <- classify(LandUse2017, cbind(c(1, 2, 3, 4, 5, 6, 7), c(0, 0, 0, 0, 1, 0, 0)))
  Intensive1999 <- classify(LandUse1999, cbind(c(1, 2, 3, 4, 5, 6, 7), c(0, 0, 0, 0, 1, 0, 0)))

  # now get buffered average values round each cell for 1km, 2km, and 3km buffers

  # 2017

  Buff1 <- focalMat(Intensive2017, 1000, "circle", fillNA = TRUE)
  Buff2 <- focalMat(Intensive2017, 2000, "circle", fillNA = TRUE)
  Buff3 <- focalMat(Intensive2017, 3000, "circle", fillNA = TRUE)
  focal(Intensive2017, Buff1, filename = paste0("input/covariates/output/htilu1km", "201706", ".tif"),
                    overwrite = TRUE, na.rm = TRUE)
  focal(Intensive2017, Buff2, filename = paste0("input/covariates/output/htilu2km", "201706", ".tif"),
                    overwrite = TRUE, na.rm = TRUE)
  focal(Intensive2017, Buff3, filename = paste0("input/covariates/output/htilu3km", "201706", ".tif"),
                    overwrite = TRUE, na.rm = TRUE)

  # 1999

  Buff1 <- focalMat(Intensive1999, 1000, "circle", fillNA = TRUE)
  Buff2 <- focalMat(Intensive1999, 2000, "circle", fillNA = TRUE)
  Buff3 <- focalMat(Intensive1999, 3000, "circle", fillNA = TRUE)
  focal(Intensive1999, Buff1, filename = paste0("input/covariates/output/htilu1km", "199906", ".tif"),
                    overwrite = TRUE, na.rm = TRUE)
  focal(Intensive1999, Buff2, filename = paste0("input/covariates/output/htilu2km", "199906", ".tif"),
                    overwrite = TRUE, na.rm = TRUE)
  focal(Intensive1999, Buff3, filename = paste0("input/covariates/output/htilu3km", "199906", ".tif"),
                    overwrite = TRUE, na.rm = TRUE)
}

# land use mask

if (!file.exists("input/covariates/output/mask/lu2017_mask.tif") | !file.exists("input/covariates/output/mask/lu1999_mask.tif")) {

  # load land-use for 1999 and 2017 - this uses the tertiary ALUM land-use category, reclassifies to 2 categories, and then converts to raster at a 30m resolution
  # note that the original land-use data was provided in a .gdb but was corrupted in some way so was first exported to a .shp file using ArcGIS Pro
  LandUse2017_mask <- vect("input/covariates/raw_data/land_use/QLD_LANDUSE_June_2019/QLD_LANDUSE_CURRENT_X.shp") %>% project(Bounding) %>% crop(Bounding) %>% rasterize(rast(Bounding, resolution = 30), field = "Tertiary")
  LandUse1999_mask <- vect("input/covariates/raw_data/land_use/QLD_LANDUSE_June_2019/QLD_LANDUSE_1999_X.shp") %>% project(Bounding) %>% crop(Bounding) %>% rasterize(rast(Bounding, resolution = 30), field = "Tertiary")
  # save rasters
  writeRaster(LandUse2017_mask, "input/covariates/output/mask/lu2017_mask.tif")
  writeRaster(LandUse1999_mask, "input/covariates/output/mask/lu1999_mask.tif")
}

# climate

# download 6 monthly climate variables since 1990 if needed from BOM

# precipitation
for (i in 1990:as.numeric(format(Sys.Date(), "%Y"))) {
  # get download URLs and file names
  URL1 <- paste0("http://opendap.bom.gov.au:8080/thredds/fileServer/agcd/precip/total/r005/06month/", as.character(i),
                 "/precip_total_r005_", as.character(i - 1), "1001_", as.character(i), "0331.nc")
  File1 <- paste0("input/covariates/raw_data/precipitation", "/precip_total_r005_", as.character(i - 1), "1001_", as.character(i), "0331.nc")
  URL2 <- paste0("http://opendap.bom.gov.au:8080/thredds/fileServer/agcd/precip/total/r005/06month/", as.character(i),
                 "/precip_total_r005_", as.character(i), "0401_", as.character(i), "0930.nc")
  File2 <- paste0("input/covariates/raw_data/precipitation", "/precip_total_r005_", as.character(i), "0401_", as.character(i), "0930.nc")

  # download files
  if (!file.exists(File1)) {
    download.file(URL1, File1, mode = "wb")
  }
  if (!file.exists(File2)) {
    download.file(URL2, File2, mode = "wb")
  }
}

# mean mean temperature
for (i in 1990:as.numeric(format(Sys.Date(), "%Y"))) {
  # get download URLs and file names
  URL1 <- paste0("http://opendap.bom.gov.au:8080/thredds/fileServer/agcd/tmean/mean/r005/06month/", as.character(i),
                 "/tmean_mean_r005_", as.character(i - 1), "1001_", as.character(i), "0331.nc")
  File1 <- paste0("input/covariates/raw_data/temp_mean", "/tmean_mean_r005_", as.character(i - 1), "1001_", as.character(i), "0331.nc")
  URL2 <- paste0("http://opendap.bom.gov.au:8080/thredds/fileServer/agcd/tmean/mean/r005/06month/", as.character(i),
                 "/tmean_mean_r005_", as.character(i), "0401_", as.character(i), "0930.nc")
  File2 <- paste0("input/covariates/raw_data/temp_mean", "/tmean_mean_r005_", as.character(i), "0401_", as.character(i), "0930.nc")

  # download files
  if (!file.exists(File1)) {
    download.file(URL1, File1, mode = "wb")
  }
  if (!file.exists(File2)) {
    download.file(URL2, File2, mode = "wb")
  }
}

# mean max temperature
for (i in 1990:as.numeric(format(Sys.Date(), "%Y"))) {
  # get download URLs and file names
  URL1 <- paste0("http://opendap.bom.gov.au:8080/thredds/fileServer/agcd/tmax/mean/r005/06month/", as.character(i),
                 "/tmax_mean_r005_", as.character(i - 1), "1001_", as.character(i), "0331.nc")
  File1 <- paste0("input/covariates/raw_data/temp_max", "/tmax_mean_r005_", as.character(i - 1), "1001_", as.character(i), "0331.nc")
  URL2 <- paste0("http://opendap.bom.gov.au:8080/thredds/fileServer/agcd/tmax/mean/r005/06month/", as.character(i),
                 "/tmax_mean_r005_", as.character(i), "0401_", as.character(i), "0930.nc")
  File2 <- paste0("input/covariates/raw_data/temp_max", "/tmax_mean_r005_", as.character(i), "0401_", as.character(i), "0930.nc")

  # download files
  if (!file.exists(File1)) {
    download.file(URL1, File1, mode = "wb")
  }
  if (!file.exists(File2)) {
    download.file(URL2, File2, mode = "wb")
  }
}

# process downloaded data

# precipitation

# process data - taking the mid month of each 6 month epoch as the date of processed data
# get all .nc files in the relevant directory
PrecFiles <- list.files("input/covariates/raw_data/precipitation", pattern = "\\.nc$")

# get each nc file, project, clip, and save
for (i in 1:length(PrecFiles)) {
  if (substr(sub(".*r005_", "", PrecFiles[i]), 5, 6) == "10") {
      Date <- paste0(as.character(as.numeric(substr(sub(".*r005_", "", PrecFiles[i]), 1, 4))), "12")
  }
  if (substr(sub(".*r005_", "", PrecFiles[i]), 5, 6) == "04") {
    Date <- paste0(as.character(as.numeric(substr(sub(".*r005_", "", PrecFiles[i]), 1, 4))), "06")
  }
  if (!file.exists(paste0("input/covariates/output/hcpre", Date, ".tif"))) {
    Rast <- rast(paste0(getwd(), "/input/covariates/raw_data/precipitation/", PrecFiles[i])) %>%
                project(rast(Bounding, resolution = 5000))
    writeRaster(Rast, paste0("input/covariates/output/hcpre", Date, ".tif"), overwrite = TRUE)
  }
}

# mean temperature

# process data - taking the mid month of each 6 month epoch as the date of processed data
# get all .nc files in the relevant directory
TmeanFiles <- list.files("input/covariates/raw_data/temp_mean", pattern = "\\.nc$")

# get each nc file, project, clip, and save
for (i in 1:length(TmeanFiles)) {
  if (substr(sub(".*r005_", "", TmeanFiles[i]), 5, 6) == "10") {
    Date <- paste0(as.character(as.numeric(substr(sub(".*r005_", "", TmeanFiles[i]), 1, 4))), "12")
  }
  if (substr(sub(".*r005_", "", TmeanFiles[i]), 5, 6) == "04") {
    Date <- paste0(as.character(as.numeric(substr(sub(".*r005_", "", TmeanFiles[i]), 1, 4))), "06")
  }
  if (!file.exists(paste0("input/covariates/output/hctmn", Date, ".tif"))) {
    Rast <- rast(paste0(getwd(), "/input/covariates/raw_data/temp_mean/", TmeanFiles[i])) %>%
                project(rast(Bounding, resolution = 5000))
    writeRaster(Rast, paste0("input/covariates/output/hctmn", Date, ".tif"), overwrite = TRUE)
  }
}

# max temperature

# process data - taking the mid month of each 6 month epoch as the date of processed data
# get all .nc files in the relevant directory
TmaxFiles <- list.files("input/covariates/raw_data/temp_max", pattern = "\\.nc$")

# get each nc file, project, clip, and save
for (i in 1:length(TmaxFiles)) {
  if (substr(sub(".*r005_", "", TmaxFiles[i]), 5, 6) == "10") {
    Date <- paste0(as.character(as.numeric(substr(sub(".*r005_", "", TmaxFiles[i]), 1, 4))), "12")
  }
  if (substr(sub(".*r005_", "", TmaxFiles[i]), 5, 6) == "04") {
    Date <- paste0(as.character(as.numeric(substr(sub(".*r005_", "", TmaxFiles[i]), 1, 4))), "06")
  }
  if (!file.exists(paste0("input/covariates/output/hctma", Date, ".tif"))) {
    Rast <- rast(paste0(getwd(), "/input/covariates/raw_data/temp_max/", TmaxFiles[i])) %>%
                project(rast(Bounding, resolution = 5000))
    writeRaster(Rast, paste0("input/covariates/output/hctma", Date, ".tif"), overwrite = TRUE)
  }
}

# download 12 monthly climate variables from 1976 - 1996 if needed from BOM

# precipitation
for (i in 1976:1996) {
  # get download URLs and file names
  URL <- paste0("http://opendap.bom.gov.au:8080/thredds/fileServer/agcd/precip/total/r005/12month/", as.character(i),
                 "/precip_total_r005_", as.character(i), "0101_", as.character(i), "1231.nc")
  File <- paste0("input/covariates/raw_data/precipitation/12month", "/precip_total_r005_", as.character(i), "0101_", as.character(i), "1231.nc")

  # download files
  if (!file.exists(File)) {
    download.file(URL, File, mode = "wb")
  }
}

# mean mean temperature
for (i in 1976:1996) {
  # get download URLs and file names
  URL <- paste0("http://opendap.bom.gov.au:8080/thredds/fileServer/agcd/tmean/mean/r005/12month/", as.character(i),
                 "/tmean_mean_r005_", as.character(i), "0101_", as.character(i), "1231.nc")
  File <- paste0("input/covariates/raw_data/temp_mean/12month", "/precip_total_r005_", as.character(i), "0101_", as.character(i), "1231.nc")

  # download files
  if (!file.exists(File)) {
    download.file(URL, File, mode = "wb")
  }
}

# process downloaded data

# precipitation

# process data - average over the 1976 to 1996 period - and then save
if (!file.exists("input/covariates/output/hcltp.tif")) {
  PrecFiles <- list.files("input/covariates/raw_data/precipitation/12month", pattern = "\\.nc$", full.names = TRUE)
  Rast <- rast(PrecFiles) %>% project(rast(Bounding, resolution = 5000)) %>% mean()
  writeRaster(Rast, "input/covariates/output/hcltp.tif", overwrite = TRUE)
}

# mean temperature

# process data - average over the 1976 to 1996 period - and then save
if (!file.exists("input/covariates/output/hcltt.tif")) {
  TmeanFiles <- list.files("input/covariates/raw_data/temp_mean/12month", pattern = "\\.nc$", full.names = TRUE)
  Rast <- rast(TmeanFiles) %>% project(rast(Bounding, resolution = 5000)) %>% mean()
  writeRaster(Rast, "input/covariates/output/hcltt.tif", overwrite = TRUE)
}

# lot size

# get dcdb feature class names
FeatureClasses <- st_layers("input/covariates/raw_data/dcdb/cadastre_superimposed.gdb")

# loop through dcdb feature classes and rasterise to 25m resolution
for (i in 1:length(FeatureClasses$name)) {
  if (!file.exists(paste0("input/covariates/output/mask/htpls", substr(FeatureClasses$name[i], start = nchar(FeatureClasses$name[i]) - 5, stop = nchar(FeatureClasses$name[i])), ".tif"))) {
    Layer <- st_read("input/covariates/raw_data/dcdb/cadastre_superimposed.gdb", layer = FeatureClasses$name[i]) %>% vect() %>% project(Bounding) %>% crop(Bounding) %>% aggregate(by = "LOTPLAN")
    Layer$Area <- expanse(Layer, unit = "ha", transform = TRUE)
    Layer <- Layer %>% rasterize(rast(Bounding, resolution = 5), field = "Area") %>% aggregate(fact = 5)

    # save rasters into "mask" folder as used for masking areas in the model predictions
    writeRaster(Layer, paste0("input/covariates/output/mask/htpls", substr(FeatureClasses$name[i], start = nchar(FeatureClasses$name[i]) - 5, stop = nchar(FeatureClasses$name[i])), ".tif"), overwrite = TRUE)

    Layer <- Layer %>% wrap()

    # set up parallel processing cluster
    cl <- makeCluster(3, type = "PSOCK")
    registerDoParallel(cl)

    # loop through buffers and calculate quantities
    foreach (j = 1:3, .packages = c("tidyverse", "terra", "sf")) %dopar% {
      # unwrap raster
      Rast <- unwrap(Layer)
      # calculate quantities in buffers
      Buff <- focalMat(Rast, j * 1000, "circle", fillNA = TRUE)
      focal(Rast, Buff, filename = paste0("input/covariates/output/htpls", j, "km", substr(FeatureClasses$name[i], start = nchar(FeatureClasses$name[i]) - 5, stop = nchar(FeatureClasses$name[i])), ".tif"), overwrite = TRUE, na.rm = TRUE)
    }

    stopCluster(cl)
  }
}

# habitat covariates

# get preclaring koala habitat and habitat category lookup table
PreClearHabitat <- vect("input/covariates/raw_data/preclear_koala_habitat/KoalaSurveyStrata_v3_ful.shp") %>% project(Bounding)

# get habitat class lookup table
LookUp <- read_csv("input/covariates/raw_data/preclear_koala_habitat/koala_habitat_matrix_rules_LUT.csv")

# join habitat lookup table and dissolve by habitat class
PreClearHabitat <- PreClearHabitat %>% merge(LookUp, by.x = "HSM_RULE", by.y = "RULEID30", all.x = TRUE) %>% makeValid() %>% aggregate("CoreDesc")
# remove all fields except habitat class field
PreClearHabitat <- PreClearHabitat[, 1]
# remove polygons that are NA
PreClearHabitat <- subset(PreClearHabitat, !is.na(PreClearHabitat$CoreDesc))

# get remnant cover feature class names
FeatureClasses <- st_layers("input/covariates/raw_data/remnant_cover")

# loop through remnant cover feature classeses
for (i in 1:length(FeatureClasses$name)) {
  if (!file.exists(paste0("input/covariates/output/hhkha", substr(FeatureClasses$name[i], start = nchar(FeatureClasses$name[i]) - 3, stop = nchar(FeatureClasses$name[i])), "06.tif"))) {
    # load relevant layer
    RemLayer <- vect(paste0("input/covariates/raw_data/remnant_cover/", FeatureClasses$name[i], ".shp")) %>% project(Bounding) %>% crop(Bounding) %>% aggregate("Cover")
    # extract only remnant cover
    RemLayer <- RemLayer %>% subset(RemLayer$Cover == "remnant")
    # clip out habitat over remnant cover to get remnant habitat
    RemHab <- PreClearHabitat %>% crop(RemLayer)
    # erase remnant habitat to get non-remnant habitat
    NonRemHab <- PreClearHabitat %>% erase(RemLayer)

    # recode the habitat types to the following
    # 1 = Remnant high suitability core, 2 = Remnant medium suitability core, 3 = Remnant low suitability non-core, 4 = Remnant rainforest and non-habitat
    # 5 = Non-remnant high suitability core, 6 = Non-remnant medium suitability core, 7 = Non-remnant low suitability non-core, 8 = Non-remnant rainforest and non-habitat
    LookUp <- as_tibble(cbind(CoreDesc = c("High suitability core", "Medium suitability core", "Low suitability non-core", "Rainforest and non-habitat"), Type = c(1, 2, 3, 4)))
    LookUp <- LookUp %>% mutate(Type = as.numeric(Type))
    RemHab <- RemHab %>% merge(LookUp, by.x = "CoreDesc", by.y = "CoreDesc", all.x = TRUE)
    LookUp <- as_tibble(cbind(CoreDesc = c("High suitability core", "Medium suitability core", "Low suitability non-core", "Rainforest and non-habitat"), Type = c(5, 6, 7, 8)))
    LookUp <- LookUp %>% mutate(Type = as.numeric(Type))
    NonRemHab <- NonRemHab %>% merge(LookUp, by.x = "CoreDesc", by.y = "CoreDesc", all.x = TRUE)
    Habitat <- rbind(RemHab, NonRemHab)

    # convert to raster with 30m resolution and save
    HabitatR <- Habitat %>% rasterize(rast(Bounding, resolution = 30), field = "Type")
    writeRaster(HabitatR, paste0("input/covariates/output/hhkha", substr(FeatureClasses$name[i], start = nchar(FeatureClasses$name[i]) - 3, stop = nchar(FeatureClasses$name[i])), "06.tif"), overwrite = TRUE)
  }
}

# forest and woody cover

# download data - UPDATE URLs when necessary - search for "national forest and sparse woody vegetation data" to find the latest version of the data (currently version 7 - 1988 to 2022)
# download northern woody cover tile time series
if (!file.exists("input/covariates/raw_data/woody_cover/sg56_woody_8822.zip") & !file.exists("input/covariates/raw_data/woody_cover/sh56_woody_8822.zip")) {
download.file("https://data.gov.au/data/dataset/48900c21-c5e2-499a-86f0-20b1b82275e2/resource/f0c169a3-afce-413c-be9e-3c27f9825a85/download/sg56_woody_8822.zip", "input/covariates/raw_data/woody_cover/sg56_woody_8822.zip", mode = "wb")
unzip("input/covariates/raw_data/woody_cover/sg56_woody_8822.zip", exdir = "input/covariates/raw_data/woody_cover/sg56/", junkpaths = TRUE)
# download southern woody cover tile time series
download.file("https://data.gov.au/data/dataset/48900c21-c5e2-499a-86f0-20b1b82275e2/resource/6870c825-e5fe-4313-9747-9fa9bce98950/download/sh56_woody_8822.zip", "input/covariates/raw_data/woody_cover/sh56_woody_8822.zip", mode = "wb")
unzip("input/covariates/raw_data/woody_cover/sh56_woody_8822.zip", exdir = "input/covariates/raw_data/woody_cover/sh56/", junkpaths = TRUE)

# get list of .tif files
SGTiffFiles <- as_tibble(list.files("input/covariates/raw_data/woody_cover/sg56/", pattern = "\\.tif$", full.names = FALSE)) %>% arrange(value)
SHTiffFiles <- as_tibble(list.files("input/covariates/raw_data/woody_cover/sh56/", pattern = "\\.tif$", full.names = FALSE)) %>% arrange(value)

# loop through each northern tile file and process
for (i in 1:length(SGTiffFiles$value)) {
  # get year
  Year <- as.numeric(substr(SGTiffFiles$value[i], nchar(SGTiffFiles$value[i]) - 5, nchar(SGTiffFiles$value[i]) - 4))
  if (Year < 88) {
    Year <- Year + 2000
  } else {
    Year <- Year + 1900
  }

  if (!file.exists(paste0("input/covariates/output/hhfwc", Year, "06.tif"))) {
    # check if there is a matching file in the southern woody cover data for this year
    if (any(SHTiffFiles$value == SGTiffFiles$value[i])) {
      # load both tiff files for this year
      SG <- rast(paste0("input/covariates/raw_data/woody_cover/sg56/", SGTiffFiles$value[i]))
      SH <- rast(paste0("input/covariates/raw_data/woody_cover/sh56/", SHTiffFiles$value[i]))
      # mosaic them together, project, clip, and reclassify
      Woody <- mosaic(SG, SH, fun = "first") %>% project(rast(Bounding, resolution = 25)) %>% crop(Bounding) %>% as.int() %>% classify(cbind(c(0, 1, 2), c(3, 2, 1)))

      # save
      writeRaster(Woody, paste0("input/covariates/output/hhfwc", Year, "06.tif"), overwrite = TRUE)
     }
    }
  }
}

# groundwater dependent ecosystems

# get ground water dependent ecosystems feature class, project and then clip to study area
if (!file.exists("input/covariates/output/hhgde.tif")) {
  GDEVect <- st_read("input/covariates/raw_data/gde/groundwater_GDE.gdb", layer = "gw_gde_hgpot") %>% vect() %>% project(Bounding) %>% crop(Bounding) %>% aggregate("HGPOT7_DESC")
  # create raster and classify
  # 1 = Brackish, saline or fluctuating salinity with intermittent connectivity, 2 = Intermittent and freshwater, 3 = Near-permanent and freshwater, 4 = Permanent and fluctuating salinity, 5 = Permanent and freshwater, 6 = Exclusion and recharge zones
  GDERast <- GDEVect %>% rasterize(rast(Bounding, resolution = 30), field = "HGPOT7_DESC")
  GDERast <- GDERast %>% classify(cbind(levels(GDERast)[[1]]$ID, c(1, 6, 2, 3, 4, 5, 6)))

  # save
  writeRaster(GDERast, "input/covariates/output/hhgde.tif", overwrite = TRUE)
}

# canopy height

# get canopy height data
if (!file.exists("input/covariates/output/hhcht.tif")) {
  CHRast <- rast("https://dap.tern.org.au/thredds/fileServer/landscapes/remote_sensing/spatial_other/jrsrp/height/alpsbk_aust_y2009_se4a2.tif") %>% project(rast(Bounding, resolution = 30))

  # save
  writeRaster(CHRast, "input/covariates/output/hhcht.tif", overwrite = TRUE)
}

# understory fraction

# get understorey fraction data
if (!file.exists("input/covariates/output/hhunf.tif")) {
  UDRast <- rast("https://dap.tern.org.au/thredds/fileServer/landscapes/remote_sensing/spatial_other/jrsrp/height/alpsbk_aust_y2009_se1a2.tif") %>% project(rast(Bounding, resolution = 30))

  # save
  writeRaster(UDRast, "input/covariates/output/hhunf.tif", overwrite = TRUE)
}

#cat("\n","\n",
#    "################################################################",
#    "################  THIS CODE HAS FINISHED  ######################",
#    "################################################################", sep = "\n")

# You can now close this script by clicking the X next to its name in the script tab.
