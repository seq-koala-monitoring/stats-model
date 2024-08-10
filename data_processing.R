# load libraries
library(tidyverse)
library(terra)
library(exactextractr)
library(sf)
library(tidyterra)

# create mask matrix

# get date intervals
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))

# get the mask for the predictions
# load grids
grid.rast <- rast("input/survey_data_500/grid_raster.tif")
grid.vect <- vect("input/survey_data_500/grid_vec.shp") %>% st_as_sf()
# load land uses
lu1999 <- rast("input/mask/lu1999_mask.tif")
lu2017 <- rast("input/mask/lu2017_mask.tif")
# get zonal statisrics for grid cells
lu1999.zonal <- exact_extract(lu1999, grid.vect, "mode") %>% as.data.frame()
lu2017.zonal <- exact_extract(lu2017, grid.vect, "mode") %>% as.data.frame()
names(lu1999.zonal) <- "Tertiary"
names(lu2017.zonal) <- "Tertiary"

lu1999.res <- resample(rast("input/mask/lu1999_mask.tif"), grid.rast, method = "near")
lu2017.res <- resample(rast("input/mask/lu2017_mask.tif"), grid.rast, method = "near")
# load lot sizes and created a named list
files <- list.files("input/mask", pattern = "htpls", full.names = T)
names <- list.files("input/mask", pattern = "htpls", full.names = F)
names <- as.character(parse_number(names))
htpls <- lapply(files, rast)
names(htpls) <- names
# summarize raster values over grid polygons (it takes some minutes to run)
htpls.grid <- lapply(htpls, function(r) {
  # mean lot size per grid cell
  g <- exact_extract(r, grid.vect, "mean")
  # rename it
  names(g) <- "mean_lot_size"
  # extract mean lot size per GridID
  gid <- grid.vect %>% cbind(g) %>% vect()
  # rename
  names(gid) <- c("GridID", "mean_lot_size")
  return(gid)})
# add raster date to the previous list
for(i in 1:length(htpls.grid)){
  htpls.grid[[i]]$date <- names(htpls.grid)[i]}
# Load look up tables to mask cells by land use only
lookup1999 <- read_csv("input/mask/land_use_lookup1999_mask.csv")
lookup2017 <- read_csv("input/mask/land_use_lookup2017_mask.csv")
# use the median lot size (ha) in Brisbane CBD and inner suburbs (2021) as threshold
# this value represents highly urbanized areas where koalas are unlikely to occur
median.lot.size.brisbane <- 0.63
# create temporal mask dataframes (it takes some minutes to run)
htpmask <- lapply(htpls.grid, function(r){
  # set reference dates
  date.1999 <- as.Date("1999-06-01")
  date.2017 <- as.Date("2017-06-01")
  # set lot size date
  date.vect <- as.Date(paste0(str_sub(unique(r$date), end = -3), "-", str_sub(unique(r$date), start = 5), "-", "01"))
  # Use either land use raster from 1999 or 2017 based on the closest date
  if(date.vect <= date.1999 |
     abs(as.numeric(date.1999 - date.vect)) <= abs(as.numeric(date.2017 - date.vect))){
    j <- cbind(as.data.frame(r), as.data.frame(lu1999.zonal)) %>% left_join(lookup1999, join_by("Tertiary" == "code")) %>% mutate(Mask = case_when(Tertiary.y %in% c("Urban residential", "Public services", "Recreation and culture") & mean_lot_size <= median.lot.size.brisbane ~ 1, mask == 1 ~ 1, .default = 0)) %>% select(GridID, Mask, date)
  } else {
    j <- cbind(as.data.frame(r), as.data.frame(lu2017.zonal)) %>% left_join(lookup2017, join_by("Tertiary" == "code")) %>% mutate(Mask = case_when(Tertiary.y %in% c("Urban residential", "Public services", "Recreation and culture") & mean_lot_size <= median.lot.size.brisbane ~ 1, mask == 1 ~ 1, .default = 0)) %>% select(GridID, Mask, date)
  }
  return(j)})

# assign a TimePeriodID to each temporal dataframe
# store TimePeriodIds
TimePeriod <- vector()
for(i in 1:length(htpmask)){
  name <- htpmask[[i]][1, "date"]
  name <- as.Date(paste0(str_sub(parse_number(name), end = -3), "-", str_sub(parse_number(name), start = 5), "-", "01"))
  sub <- DateIntervals[which(DateIntervals$start_date <= name & DateIntervals$end_date >= name),]
  htpmask[[i]]$TimePeriodID <- sub$TimePeriodID
  TimePeriod[i] <- as.character(sub$TimePeriodID)}
# name the previous list
names(htpmask) <- TimePeriod
# create an empty list to store a mask dataframe per TimePeriodID
lu.mask.list <- vector(mode = "list", length = length(unique(DateIntervals$TimePeriodID)))
names(lu.mask.list) <- as.character(unique(DateIntervals$TimePeriodID))
# create a vector to represent each TimePeriodID (from 1 to 60)
index <- as.numeric(names(htpmask))
# append each temporal mask dataframe to the corresponding slot in the empty list
# slots where the TimePeriodID does not have an exact match are filled with data from the closest matching TimePeriodID in the mask dataframe
for(i in 1:length(lu.mask.list)) {
  n <- which.min(abs(index - i))
  lu.mask.list[[i]] <- htpmask[[n]] %>% select(-date) %>% mutate(TimePeriodID = as.character(i))}
# create a matrix with GridID as rows, TimePeriodID as columns, and values as 0 (do not mask out) or 1 (mask out)
lu.mask.matrix <- lapply(lu.mask.list, function(df){
  df <- df %>% pivot_wider(values_from = Mask, names_from = TimePeriodID)}) %>% reduce(left_join, by = "GridID") %>% column_to_rownames("GridID") %>% as.matrix()
# remove temporary objects and free up memory
rm(i, n, index, files, names,lookup1999, lookup2017, htpls, lu1999.res, lu2017.res, htpls.grid, sub, lu.mask.list, TimePeriod, name, median.lot.size.brisbane, grid.rast, grid.vect)
gc()

# save the mask matrix
saveRDS(lu.mask.matrix, "input/mask/lu_mask_matrix.rds")
