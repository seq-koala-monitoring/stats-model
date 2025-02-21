# -------------------------------------------------------------------
#           DATA PRE-PROCESSING REQUIRED FOR THE MODEL
# -------------------------------------------------------------------

# This R code process all data to the format required for the Bayesian state-space model 
# used to estimate koala densities across Southeast Queensland.
#
# The basic workflow involves three main steps:
# 1) Updating the koala survey database
# 2) Extracting covariate values for the transects in the koala survey database
# 3) Preparing the data for the modelling phase
#
# Before anything, start fresh by clicking Session > Terminate R... > Yes in the pop-up (do not save anything if asked)
#
# To get started, please, select all lines by pressing Ctrl + A on a Windows PC or Command + A on a MAc. 
# Then, run these lines by pressing Ctrl + Enter on a Windows PC or Command + Return on a Mac.
#
# When you see "THIS CODE HAS FINISHED" in the Console panel (usually at the bottom left),
# you're ready to the modelling runs. Go to File > Open file... > locate the file named model_runs.R > Open.
# This file will run the Bayesian state-space models.
#
# PS: 1) Disregard any warnings on the task bar about packages that are not installed 
#     2) This code may take anywhere from a few minutes to days to run, 
#        depending on how many files need updating and your computer's specifications.

# -------------------------------------------------------------------
# install packages if required
source("code/install_packages_data_processing.R")

# load libraries
library(tidyverse)
library(terra)
library(exactextractr)
library(sf)
library(tidyterra)
library(abind)
library(mice)
library(factoextra)
library(nimble)
library(devtools)
if (!require("SEQKoalaDataPipeline")) devtools::install_github('seq-koala-monitoring/data-pipeline')
library(SEQKoalaDataPipeline)
library(readr)
library(rstudioapi)

# load parameters
source("parameters_init.txt")

# read utility functions
source("code/functions.R")

# Set secondary grid size
secondary_grid_size <- primary_grid_size * secondary_grid_multiple

# set genetic population units file path
gen_pop_file_path <- list.files(
  pattern = genetic_units,
  recursive = T,
  ignore.case = T,
  full.names = F,
  include.dirs = T
)
gen_pop_file_path <- gen_pop_file_path[
  (substring(gen_pop_file_path, nchar(gen_pop_file_path) - 3 + 1) == "shp")]
if (length(gen_pop_file_path) == 0) {
  stop(sprintf("No % s file found. Please, make sure you have a file named % s in this working directory",
               genetic_units, genetic_units))} # close {} of error message
if(length(gen_pop_file_path) > 1){
  message("Multiple files found:")
  for (i in seq_along(gen_pop_file_path)) {
    cat(i, ":", gen_pop_file_path[i], "\n")}
  choice <- as.integer(readline("Enter the number corresponding to the file you want to use: "))
  if (is.na(choice) || choice < 1 || choice > length(gen_pop_file_path)) {
    stop("Invalid selection. Process terminated. Please, rerun this function and choose one of the provided options.")}
  gen_pop_file_path <- gen_pop_file_path[choice]}

# Update database if required
if(update_database) {
  fcn_update_db()
}

# Set output path
state <- fcn_get_state()

# Get current date
current_date <- format(Sys.Date(), format="%Y%m%d")

# Prepare integration of survey data and covariate data

# Set directories
working_data_dir <- paste0(getwd(), "/input")
target_dir <- paste0(getwd(), "/input/survey_data")
fcn_set_home_dir(working_data_dir) # Home directory

## Set db path
fcn_set_db_path(list(
  `1996` = 'databases/SEQkoalaData.accdb',
  `2015` = 'databases/2015-2019 SEQKoalaDatabase DES_20231027.accdb',
  `2020` = 'databases/KoalaSurveyData2020_cur.accdb',
  `integrated` = ifelse(update_database == T,
                        'databases/Integrated_SEQKoalaDatabase_updated.accdb',
                        'databases/Integrated_SEQKoalaDatabase.accdb')
))

# Set gdb path
fcn_set_gdb_path(list(
  koala_survey_data = "KoalaSurveyData.gdb",
  total_db = ifelse(update_database == T,
                    "transects_spatial_representation/Integrated_SEQKoalaDatabase_Spatial_updated.shp",
                    "transects_spatial_representation/Integrated_SEQKoalaDatabase_Spatial.shp"),
  koala_survey_sites = "survey_sites/KoalaSurveySites_231108.shp"
))

# Set grid size
fcn_set_grid_size(grid_size = primary_grid_size)

# Set line transect buffer width in meters (if generating transects using start and end coordinate information)
fcn_set_line_transect_buffer(line_transect_buffer)

# Set covariate impute buffer distance
fcn_set_cov_impute_buffer(cov_impute_buffer)

# Set study area buffer
fcn_set_study_area_buffer(area_buffer)

# Conduct imputaiton if needed
if (use_imputation) {
  fcn_set_raster_path(list(covariates = 'covariates_impute/output'))
} else {
  fcn_set_raster_path(list(covariates = 'covariates/output'))
}

# Set output directories
out_dir <- target_dir     #paste0(target_dir, '/', paste0(current_date, "_", state$grid_size, ifelse(use_imputation, '_1500', '')))
if (!dir.exists(out_dir)) dir.create(out_dir)
if (!dir.exists(paste0(out_dir, '/cov_raster'))) dir.create(paste0(out_dir, '/cov_raster'))
if (!dir.exists(paste0(out_dir, '/cov_csv'))) dir.create(paste0(out_dir, '/cov_csv'))

# Write to grid
grid_raster <- fcn_get_grid()
terra::writeRaster(grid_raster, paste0(out_dir, "\\grid_raster.tif"), overwrite = T)
grid_vector <- terra::as.polygons(grid_raster)
terra::writeVector(grid_vector, paste0(out_dir, "\\grid_vec.shp"), overwrite = T)

# Load the survey data as tables
master <- fcn_all_tables_detect()
saveRDS(master, paste0(out_dir, '/master.rds'))
master_sf <- fcn_all_tables_sf_detect()
lapply(seq_along(master_sf), \(i) sf::st_write(master_sf[[i]], paste0(out_dir, '/master_', names(master_sf)[i], '.shp'), append=F))

# Include mean max temperature and mean total precipitation to each transect
master.sub <- master[c("line_transect", "strip_transect", "uaoa")]
download_temp_precip(master.sub)
master.sub <- extract_temp_precip(master.sub)
master <- append(master.sub, master["perp_distance"])
saveRDS(master, paste0(out_dir, '/master.rds'))

# Load covariates from the directory

# Extract covariates
if (run_cov_extraction) {
  dates <- fcn_get_date_intervals()
  cov_constant_array <- fcn_cov_array('constant', write_path = out_dir)
  fcn_cov_array_detect()
  source('code/cov_temporal_array.R')
}

# Extract and save only those in surveylocations as a separate file
# Read results back from disk (output folder)
cov_constant_array <- readr::read_rds(paste0(out_dir, "/cov_constant_array.rds"))
cov_temporal_array <- readr::read_rds(paste0(out_dir, "/cov_temporal_array.rds"))
if (use_imputation) cov_temporal_array <- fcn_impute_temporal_cov(cov_temporal_array)

# Load grid fractions as tables
grid_fractions <- fcn_all_transect_grid_fractions_detect()
grid_fractions_comb <- dplyr::bind_rows(grid_fractions, .id = 'transect')
readr::write_rds(grid_fractions_comb, paste0(out_dir, '/grid_fractions.rds'))
data.table::fwrite(grid_fractions_comb, paste0(out_dir, "/grid_fractions.csv"))

# Save grid cells in survey locations separately
cov_constant_array_surveylocations <- cov_constant_array[cov_constant_array[,1] %in% grid_fractions_comb$GridID,,]
cov_temporal_array_surveylocations <- cov_temporal_array[cov_temporal_array[,1,1] %in% grid_fractions_comb$GridID,,]
readr::write_rds(cov_constant_array_surveylocations, paste0(out_dir, "/cov_constant_array_surveylocations.rds"))
readr::write_rds(cov_temporal_array_surveylocations, paste0(out_dir, '/cov_temporal_array_surveylocations.rds'))

# Resave grid fractions for the subset with complete covariate information
grid_id_non_na <- fcn_complete_grid_id(cov_constant_array, cov_temporal_array) # Get GridID index with non-NA covariates
cov_constant_array$GridID[!(cov_constant_array$GridID %in% grid_id_non_na)]
grid_fractions_complete_cov <- fcn_all_transect_grid_fractions_detect(grid_id_vec = grid_id_non_na)
grid_fractions_comb_complete_cov <- dplyr::bind_rows(grid_fractions_complete_cov, .id = 'transect')
readr::write_rds(grid_fractions_comb_complete_cov, paste0(out_dir, '/grid_fractions_complete_cov.rds'))
data.table::fwrite(grid_fractions_comb_complete_cov, paste0(out_dir, "/grid_fractions_complete_cov.csv"))
transects_missing <- unique(grid_fractions_comb$TransectID)[(!(unique(grid_fractions_comb$TransectID) %in% unique(grid_fractions_comb_complete_cov$TransectID)))]
print(transects_missing)
GridID_missing <- grid_fractions_comb[grid_fractions_comb$TransectID %in% transects_missing,]$GridID

# Save grid cells in survey locations separately with non-NA covariate information
cov_constant_array_surveylocations <- cov_constant_array[cov_constant_array[,1] %in% grid_fractions_comb_complete_cov$GridID ,,]
cov_temporal_array_surveylocations <- cov_temporal_array[cov_temporal_array[,1,1] %in% grid_fractions_comb_complete_cov$GridID,,]

# Save full arrays
cov_constant_array[cov_constant_array[,1] %in% grid_id_non_na,,] %>%
  readr::write_rds(paste0(out_dir, '/cov_constant_array_complete_cov.rds'))
cov_temporal_array[cov_temporal_array[,1,1] %in% grid_id_non_na,,] %>%
  readr::write_rds(paste0(out_dir, '/cov_temporal_array_complete_cov.rds'))

# Check for complete cases (any missing values)
fcn_complete_cases_check(cov_constant_array_surveylocations)
fcn_complete_cases_check(cov_temporal_array_surveylocations)
readr::write_rds(cov_constant_array_surveylocations, paste0(out_dir, "/cov_constant_array_complete_cov_surveylocations.rds"))
readr::write_rds(cov_temporal_array_surveylocations, paste0(out_dir, '/cov_temporal_array_complete_cov_surveylocations.rds'))

# Produce date interval lookup table
write.csv(fcn_date_interval_lookup(), paste0(out_dir, "/date_interval_lookup.csv"))
cov_layer_df <- fcn_covariate_layer_df_detect()
write.csv(cov_layer_df[,1:5], paste0(out_dir, '/covariate_info.csv'))

# free up some memory
gc()

# Produce and save the adjacency matrix
adj_data <- fcn_adj_matrix(secondary_grid_size = secondary_grid_size)
saveRDS(adj_data, paste0(out_dir, "/adj_data_queen.rds"))
terra::writeRaster(adj_data$grid_raster_sp, paste0(out_dir, "/grid_raster_secondary.tif"), overwrite = T)
grid_vec_sp <- terra::as.polygons(adj_data$grid_raster_sp)
terra::writeVector(grid_vec_sp, paste0(out_dir, "/grid_vec_sp.shp"), overwrite=T)
adj_data <- fcn_adj_matrix(directions = 'rook')
saveRDS(adj_data, paste0(out_dir, "/adj_data_rook.rds"))

# Write lookup table of GridID to genetic populations
gen_pop_file <- sf::st_read(gen_pop_file_path)
gen_pop_lookup <- fcn_grid_intersect_feature(gen_pop_file, field = gen_pop_column_id)
saveRDS(gen_pop_lookup, paste0(out_dir, "/gen_pop_lookup.rds"))

# end
print("Data compilation step complete")

# load input data
Surveys <- readRDS("input/survey_data/master.rds")
GridFrac <- readRDS("input/survey_data/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data/cov_temporal_array_surveylocations.rds")
DateIntervals <- read_csv("input/survey_data/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))

# format data and check for multi-collinearity

# loop through order, lag, and vartrend values
for (Order in 1:2) {
  for(Lag in 0:2) {
    for (VarTrend in 0:1) {
      # set seed
      set.seed(seed)
      
      # generate formatted data
      FormattedData <- format_data(Surveys = Surveys, GridFrac = GridFrac, CovConsSurv = CovConsSurv, CovTempSurv = CovTempSurv, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup, AggGenPop = gen_pop_agg, FirstDate = FirstDate, LastDate = LastDate, Order = Order, Lag = Lag, VarTrend = VarTrend)
     
      # save data
      saveRDS(FormattedData, paste0("input/nimble_data/format_data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))
      
      # save continuous variable correlations
      write.csv(FormattedData$CorrXY, paste0("input/nimble_data/correlations/cor_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".csv"))
      
      # save PCA plots and PCAs
      saveRDS(FormattedData$SoilPCA, paste0("input/nimble_data/pca/pca_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))
      ggsave(FormattedData$SoilScree1, file = paste0("input/nimble_data/pca/soilscree1_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".jpg"), width = 20, height = 20, units = "cm", dpi = 300)
      ggsave(FormattedData$SoilScree2, file = paste0("input/nimble_data/pca/soilscree2_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".jpg"), width = 20, height = 20, units = "cm", dpi = 300)
      ggsave(FormattedData$SoilBiPlot, file = paste0("input/nimble_data/pca/soilbiplot_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".jpg"), width = 20, height = 20, units = "cm", dpi = 300)

      # free memory
      rm(FormattedData)
      gc()
    }
  }
}

# create data for nimble

# loop through order, lag, and vartrend values
for (Order in 1:2) {
  for(Lag in 0:2) {
    for (VarTrend in 0:1) {
      # set seed
      set.seed(seed)

      # load formatted data
      # save data
      FormattedData <- readRDS(paste0("input/nimble_data/format_data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

      # generate data to fit models
      FitData <- get_fit_data(Data = FormattedData, StaticVars = static_variables, DynamicVars = dynamic_variables, ObsVars = obs_variables)
      
      # save data
      saveRDS(FitData, paste0("input/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

      # free memory
      rm(FitData)
      gc()
    }
  }
}

# create mask matrix

# get date intervals
DateIntervals <- read_csv("input/survey_data/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))

# get the mask for the predictions
# load grids
grid.rast <- rast("input/survey_data/grid_raster.tif")
grid.vect <- vect("input/survey_data/grid_vec.shp") %>% st_as_sf()
# load land uses
lu1999 <- rast("input/covariates/output/mask/lu1999_mask.tif")
lu2017 <- rast("input/covariates/output/mask/lu2017_mask.tif")
# get zonal statisrics for grid cells
lu1999.zonal <- exact_extract(lu1999, grid.vect, "mode") %>% as.data.frame()
lu2017.zonal <- exact_extract(lu2017, grid.vect, "mode") %>% as.data.frame()
names(lu1999.zonal) <- "Tertiary"
names(lu2017.zonal) <- "Tertiary"

lu1999.res <- resample(rast("input/covariates/output/mask/lu1999_mask.tif"), grid.rast, method = "near")
lu2017.res <- resample(rast("input/covariates/output/mask/lu2017_mask.tif"), grid.rast, method = "near")
# load lot sizes and created a named list
files <- list.files("input/covariates/output/mask", pattern = "htpls", full.names = T)
names <- list.files("input/covariates/output/mask", pattern = "htpls", full.names = F)
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
lookup1999 <- read_csv("input/covariates/output/mask/land_use_lookup1999_mask.csv")
lookup2017 <- read_csv("input/covariates/output/mask/land_use_lookup2017_mask.csv")
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
# create a vector to represent each TimePeriodID
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
saveRDS(lu.mask.matrix, "input/covariates/output/mask/lu_mask_matrix.rds")

# end
print("Data formatting for nimble complete")
