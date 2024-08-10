# load libraries
library(tidyverse)
library(abind)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)
library(MCMCvis)
library(terra)
library(tidyterra)
library(foreach)
library(doParallel)
library(exactextractr)

# get functions
source("functions.r")

# get input data
CovCons <- readRDS("input/survey_data_500/cov_constant_array.rds")
CovTemp <- readRDS("input/survey_data_500/cov_temporal_array.rds")
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
PredData <- list(CovCons = CovCons, CovTemp = CovTemp, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup)
rm(CovCons, CovTemp, DateIntervals, GenPopLookup)

# PREDICTIONS

# ORDER = 1, LAG = 0, VARTREND = 0

# set model structure parameters and time scale parameters of the fitted model
Order = 1
Lag = 0
VarTrend = 1
FirstDate <- "1996-02-13"

# load nimble data and mcmc data
NimbleData <- readRDS(paste0("output/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))
MCMC <- readRDS(paste0("output/mcmc/sel_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))
# extract MCMC chanins and stitch MCMC chains together
MCMC <- rbind(MCMC$MCMC[[1]], MCMC$MCMC[[2]], MCMC$MCMC[[3]])

# sub sample MCMC chains if needed
MCMC <- MCMC[sample(1:dim(MCMC)[1], 1000), ]

# load feature class of small grids
Grid <- vect("input/survey_data_500/grid_vec.shp")

# make cluster
#cl <- makeCluster(3)
#registerDoParallel(cl)

# export packages to cluster
#clusterEvalQ(cl,{
#  library(tidyverse)
#  library(terra)
#  library(tidyterra)
#  library(abind)
#  })

# loop through years to generate predictions
Predictions <- foreach(i = 1996:2023) %do% {

  # format data for predictions
  Data <- get_prediction_data(Year = i, NimbleData, PredData, RainForestMask = TRUE)

  # generate predictions
  Preds <- get_predictions(MCMC, Data)

  # set NAs to -9999 for shapefile to reflect masked grid cells
  PredsShp <- Preds$Spatial %>% mutate(Expected = ifelse((Expected == 0) & (LowerCI == 0) & (UpperCI == 0), -9999, Expected), LowerCI = ifelse((Expected == 0) & (LowerCI == 0) & (UpperCI == 0), -9999, LowerCI), UpperCI = ifelse((Expected == 0) & (LowerCI == 0) & (UpperCI == 0), -9999, UpperCI), SD = ifelse((Expected == 0) & (LowerCI == 0) & (UpperCI == 0), -9999, SD))

  # merge predictions with spatial grid
  PredsGrid <- Grid %>% merge(PredsShp, all.x = TRUE, by.x = "GridID", by.y = "GridID")

  # remove grids outside of genetic populations or with no predictions due to missing data
  PredsGrid <- PredsGrid %>% filter(!is.na(Expected))

  # save vector
  writeVector(PredsGrid, paste0("output/predictions/predictions", i, ".shp"), overwrite = TRUE)

  Preds
}

#stop cluster
stopCluster(cl)


