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

# get functions
source("functions.r")

# load input data
Surveys <- readRDS("input/survey_data_500/master.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
CovCons <- readRDS("input/survey_data_500/cov_constant_array.rds")
CovTemp <- readRDS("input/survey_data_500/cov_temporal_array.rds")
PredData <- list(CovCons = CovCons, CovTemp = CovTemp, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup)
rm(CovCons, CovTemp, DateIntervals, GenPopLookup)

# set up data frame to store model selection results
ModelWAICs <- tibble(Order = rep(1, 6), Lag = rep(NA, 6), VarTrend = rep(NA, 6), WAIC = rep(NA, 6))

# loop through models
for (Lag in c(0, 1, 2)) {
  for (VarTrend in c(0,1)) {

    # get model
    Model <- readRDS(paste0("output/mcmc/sat_order1_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))
    ModelWAICs[Lag * 2 + VarTrend + 1, "Lag"] <- Lag
    ModelWAICs[Lag * 2 + VarTrend + 1, "VarTrend"] <- VarTrend
    ModelWAICs[Lag * 2 + VarTrend + 1, "WAIC"] <- lapply(Model$WAIC, FUN = function(x) {return(x$WAIC)}) %>% unlist() %>% mean()
  }
}

# get the best model
BestIndex <- which(ModelWAICs$WAIC == min(ModelWAICs$WAIC))
BestModelSat <- readRDS(paste0("output/mcmc/sat_order1_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".rds"))
BestModelSel <- readRDS(paste0("output/mcmc/sel_order1_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".rds"))

# check for convergence
MCMCsummary(BestModelSat$MCMC)
MCMCtrace(BestModelSat$MCMC, filenane = "output/mcmc/figures/best_sat_trace.jpg")
MCMCsummary(BestModelSel$MCMC)
MCMCtrace(BestModelSel$MCMC, filenane = "output/mcmc/figures/best_sel_trace.jpg"))

# do posterior predictive checks





# generate parameter estimates table



# generate predictions

# extract MCMC chains and stitch MCMC chains together
MCMC <- rbind(BestModelSel$MCMC[[1]], BestModelSat$MCMC[[2]], BestModelSat$MCMC[[3]])

# sub sample MCMC chains if needed
MCMC <- MCMC[sample(1:dim(MCMC)[1], 10000), ]

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

# set up list to store predictions
PredsList <- list()

# loop through years to generate predictions
Predictions <- foreach(i = year(FirstDate):year(LastDate)) %do% {

  # format data for predictions
  Data <- get_prediction_data(Year = i, BestModelSel$Data, PredData, RainForestMask = TRUE)

  # generate predictions
  Preds <- get_predictions(MCMC, Data)

  # set masked grid cells to NA
  PredsShp <- Preds$Spatial %>% mutate(Expected = ifelse((Expected == 0) & (LowerCI == 0) & (UpperCI == 0), NA, Expected), LowerCI = ifelse((Expected == 0) & (LowerCI == 0) & (UpperCI == 0), NA, LowerCI), UpperCI = ifelse((Expected == 0) & (LowerCI == 0) & (UpperCI == 0), NA, UpperCI), SD = ifelse((Expected == 0) & (LowerCI == 0) & (UpperCI == 0), NA, SD))

  # merge predictions with spatial grid
  PredsGrid <- Grid %>% merge(PredsShp, all.x = TRUE, by.x = "GridID", by.y = "GridID")

  # remove NA values
  PredsGrid <- PredsGrid %>% filter(!is.na(Expected))

  # save vector
  writeVector(PredsGrid, paste0("output/predictions/sel_order1_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, "_year", i, ".shp"), overwrite = TRUE)

  PredsList <- append(PredsList, list(Preds))
}

# save predictions
saveRDS(PredsList, paste0("output/predictions/sel_order1_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".rds"))

Years <- year(FirstDate):year(LastDate)
Expected <- lapply(PredsList, FUN = function(x){return(x$TotalSumMean)}) %>% unlist()
Lower <- lapply(PredsList, FUN = function(x){return(x$TotalSumLower)}) %>% unlist()
Upper <- lapply(PredsList, FUN = function(x){return(x$TotalSumUpper)}) %>% unlist()
Abundances <- tibble(Year = Years, Expected = Expected, Lower = Lower, Upper = Upper)

Plot1 <- ggplot(Abundances, aes(x = Year, y = Expected)) + geom_line() + geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1, color = "red") + labs(x = "Year", y = "Abundance") + theme_minimal() + theme(legend.position = "none", axis.text = element_text(size = 16),  axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18, vjust = -1)) + scale_y_continuous(labels = scales::comma)

ggsave(Plot1, file = "figures/trend.jpg", width = 20, height = 15, units = "cm", dpi = 300)
