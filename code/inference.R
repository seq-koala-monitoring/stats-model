# THIS SCRIPT DOES THE INFERENCE BY SELECTING THE BEST MODEL, CHECKS MODEL ADEQUACY, AND GENERATES PREDICTIONS 

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
library(patchwork)

# load parameters
source("parameters_init.txt")

# get functions
source("code/functions.r")

# read nimble code
source("code/nimble_code.R")

# load input data
Surveys <- readRDS("input/survey_data/master.rds")
DateIntervals <- read_csv("input/survey_data/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
CovCons <- readRDS("input/survey_data/cov_constant_array.rds")
CovTemp <- readRDS("input/survey_data/cov_temporal_array.rds")
Mask <- readRDS("input/mask/lu_mask_matrix.rds")
PredData <- list(CovCons = CovCons, CovTemp = CovTemp, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup, Mask = Mask)
rm(CovCons, CovTemp, DateIntervals, GenPopLookup, Mask)
gc()

# set orders, lags, and vartrends want to consider
Orders <- c(1)
Lags <- c(0, 1, 2)
VarTrends <- c(0, 1)

# set up data frame to store model selection results
ModelWAICs <- tibble(Order = rep(NA, length(Orders) * length(Lags) * length(VarTrends)), Lag = rep(NA, length(Orders) * length(Lags) * length(VarTrends)), VarTrend = rep(NA, length(Orders) * length(Lags) * length(VarTrends)), WAIC = rep(NA, length(Orders) * length(Lags) * length(VarTrends)))

# loop through models
for (Order in Orders) {
  for (Lag in Lags) {
    for (VarTrend in VarTrends) {
      # get model
      Model <- readRDS(paste0("output/mcmc/sat_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))
      ModelWAICs[(Order - 1) * length(Orders) * length(Lags) + Lag * length(VarTrends) + VarTrend + 1, "Order"] <- Order
      ModelWAICs[(Order - 1) * length(Orders) * length(Lags) + Lag * length(VarTrends) + VarTrend + 1, "Lag"] <- Lag
      ModelWAICs[(Order - 1) * length(Orders) * length(Lags) + Lag * length(VarTrends) + VarTrend + 1, "VarTrend"] <- VarTrend
      ModelWAICs[(Order - 1) * length(Orders) * length(Lags) + Lag * length(VarTrends) + VarTrend + 1, "WAIC"] <- lapply(Model$WAIC, FUN = function(x) {return(x$WAIC)}) %>% unlist() %>% mean()
    }
  }
}

# write table
write_csv(ModelWAICs, paste0("output/inference/waics", "_firstdate", FirstDate, ".csv"))

# get the best model
BestIndex <- which(ModelWAICs$WAIC == min(ModelWAICs$WAIC))
BestModelSat <- readRDS(paste0("output/mcmc/sat_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".rds"))
BestModelSel <- readRDS(paste0("output/mcmc/sel_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".rds"))

# check for convergence
write_csv(MCMCsummary(BestModelSat$MCMC), paste0("output/assessment/sat_convergence_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".csv"))
MCMCtrace(BestModelSat$MCMC, filename = paste0("output/assessment/sat_trace_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".pdf"))
write_csv(MCMCsummary(BestModelSel$MCMC), paste0("output/assessment/sel_convergence_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".csv"))
MCMCtrace(BestModelSel$MCMC, filename = paste0("output/assessment/sel_trace_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".pdf"))

# do posterior predictive checks for the best model

# get the data for fitting models
FitData <- readRDS(paste0("input/nimble_data/data_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".rds"))

# fit the model for posterior predictive checks

# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(coda)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = c(seed, seed + 20, seed + 40), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("Strip_res_obs", "Strip_res_sim", "Strip_var", "AoA_res_obs", "AoA_res_sim", "AoA_var", "Line_res_obs", "Line_res_sim", "Line_var") , Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model_check)

#stop cluster
stopCluster(cl)

# combine samples and save
PostCheck <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)

# save posterior predictive checks model
saveRDS(PostCheck, paste0("output/assessment/postcheck_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".rds"))

# free up memory
rm(FitData)
rm(Samples)
gc()

# load posterior predictive checks model if necessary
PostCheck <- readRDS(paste0("output/assessment/postcheck_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".rds"))

# get the MCMC runs for the posterior predictive checks
PostCheckMCMC <- rbind(PostCheck$MCMC[[1]], PostCheck$MCMC[[2]], PostCheck$MCMC[[3]])  

# get the residuals, etc for each data type
StripResObs <- select(as_tibble(PostCheckMCMC), contains("Strip_res_obs[")) %>% as.matrix()
StripResSim <- select(as_tibble(PostCheckMCMC), contains("Strip_res_sim[")) %>% as.matrix()
StripVar <- select(as_tibble(PostCheckMCMC), contains("Strip_var[")) %>% as.matrix()
AoAResObs <- select(as_tibble(PostCheckMCMC), contains("AoA_res_obs[")) %>% as.matrix()
AoAResSim <- select(as_tibble(PostCheckMCMC), contains("AoA_res_sim[")) %>% as.matrix()
AoAVar <- select(as_tibble(PostCheckMCMC), contains("AoA_var[")) %>% as.matrix()
LineResObs <- select(as_tibble(PostCheckMCMC), contains("Line_res_obs[")) %>% as.matrix()
LineResSim <- select(as_tibble(PostCheckMCMC), contains("Line_res_sim[")) %>% as.matrix()
LineVar <- select(as_tibble(PostCheckMCMC), contains("Line_var[")) %>% as.matrix()

# free up memory
rm(PostCheck, PostCheckMCMC)
gc()

# do posterior predictive checks
TStripObs <- ((StripResObs ^ 2)) %>% apply(1, sum)
TStripSim <- ((StripResSim ^ 2)) %>% apply(1, sum)
pStrip <- sum(TStripSim >= TStripObs) / length(TStripSim)
TAoAObs <- ((AoAResObs ^ 2) / 1) %>% apply(1, sum)
TAoASim <- ((AoAResSim ^ 2) / 1) %>% apply(1, sum)
pAoA <- sum(TAoASim >= TAoAObs) / length(TAoASim)
TLineObs <- ((LineResObs ^ 2) / 1) %>% apply(1, sum)
TLineSim <- ((LineResSim ^ 2) / 1) %>% apply(1, sum)
pLine <- sum(TLineSim >= TLineObs) / length(TLineSim)

# create quantile-quantile plots
# strip transects
QQStrip <- qq.plot.ci(StripResSim, StripResObs)
# all of areas searches
QQAoA <- qq.plot.ci(AoAResSim, AoAResObs)
# line transects
QQLine <- qq.plot.ci(LineResSim, LineResObs)

minX <- min(QQStrip$Simulated)
maxX <- max(QQStrip$Simulated)
minY <- min(c(QQStrip$UpperSim,QQStrip$LowerSim,QQStrip$UpperObs,QQStrip$LowerObs))
maxY <- max(c(QQStrip$UpperSim,QQStrip$LowerSim,QQStrip$UpperObs,QQStrip$LowerObs))
PlotQQStrip <- ggplot(data = QQStrip, aes(x = Simulated, y = Observed, ymin = LowerObs, ymax = UpperObs)) + geom_ribbon(alpha = 0.2) + geom_point(shape=21, size=2, fill = "blue", colour = "blue") + labs(x = "Simulated Quantiles", y = "Observed Quantiles") + theme_minimal() + geom_abline(intercept = 0, slope = 1, size = 1, linetype = 2) + geom_line(data = QQStrip, aes(x = Simulated, y = UpperSim), size = 1, color = "red") + geom_line(data = QQStrip, aes(x = Simulated, y = LowerSim), size = 1, color = "red") + scale_y_continuous(limits = c(minY, maxY)) + scale_x_continuous(limits = c(minX, maxX)) + theme(legend.position = "none", axis.text = element_text(size = 18),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + ggtitle("Strip Transects") + theme(plot.title = element_text(size=25)) + theme(plot.margin = unit(c(1,1,1,1), "cm"))

minX <- min(QQAoA$Simulated)
maxX <- max(QQAoA$Simulated)
minY <- min(c(QQAoA$UpperSim,QQAoA$LowerSim,QQAoA$UpperObs,QQAoA$LowerObs))
maxY <- max(c(QQAoA$UpperSim,QQAoA$LowerSim,QQAoA$UpperObs,QQAoA$LowerObs))
PlotQQAoA <- ggplot(data = QQAoA, aes(x = Simulated, y = Observed, ymin = LowerObs, ymax = UpperObs)) + geom_ribbon(alpha = 0.2) + geom_point(shape=21, size=2, fill = "blue", colour = "blue") + labs(x = "Simulated Quantiles", y = "Observed Quantiles") + theme_minimal() + geom_abline(intercept = 0, slope = 1, size = 1, linetype = 2) + geom_line(data = QQAoA, aes(x = Simulated, y = UpperSim), size = 1, color = "red") + geom_line(data = QQAoA, aes(x = Simulated, y = LowerSim), size = 1, color = "red") + scale_y_continuous(limits = c(minY, maxY)) + scale_x_continuous(limits = c(minX, maxX)) + theme(legend.position = "none", axis.text = element_text(size = 18),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + ggtitle("All-of-area Searches") + theme(plot.title = element_text(size=25)) + theme(plot.margin = unit(c(1,1,1,1), "cm"))

minX <- min(QQLine$Simulated)
maxX <- max(QQLine$Simulated)
minY <- min(c(QQLine$UpperSim,QQLine$LowerSim,QQLine$UpperObs,QQLine$LowerObs))
maxY <- max(c(QQLine$UpperSim,QQLine$LowerSim,QQLine$UpperObs,QQLine$LowerObs))
PlotQQLine <- ggplot(data = QQLine, aes(x = Simulated, y = Observed, ymin = LowerObs, ymax = UpperObs)) + geom_ribbon(alpha = 0.2) + geom_point(shape=21, size=2, fill = "blue", colour = "blue") + labs(x = "Simulated Quantiles", y = "Observed Quantiles") + theme_minimal() + geom_abline(intercept = 0, slope = 1, size = 1, linetype = 2) + geom_line(data = QQLine, aes(x = Simulated, y = UpperSim), size = 1, color = "red") + geom_line(data = QQLine, aes(x = Simulated, y = LowerSim), size = 1, color = "red") + scale_y_continuous(limits = c(minY, maxY)) + scale_x_continuous(limits = c(minX, maxX)) + theme(legend.position = "none", axis.text = element_text(size = 18),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + ggtitle("Line Transects") + theme(plot.title = element_text(size=25)) + theme(plot.margin = unit(c(1,1,1,1), "cm"))

QQPlots <- PlotQQStrip + PlotQQAoA + PlotQQLine + plot_layout(nrow = 2, ncol = 2) + plot_spacer()

# NEED TO ADD A LEGEND TO THE PLOTS - TO DO

# save plots
ggsave(QQPlots, file = paste0("output/assessment/qqplots_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".jpg"), width = 40, height = 40, units = "cm", dpi = 300)
ggsave(QQPlots, file = paste0("output/assessment/qqplots_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".pdf"), width = 40, height = 40, units = "cm", dpi = 300)

# generate table of parameter estimates and variable selection frequences from best model

# categorical variable definitions
# hhgde: 1 = Intermittent and freshwater [reference category], 2 = Brackish, saline or fluctuating salinity with intermittent connectivity, 3 = Near-permanent and freshwater, 4 = Permanent and fluctuating salinity, 5 = Permanent and freshwater, 6 = Exclusion and recharge zones
# hseas: 0 = non-breeding season (winter - April to September) [reference category], 1 = breading season (summer - October to March)
# hhkha: 1 = Remnant high suitability core [reference category], 2 = Remnant medium suitability core, 3 = Remnant low suitability non-core, 4 = Remnant rainforest and non-habitat, 5 = Non-remnant high suitability core, 6 = Non-remnant medium suitability core, 7 = Non-remnant low suitability non-core, 8 = Non-remnant rainforest and non-habitat
# htlus: 1 = Natural [reference category], 2 = Production from natural, 3 = Intensive, 4 = Other (plantation, agriculture, wetland, water)

# extract MCMC chains and stitch MCMC chains together
MCMC <- rbind(BestModelSel$MCMC[[1]], BestModelSel$MCMC[[2]], BestModelSel$MCMC[[3]])

# get the names of the variables
DensPredNames <- c(BestModelSel$Data$NamesX, BestModelSel$Data$NamesY)
ObsPredNames <- c(BestModelSel$Data$NamesZ[2:length(BestModelSel$Data$NamesZ)])
AllPredNames <- c(DensPredNames, ObsPredNames)

# get the variables from the MCMC chains
DensPred_MCMC <- select(as_tibble(MCMC),contains("beta_d[")) 
ObsPred_MCMC <- select(as_tibble(MCMC),contains("beta_shn[")) %>% select(-`beta_shn[1]`)
AllPred_MCMC <- bind_cols(DensPred_MCMC, ObsPred_MCMC)
colnames(AllPred_MCMC) <- AllPredNames

# get the mean, lower, and upper quantiles, and variable selection frequencies
Mean <- AllPred_MCMC %>% summarise_all(mean) %>% t()
colnames(Mean) <- "Mean"
Lower <- AllPred_MCMC %>% summarise_all(~quantile(.x, 0.025)) %>% t()
colnames(Lower) <- "Lower"
Upper <- AllPred_MCMC %>% summarise_all(~quantile(.x, 0.975)) %>% t()
colnames(Upper) <- "Upper"
SelFreq <- AllPred_MCMC %>% summarise_all(function(x) {mean(ifelse(x == 0, 0, 1))}) %>% t()
colnames(SelFreq) <- "SelFreq"
Summary_Stats <- bind_cols(Mean, Lower, Upper, SelFreq) %>% mutate(Param = AllPredNames) %>% select(Param, Mean, Lower, Upper, SelFreq)
# write table
write_csv(Summary_Stats, paste0("output/inference/estimates_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".csv"))

# generate predictions

# set the seed
set.seed(seed)

# sub sample MCMC chains if needed
MCMC <- MCMC[sample(1:dim(MCMC)[1], 1000), ]

# load feature class of small grids
Grid <- vect("input/survey_data/grid_vec.shp")

# spatio-temporal predictions masking out rainforest

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
  Data <- get_prediction_data(Year = i, BestModelSel$Data, PredData, RainForestMask = RainMask)

  # generate predictions
  Preds <- get_predictions(MCMC, Data)

  # set masked grid cells to NA (-9999 so it can be recognised in a shapefile)
  PredsShp <- Preds$Spatial %>% mutate(Expected = ifelse((Expected == 0) & (LowerCI == 0) & (UpperCI == 0), -9999, Expected), LowerCI = ifelse((Expected == 0) & (LowerCI == 0) & (UpperCI == 0), -9999, LowerCI), UpperCI = ifelse((Expected == 0) & (LowerCI == 0) & (UpperCI == 0), -9999, UpperCI), SD = ifelse((Expected == 0) & (LowerCI == 0) & (UpperCI == 0), -9999, SD), CV = ifelse((Expected == 0) & (LowerCI == 0) & (UpperCI == 0), -9999, CV))

  # merge predictions with spatial grid
  PredsGrid <- Grid %>% merge(PredsShp, all.x = TRUE, by.x = "GridID", by.y = "GridID")

  # remove NA values
  PredsGrid <- PredsGrid %>% filter(!is.na(Expected))

  # save vector
  writeVector(PredsGrid, paste0("output/predictions/sel_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, "_year", i, ".shp"), overwrite = TRUE)
  
  # append new predictions to list
  PredsList <- append(PredsList, list(Preds))

  # free up memory
  rm(Data, Preds, PredsShp, PredsGrid)
  gc()
}

# save predictions

saveRDS(PredsList, paste0("output/predictions/sel_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".rds"))

# get the change between two dates

# 2015 to end

# specify dates

Year1 <- 2015
Year2 <- year(LastDate)

# get prediction data for the two years
Data1 <- get_prediction_data(Year = Year1, BestModelSel$Data, PredData, RainForestMask = RainMask)
Data2 <- get_prediction_data(Year = Year2, BestModelSel$Data, PredData, RainForestMask = RainMask)

# generate change estimates
Change <- get_change(MCMC, Data1, Data2)

# save change
saveRDS(Change, paste0("output/predictions/change_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_", Year1, "_", Year2, ".rds")) 

# free up memory
rm(Change)
gc()

# 2018 to end

# specify dates

Year1 <- 2018
Year2 <- year(LastDate)

# get prediction data for the two years
Data1 <- get_prediction_data(Year = Year1, BestModelSel$Data, PredData, RainForestMask = RainMask)
Data2 <- get_prediction_data(Year = Year2, BestModelSel$Data, PredData, RainForestMask = RainMask)

# generate change estimates
Change <- get_change(MCMC, Data1, Data2)

# save change
saveRDS(Change, paste0("output/predictions/change_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_", Year1, "_", Year2, ".rds")) 

# free up memory
rm(Change)
gc()

# 2020 to end

# specify dates

Year1 <- 2020
Year2 <- year(LastDate)

# get prediction data for the two years
Data1 <- get_prediction_data(Year = Year1, BestModelSel$Data, PredData, RainForestMask = RainMask)
Data2 <- get_prediction_data(Year = Year2, BestModelSel$Data, PredData, RainForestMask = RainMask)

# generate change estimates
Change <- get_change(MCMC, Data1, Data2)

# save change
saveRDS(Change, paste0("output/predictions/change_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_", Year1, "_", Year2, ".rds")) 

# free up memory
rm(Change)
gc()

# create some plots

# trend plots 2020 to end

# load outputs if necessary
PredsList <- readRDS(paste0("output/predictions/sel_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".rds"))

Years <- 2020:year(LastDate)
Expected <- lapply(PredsList[(length(PredsList) - 3):length(PredsList)], FUN = function(x){return(x$Total$Mean)}) %>% unlist()
Lower <- lapply(PredsList[(length(PredsList) - 3):length(PredsList)], FUN = function(x){return(x$Total$Lower)}) %>% unlist()
Upper <- lapply(PredsList[(length(PredsList) - 3):length(PredsList)], FUN = function(x){return(x$Total$Upper)}) %>% unlist()
Abundances <- tibble(Year = Years, Expected = Expected, Lower = Lower, Upper = Upper)

ExpectedNC <- lapply(PredsList[(length(PredsList) - 3):length(PredsList)], FUN = function(x){return(x$NC$Mean)}) %>% unlist()
LowerNC <- lapply(PredsList[(length(PredsList) - 3):length(PredsList)], FUN = function(x){return(x$NC$Lower)}) %>% unlist()
UpperNC <- lapply(PredsList[(length(PredsList) - 3):length(PredsList)], FUN = function(x){return(x$NC$Upper)}) %>% unlist()
AbundancesNC <- tibble(Year = Years, Expected = ExpectedNC, Lower = LowerNC, Upper = UpperNC)

ExpectedWI <- lapply(PredsList[(length(PredsList) - 3):length(PredsList)], FUN = function(x){return(x$WI$Mean)}) %>% unlist()
LowerWI <- lapply(PredsList[(length(PredsList) - 3):length(PredsList)], FUN = function(x){return(x$WI$Lower)}) %>% unlist()
UpperWI <- lapply(PredsList[(length(PredsList) - 3):length(PredsList)], FUN = function(x){return(x$WI$Upper)}) %>% unlist()
AbundancesWI <- tibble(Year = Years, Expected = ExpectedWI, Lower = LowerWI, Upper = UpperWI)

ExpectedSC <- lapply(PredsList[(length(PredsList) - 3):length(PredsList)], FUN = function(x){return(x$SC$Mean)}) %>% unlist()
LowerSC <- lapply(PredsList[(length(PredsList) - 3):length(PredsList)], FUN = function(x){return(x$SC$Lower)}) %>% unlist()
UpperSC <- lapply(PredsList[(length(PredsList) - 3):length(PredsList)], FUN = function(x){return(x$SC$Upper)}) %>% unlist()
AbundancesSC <- tibble(Year = Years, Expected = ExpectedSC, Lower = LowerSC, Upper = UpperSC)

PlotT <- ggplot(Abundances, aes(x = Year, y = Expected, ymin = Lower, ymax = Upper)) + geom_ribbon(alpha = 0.2, aes(fill = "95% Credible Interval")) + geom_line(aes(colour = "Expected")) + geom_point(shape=21, size=2, fill = "blue", colour = "blue") + theme_minimal() + labs(x = "Year", y = "Number of Koalas") + theme(axis.text = element_text(size = 16),  axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18, vjust = -1)) + scale_y_continuous(labels = scales::comma) + ggtitle("Total") + theme(plot.title = element_text(size=22)) + theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_colour_manual("", values = "black") + scale_fill_manual("", values = "grey12") + theme(legend.text=element_text(size=18))

PlotNC <- ggplot(AbundancesNC, aes(x = Year, y = Expected, ymin = Lower, ymax = Upper)) + geom_ribbon(alpha = 0.2, aes(fill = "95% Credible Interval")) + geom_line(aes(colour = "Expected")) + geom_point(shape=21, size=2, fill = "blue", colour = "blue") + theme_minimal() + labs(x = "Year", y = "Number of Koalas") + theme(axis.text = element_text(size = 16),  axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18, vjust = -1)) + scale_y_continuous(labels = scales::comma) + ggtitle("Northern Coast") + theme(plot.title = element_text(size=22)) + theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_colour_manual("", values = "black") + scale_fill_manual("", values = "grey12") + theme(legend.text=element_text(size=18))

PlotWI <- ggplot(AbundancesWI, aes(x = Year, y = Expected, ymin = Lower, ymax = Upper)) + geom_ribbon(alpha = 0.2, aes(fill = "95% Credible Interval")) + geom_line(aes(colour = "Expected")) + geom_point(shape=21, size=2, fill = "blue", colour = "blue") + theme_minimal() + labs(x = "Year", y = "Number of Koalas") + theme(axis.text = element_text(size = 16),  axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18, vjust = -1)) + scale_y_continuous(labels = scales::comma) + ggtitle("Western Inland") + theme(plot.title = element_text(size=22)) + theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_colour_manual("", values = "black") + scale_fill_manual("", values = "grey12") + theme(legend.text=element_text(size=18))

PlotSC <- ggplot(AbundancesSC, aes(x = Year, y = Expected, ymin = Lower, ymax = Upper)) + geom_ribbon(alpha = 0.2, aes(fill = "95% Credible Interval")) + geom_line(aes(colour = "Expected")) + geom_point(shape=21, size=2, fill = "blue", colour = "blue") + theme_minimal() + labs(x = "Year", y = "Number of Koalas") + theme(axis.text = element_text(size = 16),  axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18, vjust = -1)) + scale_y_continuous(labels = scales::comma) + ggtitle("Southern Coast") + theme(plot.title = element_text(size=22)) + theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_colour_manual("", values = "black") + scale_fill_manual("", values = "grey12") + theme(legend.text=element_text(size=18))

TrendPlot_2020_2023 <- PlotT + PlotNC + PlotWI + PlotSC + plot_layout(nrow = 2, ncol = 2, guides = "collect") & theme(legend.position = 'bottom')

ggsave(TrendPlot_2020_2023, file = paste0("output/inference/figures/trend_2020_", year(LastDate), ".jpg"), width = 40, height = 30, units = "cm", dpi = 300)

# trend plots 1996 to end

# load outputs if necessary
PredsList <- readRDS(paste0("output/predictions/sel_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_firstdate", FirstDate, ".rds"))

Years <- 1996:year(LastDate)
Expected <- lapply(PredsList[1:length(PredsList)], FUN = function(x){return(x$Total$Mean)}) %>% unlist()
Lower <- lapply(PredsList[1:length(PredsList)], FUN = function(x){return(x$Total$Lower)}) %>% unlist()
Upper <- lapply(PredsList[1:length(PredsList)], FUN = function(x){return(x$Total$Upper)}) %>% unlist()
Abundances <- tibble(Year = Years, Expected = Expected, Lower = Lower, Upper = Upper)

PlotT <- ggplot(Abundances, aes(x = Year, y = Expected, ymin = Lower, ymax = Upper)) + geom_ribbon(alpha = 0.2, aes(fill = "95% Credible Interval")) + geom_line(aes(colour = "Expected")) + geom_point(shape=21, size=2, fill = "blue", colour = "blue") + theme_minimal() + labs(x = "Year", y = "Number of Koalas") + theme(axis.text = element_text(size = 16),  axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18, vjust = -1)) + scale_y_continuous(labels = scales::comma) + ggtitle("Total") + theme(plot.title = element_text(size=22)) + theme(plot.margin = unit(c(1,1,1,1), "cm")) + scale_colour_manual("", values = "black") + scale_fill_manual("", values = "grey12") + theme(legend.text=element_text(size=18), legend.position = "bottom")

ggsave(PlotT, file = paste0("output/inference/figures/trend_1996_", year(LastDate), ".jpg"), width = 40, height = 30, units = "cm", dpi = 300)

# violin plots of change since 2020

# load outputs if needed
Change_2015_end <- readRDS(paste0("output/predictions/change_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_", 2015, "_", year(LastDate), ".rds"))
Change_2018_end <- readRDS(paste0("output/predictions/change_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_", 2018, "_", year(LastDate), ".rds"))
Change_2020_end <- readRDS(paste0("output/predictions/change_order", ModelWAICs$Order[BestIndex], "_lag", ModelWAICs$Lag[BestIndex], "_vartrend", ModelWAICs$VarTrend[BestIndex], "_", 2020, "_", year(LastDate), ".rds"))

# change 2020 - end
ChTotal <- as_tibble(Change_2020_end$DistTot) %>% rename(Change = value) %>% mutate(Region = "Total", Change = (Change - 1) * 100)
ChNC <- as_tibble(Change_2020_end$DistNC) %>% rename(Change = value) %>% mutate(Region = "Northern Coast", Change = (Change - 1) * 100)
ChWI <- as_tibble(Change_2020_end$DistWI) %>% rename(Change = value) %>% mutate(Region = "Western Inland", Change = (Change - 1) * 100)
ChSC <- as_tibble(Change_2020_end$DistSC) %>% rename(Change = value) %>% mutate(Region = "Southern Coast", Change = (Change - 1) * 100)
ChAll <- bind_rows(ChTotal, ChNC, ChWI, ChSC) %>% mutate(Region = factor(Region, levels = c("Total", "Northern Coast", "Western Inland", "Southern Coast")))

# Cat 1 and Cat 2 scenario
Plot <- ggplot(ChAll, aes(x = Region, y = Change, fill = Region)) + geom_violin(color = NA) + theme_minimal() + theme(legend.position = "none") + geom_hline(yintercept = 0) + labs(x = "Sub-region", y = "Percentage Change") + theme(axis.text = element_text(size = 16),  axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18, vjust = -1)) + scale_y_continuous(limits = c(-100, 250), breaks = seq(-100, 250, by = 25))

ggsave(Plot, file = "output/inference/figures/change_violin_2020_end.jpg", width = 40, height = 30, units = "cm", dpi = 300)
