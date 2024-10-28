# THIS SCRIPT FITS THE MODELS 

# load libraries
library(tidyverse)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)
library(MCMCvis)
library(corrplot)

# read utility functions
source("functions.R")

# read nimble code
source("nimble_code.R")

# get survey data
Surveys <- readRDS("input/survey_data_500/master.rds")

# get first date
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))

# free up memory
rm(Surveys)
gc()

# fit the saturated models

for (Order in 1:1) {
  for (Lag in 0:2) {
    for (VarTrend in 0:1) {
      
      # get the data for fitting models
      FitData <- readRDS(paste0("input/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

      # fit the model
      Samples <- get_mcmc(FitData = FitData, Type = "sat", Iter = 15000, Burnin = 5000, Thin = 1)

      # save samples
      saveRDS(Samples, paste0("output/mcmc/", "sat_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

      # free up memory
      rm(FitData)
      rm(Samples)
      gc()
    }
  }
}


















# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(coda)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

# fit model
if (VarTrend == 1) {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = sample.int(100, 3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_std", "std", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
} else {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = sample.int(100, 3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_td", "td", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
}

#stop cluster
stopCluster(cl)

# combine samples and save
Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
saveRDS(Samples_Combined, paste0("output/mcmc/sat_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))











# SAT MODELS

# ORDER = 1, LAG = 0, VARTREND = 0

# load libraries
library(tidyverse)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)

# read utility functions
source("functions.R")

# read nimble code
source("nimble_code.R")

# load input data
Surveys <- readRDS("input/survey_data_500/master.rds")
GridFrac <- readRDS("input/survey_data_500/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data_500/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data_500/cov_temporal_array_surveylocations.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
Order = 1
Lag = 0
VarTrend = 0

set.seed(20)

FitData <- get_fit_data(Surveys = Surveys, GridFrac = GridFrac, CovConsSurv = CovConsSurv, CovTempSurv = CovTempSurv, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup, Order = Order, Lag = Lag, VarTrend = VarTrend, FirstDate = FirstDate, LastDate = LastDate, StaticVars = c("htslo", "hspc1", "hspc2", "hcltp", "hcltt", "hhgde"), DynamicVars = c("hhfwc", "hhpgr2km", "htpls2km", "hcpre", "hctmn", "hseas", "hhkha", "htlus"))

saveRDS(FitData, paste0("output/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(coda)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

# fit model
if (VarTrend == 1) {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = sample.int(100, 3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_std", "std", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
} else {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = sample.int(100, 3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_td", "td", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
}

#stop cluster
stopCluster(cl)

# combine samples and save
Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
saveRDS(Samples_Combined, paste0("output/mcmc/sat_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# ORDER = 1, LAG = 0, VARTREND = 1

# load libraries
library(tidyverse)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)
library(MCMCvis)
library(corrplot)
library(abind)
library(mice)
library(factoextra)

# read utility functions
source("functions.R")

# read nimble code
source("nimble_code.R")

# load input data
Surveys <- readRDS("input/survey_data_500/master.rds")
GridFrac <- readRDS("input/survey_data_500/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data_500/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data_500/cov_temporal_array_surveylocations.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
Order = 1
Lag = 0
VarTrend = 1

set.seed(20)

FitData <- get_fit_data(Surveys = Surveys, GridFrac = GridFrac, CovConsSurv = CovConsSurv, CovTempSurv = CovTempSurv, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup, Order = Order, Lag = Lag, VarTrend = VarTrend, FirstDate = FirstDate, LastDate = LastDate, StaticVars = c("htslo", "hspc1", "hspc2", "hcltp", "hcltt", "hhgde"), DynamicVars = c("hhfwc", "hhpgr2km", "htpls2km", "hcpre", "hctmn", "hseas", "hhkha", "htlus"))

saveRDS(FitData, paste0("output/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(coda)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

# fit model
if (VarTrend == 1) {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_std", "std", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
} else {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_td", "td", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
}

#stop cluster
stopCluster(cl)

# combine samples and save
Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
saveRDS(Samples_Combined, paste0("output/mcmc/sat_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# ORDER = 1, LAG = 1, VARTREND = 0

# load libraries
library(tidyverse)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)
library(MCMCvis)
library(corrplot)
library(abind)
library(mice)
library(factoextra)

# read utility functions
source("functions.R")

# read nimble code
source("nimble_code.R")

# load input data
Surveys <- readRDS("input/survey_data_500/master.rds")
GridFrac <- readRDS("input/survey_data_500/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data_500/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data_500/cov_temporal_array_surveylocations.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
Order = 1
Lag = 1
VarTrend = 0

set.seed(20)

FitData <- get_fit_data(Surveys = Surveys, GridFrac = GridFrac, CovConsSurv = CovConsSurv, CovTempSurv = CovTempSurv, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup, Order = Order, Lag = Lag, VarTrend = VarTrend, FirstDate = FirstDate, LastDate = LastDate, StaticVars = c("htslo", "hspc1", "hspc2", "hcltp", "hcltt", "hhgde"), DynamicVars = c("hhfwc", "hhpgr2km", "htpls2km", "hcpre", "hctmn", "hseas", "hhkha", "htlus"))

saveRDS(FitData, paste0("output/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(coda)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

# fit model
if (VarTrend == 1) {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_std", "std", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
} else {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_td", "td", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
}

#stop cluster
stopCluster(cl)

# combine samples and save
Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
saveRDS(Samples_Combined, paste0("output/mcmc/sat_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# ORDER = 1, LAG = 1, VARTREND = 1

# load libraries
library(tidyverse)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)
library(MCMCvis)
library(corrplot)
library(abind)
library(mice)
library(factoextra)

# read utility functions
source("functions.R")

# read nimble code
source("nimble_code.R")

# load input data
Surveys <- readRDS("input/survey_data_500/master.rds")
GridFrac <- readRDS("input/survey_data_500/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data_500/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data_500/cov_temporal_array_surveylocations.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
Order = 1
Lag = 1
VarTrend = 1

set.seed(20)

FitData <- get_fit_data(Surveys = Surveys, GridFrac = GridFrac, CovConsSurv = CovConsSurv, CovTempSurv = CovTempSurv, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup, Order = Order, Lag = Lag, VarTrend = VarTrend, FirstDate = FirstDate, LastDate = LastDate, StaticVars = c("htslo", "hspc1", "hspc2", "hcltp", "hcltt", "hhgde"), DynamicVars = c("hhfwc", "hhpgr2km", "htpls2km", "hcpre", "hctmn", "hseas", "hhkha", "htlus"))

saveRDS(FitData, paste0("output/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(coda)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

# fit model
if (VarTrend == 1) {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_std", "std", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
} else {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_td", "td", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
}

#stop cluster
stopCluster(cl)

# combine samples and save
Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
saveRDS(Samples_Combined, paste0("output/mcmc/sat_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# ORDER = 1, LAG = 2, VARTREND = 0

# load libraries
library(tidyverse)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)
library(MCMCvis)
library(corrplot)
library(abind)
library(mice)
library(factoextra)

# read utility functions
source("functions.R")

# read nimble code
source("nimble_code.R")

# load input data
Surveys <- readRDS("input/survey_data_500/master.rds")
GridFrac <- readRDS("input/survey_data_500/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data_500/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data_500/cov_temporal_array_surveylocations.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
Order = 1
Lag = 2
VarTrend = 0

set.seed(20)

FitData <- get_fit_data(Surveys = Surveys, GridFrac = GridFrac, CovConsSurv = CovConsSurv, CovTempSurv = CovTempSurv, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup, Order = Order, Lag = Lag, VarTrend = VarTrend, FirstDate = FirstDate, LastDate = LastDate, StaticVars = c("htslo", "hspc1", "hspc2", "hcltp", "hcltt", "hhgde"), DynamicVars = c("hhfwc", "hhpgr2km", "htpls2km", "hcpre", "hctmn", "hseas", "hhkha", "htlus"))

saveRDS(FitData, paste0("output/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(coda)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

# fit model
if (VarTrend == 1) {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_std", "std", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
} else {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_td", "td", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
}

#stop cluster
stopCluster(cl)

# combine samples and save
Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
saveRDS(Samples_Combined, paste0("output/mcmc/sat_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# ORDER = 1, LAG = 2, VARTREND = 1

# load libraries
library(tidyverse)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)
library(MCMCvis)
library(corrplot)
library(abind)
library(mice)
library(factoextra)

# read utility functions
source("functions.R")

# read nimble code
source("nimble_code.R")

# load input data
Surveys <- readRDS("input/survey_data_500/master.rds")
GridFrac <- readRDS("input/survey_data_500/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data_500/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data_500/cov_temporal_array_surveylocations.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
Order = 1
Lag = 2
VarTrend = 1

set.seed(20)

FitData <- get_fit_data(Surveys = Surveys, GridFrac = GridFrac, CovConsSurv = CovConsSurv, CovTempSurv = CovTempSurv, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup, Order = Order, Lag = Lag, VarTrend = VarTrend, FirstDate = FirstDate, LastDate = LastDate, StaticVars = c("htslo", "hspc1", "hspc2", "hcltp", "hcltt", "hhgde"), DynamicVars = c("hhfwc", "hhpgr2km", "htpls2km", "hcpre", "hctmn", "hseas", "hhkha", "htlus"))

saveRDS(FitData, paste0("output/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(coda)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

# fit model
if (VarTrend == 1) {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_std", "std", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
} else {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_td", "td", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
}

#stop cluster
stopCluster(cl)

# combine samples and save
Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
saveRDS(Samples_Combined, paste0("output/mcmc/sat_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# SEL MODELS

# ORDER = 1, LAG = 0, VARTREND = 0

# load libraries
library(tidyverse)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)
library(MCMCvis)
library(corrplot)
library(abind)
library(mice)
library(factoextra)

# read utility functions
source("functions.R")

# read nimble code
source("nimble_code.R")

# load input data
Surveys <- readRDS("input/survey_data_500/master.rds")
GridFrac <- readRDS("input/survey_data_500/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data_500/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data_500/cov_temporal_array_surveylocations.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
Order = 1
Lag = 0
VarTrend = 0

set.seed(20)

FitData <- get_fit_data(Surveys = Surveys, GridFrac = GridFrac, CovConsSurv = CovConsSurv, CovTempSurv = CovTempSurv, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup, Order = Order, Lag = Lag, VarTrend = VarTrend, FirstDate = FirstDate, LastDate = LastDate, StaticVars = c("htslo", "hspc1", "hspc2", "hcltp", "hcltt", "hhgde"), DynamicVars = c("hhfwc", "hhpgr2km", "htpls2km", "hcpre", "hctmn", "hseas", "hhkha", "htlus"))

saveRDS(FitData, paste0("output/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(coda)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

# fit model
if (VarTrend == 1) {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sel_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_std", "std", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
} else {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sel_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_td", "td", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
}

#stop cluster
stopCluster(cl)

# combine samples and save
Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
saveRDS(Samples_Combined, paste0("output/mcmc/sel_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# ORDER = 1, LAG = 0, VARTREND = 1

# load libraries
library(tidyverse)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)
library(MCMCvis)
library(corrplot)
library(abind)
library(mice)
library(factoextra)

# read utility functions
source("functions.R")

# read nimble code
source("nimble_code.R")

# load input data
Surveys <- readRDS("input/survey_data_500/master.rds")
GridFrac <- readRDS("input/survey_data_500/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data_500/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data_500/cov_temporal_array_surveylocations.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
Order = 1
Lag = 0
VarTrend = 1

set.seed(20)

FitData <- get_fit_data(Surveys = Surveys, GridFrac = GridFrac, CovConsSurv = CovConsSurv, CovTempSurv = CovTempSurv, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup, Order = Order, Lag = Lag, VarTrend = VarTrend, FirstDate = FirstDate, LastDate = LastDate, StaticVars = c("htslo", "hspc1", "hspc2", "hcltp", "hcltt", "hhgde"), DynamicVars = c("hhfwc", "hhpgr2km", "htpls2km", "hcpre", "hctmn", "hseas", "hhkha", "htlus"))

saveRDS(FitData, paste0("output/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(coda)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

# fit model
if (VarTrend == 1) {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sel_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_std", "std", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
} else {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sel_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_td", "td", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
}

#stop cluster
stopCluster(cl)

# combine samples and save
Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
saveRDS(Samples_Combined, paste0("output/mcmc/sel_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# ORDER = 1, LAG = 1, VARTREND = 0

# load libraries
library(tidyverse)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)
library(MCMCvis)
library(corrplot)
library(abind)
library(mice)
library(factoextra)

# read utility functions
source("functions.R")

# read nimble code
source("nimble_code.R")

# load input data
Surveys <- readRDS("input/survey_data_500/master.rds")
GridFrac <- readRDS("input/survey_data_500/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data_500/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data_500/cov_temporal_array_surveylocations.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
Order = 1
Lag = 1
VarTrend = 0

set.seed(20)

FitData <- get_fit_data(Surveys = Surveys, GridFrac = GridFrac, CovConsSurv = CovConsSurv, CovTempSurv = CovTempSurv, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup, Order = Order, Lag = Lag, VarTrend = VarTrend, FirstDate = FirstDate, LastDate = LastDate, StaticVars = c("htslo", "hspc1", "hspc2", "hcltp", "hcltt", "hhgde"), DynamicVars = c("hhfwc", "hhpgr2km", "htpls2km", "hcpre", "hctmn", "hseas", "hhkha", "htlus"))

saveRDS(FitData, paste0("output/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(coda)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

# fit model
if (VarTrend == 1) {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sel_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_std", "std", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
} else {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sel_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_td", "td", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
}

#stop cluster
stopCluster(cl)

# combine samples and save
Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
saveRDS(Samples_Combined, paste0("output/mcmc/sel_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# ORDER = 1, LAG = 1, VARTREND = 1

# load libraries
library(tidyverse)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)
library(MCMCvis)
library(corrplot)
library(abind)
library(mice)
library(factoextra)

# read utility functions
source("functions.R")

# read nimble code
source("nimble_code.R")

# load input data
Surveys <- readRDS("input/survey_data_500/master.rds")
GridFrac <- readRDS("input/survey_data_500/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data_500/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data_500/cov_temporal_array_surveylocations.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
Order = 1
Lag = 1
VarTrend = 1

set.seed(20)

FitData <- get_fit_data(Surveys = Surveys, GridFrac = GridFrac, CovConsSurv = CovConsSurv, CovTempSurv = CovTempSurv, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup, Order = Order, Lag = Lag, VarTrend = VarTrend, FirstDate = FirstDate, LastDate = LastDate, StaticVars = c("htslo", "hspc1", "hspc2", "hcltp", "hcltt", "hhgde"), DynamicVars = c("hhfwc", "hhpgr2km", "htpls2km", "hcpre", "hctmn", "hseas", "hhkha", "htlus"))

saveRDS(FitData, paste0("output/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(coda)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

# fit model
if (VarTrend == 1) {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sel_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_std", "std", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
} else {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sel_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_td", "td", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
}

#stop cluster
stopCluster(cl)

# combine samples and save
Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
saveRDS(Samples_Combined, paste0("output/mcmc/sel_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# ORDER = 1, LAG = 2, VARTREND = 0

# load libraries
library(tidyverse)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)
library(MCMCvis)
library(corrplot)
library(abind)
library(mice)
library(factoextra)

# read utility functions
source("functions.R")

# read nimble code
source("nimble_code.R")

# load input data
Surveys <- readRDS("input/survey_data_500/master.rds")
GridFrac <- readRDS("input/survey_data_500/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data_500/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data_500/cov_temporal_array_surveylocations.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
Order = 1
Lag = 2
VarTrend = 0

set.seed(20)

FitData <- get_fit_data(Surveys = Surveys, GridFrac = GridFrac, CovConsSurv = CovConsSurv, CovTempSurv = CovTempSurv, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup, Order = Order, Lag = Lag, VarTrend = VarTrend, FirstDate = FirstDate, LastDate = LastDate, StaticVars = c("htslo", "hspc1", "hspc2", "hcltp", "hcltt", "hhgde"), DynamicVars = c("hhfwc", "hhpgr2km", "htpls2km", "hcpre", "hctmn", "hseas", "hhkha", "htlus"))

saveRDS(FitData, paste0("output/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(coda)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

# fit model
if (VarTrend == 1) {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sel_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_std", "std", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
} else {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sel_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_td", "td", "nbr"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
}

#stop cluster
stopCluster(cl)

# combine samples and save
Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
saveRDS(Samples_Combined, paste0("output/mcmc/sel_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# ORDER = 1, LAG = 2, VARTREND = 1

# load libraries
library(tidyverse)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)
library(MCMCvis)
library(corrplot)
library(abind)
library(mice)
library(factoextra)

# read utility functions
source("functions.R")

# read nimble code
source("nimble_code.R")

# load input data
Surveys <- readRDS("input/survey_data_500/master.rds")
GridFrac <- readRDS("input/survey_data_500/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data_500/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data_500/cov_temporal_array_surveylocations.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
Order = 1
Lag = 2
VarTrend = 1

set.seed(20)

FitData <- get_fit_data(Surveys = Surveys, GridFrac = GridFrac, CovConsSurv = CovConsSurv, CovTempSurv = CovTempSurv, DateIntervals = DateIntervals, GenPopLookup = GenPopLookup, Order = Order, Lag = Lag, VarTrend = VarTrend, FirstDate = FirstDate, LastDate = LastDate, StaticVars = c("htslo", "hspc1", "hspc2", "hcltp", "hcltt", "hhgde"), DynamicVars = c("hhfwc", "hhpgr2km", "htpls2km", "hcpre", "hctmn", "hseas", "hhkha", "htlus"))

saveRDS(FitData, paste0("input/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(coda)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

# fit model
if (VarTrend == 1) {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sel_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_std", "std", "nbrst", "nbraa", "nbrli"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
} else {
	Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sel_model, Seeds = sample.int(100,3), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_td", "td", "nbrst", "nbraa", "nbrli"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
}

#stop cluster
stopCluster(cl)

# combine samples and save
Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
saveRDS(Samples_Combined, paste0("output/mcmc/sel_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

