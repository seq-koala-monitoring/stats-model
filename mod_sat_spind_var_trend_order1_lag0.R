# load libraries
library(tidyverse)
library(nimble)
library(extraDistr)
library(parallel)

# read utility functions
source("functions.R")

# read nimble code
source("nimble_code.R")

# Model: saturated model with independent spatial random-effect spatial variation in trends and with order = 1 and lag = 0 with 3 chains

# get data inputs
NimbleInputs <- readRDS("output/nimble_data/data_order1_lag0.rds")

# set number of iterations, burnin, number chains, and thin interval
N.iter <- 10000
N.burnin <- 4000
N.thin <- 1

# set seeds for creating on initial values and for MCMC sampling for the three chains
Seeds <- c(666, 777, 888)

# make cluster
cl <- makeCluster(3)

# export packages and nimble functions to cluster
clusterEvalQ(cl,{
  library(tidyverse)
  library(nimble)
  library(extraDistr)})
clusterExport(cl, {"cdfhnorm"})

# fit model in parallel (3 cores, 3 chains)
Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_spind_var_trend, Seeds = Seeds, Data = NimbleInputs, Code = code_sat_spind_var_trend, Iter = N.iter, Burnin = N.burnin, Thin = N.thin, Monitors = c("beta_d", "beta_shn", "mean_mud"))

stopCluster(cl)

saveRDS(Samples, "output/mcmc/sat_spind_var_trend_order1_lag0.rds")
