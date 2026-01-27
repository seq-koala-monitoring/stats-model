# THIS SCRIPT FITS THE MODELS

# clean global environment
rm(list=ls())
try(dev.off(dev.list()["RStudioGD"]), silent=TRUE)
gc()

# redirect messages to a log file
log_file <- "output/log_model_runs.txt"
cat("\n==== Log started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "====\n\n",
    file = log_file, append = FALSE)
log_con <- file(log_file, open = "a")
sink(log_con)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

# install and load packages
packages <- c("tidyverse", "nimble", "coda", "extraDistr", "parallel", "MCMCvis", "corrplot")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, quiet = T)
invisible(lapply(packages, library, character.only = TRUE))

# load parameters
source("parameters_init.txt")

# read utility functions
source("code/functions.R")

# read nimble code
source("code/nimble_code.R")

# get survey data
Surveys <- readRDS("input/survey_data/master.rds")

# get first date
if(is.null(FirstDate)){
  FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
} else {
  FirstDate <- as.Date(FirstDate, format = "%d/%m/%Y")
}
# free up memory
rm(Surveys)
gc()

# fit the saturated models

for (Order in 1:1) {
  for (Lag in 0:2) {
    for (VarTrend in 0:1) {
      # get the data for fitting models
      FitData <- readRDS(paste0("input/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))
      
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
        Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = c(seed, seed + 20, seed + 40), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_sd", "sd", "sigma_std", "std", "nbrst", "nbraa", "nbrli"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
      } else {
        Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sat_model, Seeds = c(seed, seed + 20, seed + 40), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_sd", "sd", "sigma_td", "td", "nbrst", "nbraa", "nbrli"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
      }

      #stop cluster
      stopCluster(cl)

      # combine samples and save
      Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
      saveRDS(Samples_Combined, paste0("output/mcmc/sat_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

      # free up memory
      rm(FitData)
      rm(Samples)
      rm(Samples_Combined)
      gc()
    }
  }
}

# fit the variable selection models

for (Order in 1:1) {
  for (Lag in 0:2) {
    for (VarTrend in 0:1) {
      # get the data for fitting models
      FitData <- readRDS(paste0("input/nimble_data/data_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))
      
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
        Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sel_model, Seeds = c(seed, seed + 20, seed + 40), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_sd", "sd", "sigma_std", "std", "nbrst", "nbraa", "nbrli"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
      } else {
        Samples <- parLapply(cl = cl, X = 1:3, fun = fit_sel_model, Seeds = c(seed, seed + 20, seed + 40), Iter = 15000, Burnin = 5000, Thin = 1, Monitors = c("beta_d", "beta_shn", "sigma_sd", "sd", "sigma_td", "td", "nbrst", "nbraa", "nbrli"), Calculate = FALSE, EnableWAIC = TRUE, Data = FitData, Code = nimble_sat_model)
      }

      #stop cluster
      stopCluster(cl)

      # combine samples and save
      Samples_Combined <- list(MCMC = list(Samples[[1]]$Samples$samples, Samples[[2]]$Samples$samples, Samples[[3]]$Samples$samples), WAIC = list(Samples[[1]]$Samples$WAIC, Samples[[2]]$Samples$WAIC, Samples[[3]]$Samples$WAIC), Data = Samples[[1]]$Data, Code = Samples[[1]]$Code)
      saveRDS(Samples_Combined, paste0("output/mcmc/sel_order", Order, "_lag", Lag, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))

      # free up memory
      rm(FitData)
      rm(Samples)
      rm(Samples_Combined)
      gc()
    }
  }
}

# end
print("model runs complete")

# reset sink and close connection
sink(type = "message")
sink()
close(log_con)
