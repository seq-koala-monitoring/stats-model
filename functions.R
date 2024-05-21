# functions

# nimble function to get the cdf of f(x) for distance sampling
cdfhnorm <- nimbleRcall(function(q = double(0), sigma = double(0)){}, Rfun = 'phnorm', returnType = double(0))

# 1. fit model with independent spatial random effects and no spatial variation in trends (no model selection
fit_sat_spind_const_trend <- function(X, Seeds, Data, Code, Iter, Burnin, Thin, Monitors) {
# note that Monitors must at least contain "beta_d" and "beta_shn" as a minimum

  # set seed for creating initial values
  set.seed(Seeds[X])

  # get initial values for parameters
  Inits <- list(mu_sd = runif(n = 1, -2, 2), sigma_sd = runif(n = 1,0, 5), sd = rep(0, Data$Constants$NGPops), sigma_td = runif(n = 1,0, 5), td = rep(0, Data$Constants$NTime), beta_d = runif(n = (Data$Constants$NX + Data$Constants$NY), -2, 2), sigma_d = runif(n = 1,0, 5), d = matrix(10, nrow = Data$Constants$NSGrids, ncol = (Data$Constants$LastDateID - Data$Constants$FirstDateID + 1)), PDists = ifelse(is.na(Data$Data$PDists), 15, NA), int_shn = runif(n = 1, -2, 2), beta_shn = runif(n = Data$Constants$NZ, -2, 2))

  # set up nimble mcmc, samplers and compile
  NimbleModel <- nimbleModel(code = Code, constants = Data$Constants, data = Data$Data, inits = Inits, calculate = TRUE)
  CNimbleModel <- compileNimble(NimbleModel)
  NimbleModelConf <- configureMCMC(NimbleModel, monitors = Monitors, enableWAIC = TRUE)
  # change default samplers to a slice sampler to improve mixing
  NimbleModelConf$removeSamplers(c("beta_d", "beta_shn"))
  NimbleModelConf$addSampler(c("beta_d", "beta_shn"), type = "AF_slice")
  NimbleModelMCMC <- buildMCMC(NimbleModelConf)
  CNimbleModelMCMC <- compileNimble(NimbleModelMCMC, project = NimbleModel, resetFunctions = TRUE)
  Samples <- runMCMC(mcmc = CNimbleModelMCMC, niter = Iter, nburnin = Burnin, nchains = 1, thin = Thin, setSeed = Seeds[X], WAIC = TRUE)

  # return samples
  return(Samples)
}

# 2. fit model with independent spatial random effects and spatial variation in trends
fit_sat_spind_var_trend <- function(X, Seeds, Data, Code, Iter, Burnin, Thin, Monitors) {
# note that Monitors must at least contain "beta_d" and "beta_shn" as a minimum

  # set seed for creating initial values
  set.seed(Seeds[X])

  # get initial values for parameters
  Inits <- list(mu_sd = runif(n = 1, -2, 2), sigma_sd = runif(n = 1,0, 5), sd = rep(0, Data$Constants$NGPops), sigma_std = runif(n = 1,0, 5), std = matrix(0, nrow = Data$Constants$NGPops, ncol = Data$Constants$NTime), beta_d = runif(n = (Data$Constants$NX + Data$Constants$NY), -2, 2), sigma_d = runif(n = 1,0, 5), d = matrix(10, nrow = Data$Constants$NSGrids, ncol = (Data$Constants$LastDateID - Data$Constants$FirstDateID + 1)), PDists = ifelse(is.na(Data$Data$PDists), 15, NA), int_shn = runif(n = 1, -2, 2), beta_shn = runif(n = Data$Constants$NZ, -2, 2))

  # set up nimble mcmc, samplers and compile
  NimbleModel <- nimbleModel(code = Code, constants = Data$Constants, data = Data$Data, inits = Inits, calculate = TRUE)
  CNimbleModel <- compileNimble(NimbleModel)
  NimbleModelConf <- configureMCMC(NimbleModel, monitors = Monitors, enableWAIC = TRUE)
  # change default samplers to a slice sampler to improve mixing
  NimbleModelConf$removeSamplers(c("beta_d", "beta_shn"))
  NimbleModelConf$addSampler(c("beta_d", "beta_shn"), type = "AF_slice")
  NimbleModelMCMC <- buildMCMC(NimbleModelConf)
  CNimbleModelMCMC <- compileNimble(NimbleModelMCMC, project = NimbleModel, resetFunctions = TRUE)
  Samples <- runMCMC(mcmc = CNimbleModelMCMC, niter = Iter, nburnin = Burnin, nchains = 1, thin = Thin, setSeed = Seeds[X], WAIC = TRUE)

  # return samples
  return(Samples)
}

# 3. fit model with car process spatial random effects and no spatial variation in trends
fit_sat_spcar_const_trend <- function(X, Seeds, Data, Code, Iter, Burnin, Thin, Monitors) {
# note that Monitors must at least contain "beta_d" and "beta_shn" as a minimum

  # set seed for creating initial values
  set.seed(Seeds[X])

  # get initial values for parameters
  Inits <- list(sigma_sd = runif(n = 1,0, 5), sd = rep(0, Data$Constants$NLGrids), sigma_td = runif(n = 1,0, 5), td = rep(0, Data$Constants$NTime), beta_d = runif(n = (Data$Constants$NX + Data$Constants$NY), -2, 2), sigma_d = runif(n = 1,0, 5), d = matrix(10, nrow = Data$Constants$NSGrids, ncol = (Data$Constants$LastDateID - Data$Constants$FirstDateID + 1)), PDists = ifelse(is.na(Data$Data$PDists), 15, NA), int_shn = runif(n = 1, -2, 2), beta_shn = runif(n = Data$Constants$NZ, -2, 2))

  # set up nimble mcmc, samplers and compile
  NimbleModel <- nimbleModel(code = Code, constants = Data$Constants, data = Data$Data, inits = Inits, calculate = TRUE)
  CNimbleModel <- compileNimble(NimbleModel)
  NimbleModelConf <- configureMCMC(NimbleModel, monitors = Monitors, enableWAIC = TRUE)
  # change default samplers to a slice sampler to improve mixing
  NimbleModelConf$removeSamplers(c("beta_d", "beta_shn"))
  NimbleModelConf$addSampler(c("beta_d", "beta_shn"), type = "AF_slice")
  NimbleModelMCMC <- buildMCMC(NimbleModelConf)
  CNimbleModelMCMC <- compileNimble(NimbleModelMCMC, project = NimbleModel, resetFunctions = TRUE)
  Samples <- runMCMC(mcmc = CNimbleModelMCMC, niter = Iter, nburnin = Burnin, nchains = 1, thin = Thin, setSeed = Seeds[X], WAIC = TRUE)

  # return samples
  return(Samples)
}

# 4. fit model with car process spatial random effects and spatial variation in trends
fit_sat_spcar_var_trend <- function(X, Seeds, Data, Code, Iter, Burnin, Thin, Monitors) {
# note that Monitors must at least contain "beta_d" and "beta_shn" as a minimum

  # set seed for creating initial values
  set.seed(Seeds[X])

  # get initial values for parameters
  Inits <- list(sigma_sd = runif(n = 1,0, 5), sd = rep(0, Data$Constants$NLGrids), sigma_std = runif(n = 1,0, 5), std = matrix(0, nrow = Data$Constants$NGPops, ncol = Data$Constants$NTime), beta_d = runif(n = (Data$Constants$NX + Data$Constants$NY), -2, 2), sigma_d = runif(n = 1,0, 5), d = matrix(10, nrow = Data$Constants$NSGrids, ncol = (Data$Constants$LastDateID - Data$Constants$FirstDateID + 1)), PDists = ifelse(is.na(Data$Data$PDists), 15, NA), int_shn = runif(n = 1, -2, 2), beta_shn = runif(n = Data$Constants$NZ, -2, 2))

  # set up nimble mcmc, samplers and compile
  NimbleModel <- nimbleModel(code = Code, constants = Data$Constants, data = Data$Data, inits = Inits, calculate = TRUE)
  CNimbleModel <- compileNimble(NimbleModel)
  NimbleModelConf <- configureMCMC(NimbleModel, monitors = Monitors, enableWAIC = TRUE)
  # change default samplers to a slice sampler to improve mixing
  NimbleModelConf$removeSamplers(c("beta_d", "beta_shn"))
  NimbleModelConf$addSampler(c("beta_d", "beta_shn"), type = "AF_slice")
  NimbleModelMCMC <- buildMCMC(NimbleModelConf)
  CNimbleModelMCMC <- compileNimble(NimbleModelMCMC, project = NimbleModel, resetFunctions = TRUE)
  Samples <- runMCMC(mcmc = CNimbleModelMCMC, niter = Iter, nburnin = Burnin, nchains = 1, thin = Thin, setSeed = Seeds[X], WAIC = TRUE)

  # return samples
  return(Samples)
}

# 5. fit model with independent spatial random effects and no spatial variation in trends (with model selection)
fit_sel_spind_const_trend <- function(X, Seeds, Data, Code, Iter, Burnin, Thin, Monitors) {
# note that Monitors must at least contain "beta_d", "beta_shn", and "z" as a minimum

  # set seed for creating initial values
  set.seed(Seeds[X])

  # get initial values for parameters
  Inits <- list(mu_sd = runif(n = 1, -2, 2), sigma_sd = runif(n = 1,0, 5), sd = rep(0, Data$Constants$NGPops), sigma_td = runif(n = 1,0, 5), td = rep(0, Data$Constants$NTime), beta_d = runif(n = (Data$Constants$NX + Data$Constants$NY), -2, 2), sigma_d = runif(n = 1,0, 5), d = matrix(10, nrow = Data$Constants$NSGrids, ncol = (Data$Constants$LastDateID - Data$Constants$FirstDateID + 1)), PDists = ifelse(is.na(Data$Data$PDists), 15, NA), int_shn = runif(n = 1, -2, 2), beta_shn = runif(n = Data$Constants$NZ, -2, 2), psi = runif(n = 1, 0, 1), z = rbinom(n = (Data$Constants$NX + Data$Constants$NY), size = 1, prob = 0.5))

  # set up nimble mcmc, samplers and compile
  NimbleModel <- nimbleModel(code = Code, constants = Data$Constants, data = Data$Data, inits = Inits, calculate = TRUE)
  CNimbleModel <- compileNimble(NimbleModel)
  NimbleModelConf <- configureMCMC(NimbleModel, monitors = Monitors)
  NimbleModelMCMC <- buildMCMC(NimbleModelConf)
  CNimbleModelMCMC <- compileNimble(NimbleModelMCMC, project = NimbleModel, resetFunctions = TRUE)
  Samples <- runMCMC(mcmc = CNimbleModelMCMC, niter = Iter, nburnin = Burnin, nchains = 1, thin = Thin, setSeed = Seeds[X])

  # return samples
  return(Samples)
}
