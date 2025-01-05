# functions

# nimble function to get the cdf of f(x) for distance sampling
cdfhnorm <- nimbleRcall(function(q = double(0), sigma = double(0)){}, Rfun = 'phnorm', returnType = double(0))

# fit saturated model

# function to fit a saturated model for a single chain - set up to be used with parLapply
fit_sat_model <- function(X, Seeds, Iter, Burnin, Thin, Monitors, Calculate = FALSE, EnableWAIC = FALSE, Data, Code) {

  # get Order
  Order <- Data$Order

  # get Lag
  Lag <- Data$Lag

  # get VarTrend
  VarTrend <- Data$VarTrend

  # set up initial values
  if (Order == 2) {
    if (VarTrend == 1) {
      Inits <- list(sigma_std = runif(n = 1,0, 5), std = matrix(0, nrow = Data$Constants$NGPops, ncol = Data$Constants$NTime2), beta_d = runif(n = (Data$Constants$NX + Data$Constants$NY), -2, 2), PDists = ifelse(is.na(Data$Data$PDists), 15, NA), beta_shn = runif(n = Data$Constants$NZ, -2, 2), nbr = runif(n = 1,0, 5))
    } else {
      Inits <- list(sigma_td = runif(n = 1,0, 5), td = rep(0, Data$Constants$NTime2), beta_d = runif(n = (Data$Constants$NX + Data$Constants$NY), -2, 2), PDists = ifelse(is.na(Data$Data$PDists), 15, NA), beta_shn = runif(n = Data$Constants$NZ, -2, 2), nbr = runif(n = 1,0, 5))
    }
  } else {
    if (VarTrend == 1) {
      Inits <- list(sigma_std = runif(n = 1,0, 5), std = matrix(0, nrow = Data$Constants$NGPops, ncol = Data$Constants$NTime), beta_d = runif(n = (Data$Constants$NX + Data$Constants$NY), -2, 2), PDists = ifelse(is.na(Data$Data$PDists), 15, NA), beta_shn = runif(n = Data$Constants$NZ, -2, 2), nbr = runif(n = 1,0, 5))
    } else {
      Inits <- list(sigma_td = runif(n = 1,0, 5), td = rep(0, Data$Constants$NTime), beta_d = runif(n = (Data$Constants$NX + Data$Constants$NY), -2, 2), PDists = ifelse(is.na(Data$Data$PDists), 15, NA), beta_shn = runif(n = Data$Constants$NZ, -2, 2), nbr = runif(n = 1,0, 5))
    }
  }

  # set up nimble mcmc, samplers and compile
  NimbleModel <- nimbleModel(code = Code, constants = Data$Constants, data = Data$Data, inits = Inits, calculate = Calculate)
  CNimbleModel <- compileNimble(NimbleModel)
  NimbleModelConf <- configureMCMC(NimbleModel, monitors = Monitors, enableWAIC = EnableWAIC)
  NimbleModelMCMC <- buildMCMC(NimbleModelConf)
  CNimbleModelMCMC <- compileNimble(NimbleModelMCMC, project = NimbleModel, resetFunctions = TRUE)
  Samples <- runMCMC(mcmc = CNimbleModelMCMC, niter = Iter, nburnin = Burnin, nchains = 1, thin = Thin, setSeed = Seeds[X], WAIC = EnableWAIC)

  # return samples
  return(list(Samples = Samples, Data = Data, Code = Code))
}

# function to fit a model wtih variable selection for a single chain - set up to be used with parLapply
fit_sel_model <- function(X, Seeds, Iter, Burnin, Thin, Monitors, Calculate = FALSE, EnableWAIC = FALSE, Data, Code) {

  # get Order
  Order <- Data$Order

  # get Lag
  Lag <- Data$Lag

  # get VarTrend
  VarTrend <- Data$VarTrend

  # set up initial values
  if (Order == 2) {
    if (VarTrend == 1) {
      Inits <- list(sigma_std = runif(n = 1,0, 5), std = matrix(0, nrow = Data$Constants$NGPops, ncol = Data$Constants$NTime2), beta_d = runif(n = (Data$Constants$NX + Data$Constants$NY), -2, 2), PDists = ifelse(is.na(Data$Data$PDists), 15, NA), beta_shn = runif(n = Data$Constants$NZ, -2, 2), nbr = runif(n = 1,0, 5))
    } else {
      Inits <- list(sigma_td = runif(n = 1,0, 5), td = rep(0, Data$Constants$NTime2), beta_d = runif(n = (Data$Constants$NX + Data$Constants$NY), -2, 2), PDists = ifelse(is.na(Data$Data$PDists), 15, NA), beta_shn = runif(n = Data$Constants$NZ, -2, 2), nbr = runif(n = 1,0, 5))
    }
  } else {
    if (VarTrend == 1) {
      Inits <- list(sigma_std = runif(n = 1,0, 5), std = matrix(0, nrow = Data$Constants$NGPops, ncol = Data$Constants$NTime), beta_d = runif(n = (Data$Constants$NX + Data$Constants$NY), -2, 2), PDists = ifelse(is.na(Data$Data$PDists), 15, NA), beta_shn = runif(n = Data$Constants$NZ, -2, 2), nbr = runif(n = 1,0, 5))
    } else {
      Inits <- list(sigma_td = runif(n = 1,0, 5), td = rep(0, Data$Constants$NTime), beta_d = runif(n = (Data$Constants$NX + Data$Constants$NY), -2, 2), PDists = ifelse(is.na(Data$Data$PDists), 15, NA), beta_shn = runif(n = Data$Constants$NZ, -2, 2), nbr = runif(n = 1,0, 5))
    }
  }

  # set up nimble mcmc, samplers and compile
  NimbleModel <- nimbleModel(code = Code, constants = Data$Constants, data = Data$Data, inits = Inits, calculate = Calculate)
  CNimbleModel <- compileNimble(NimbleModel)
  NimbleModelConf <- configureMCMC(NimbleModel, monitors = Monitors, enableWAIC = EnableWAIC)
  #NimbleModelConf$addMonitors("z_d")
  #NimbleModelConf$addMonitors("z_shn")
  #configureRJ(NimbleModelConf, targetNodes = c("beta_d"), indicatorNodes = c("z_d"))
  #configureRJ(NimbleModelConf, targetNodes = c("beta_shn"), indicatorNodes = c("z_shn"))
  configureRJ(NimbleModelConf, targetNodes = c(paste0("beta_d[1:", (Data$Constants$NX + Data$Constants$NY), "]" )), priorProb = 0.5)
  configureRJ(NimbleModelConf, targetNodes = c(paste0("beta_shn[2:", Data$Constants$NZ, "]")), priorProb = 0.5)
  NimbleModelMCMC <- buildMCMC(NimbleModelConf)
  CNimbleModelMCMC <- compileNimble(NimbleModelMCMC, project = NimbleModel, resetFunctions = TRUE)
  Samples <- runMCMC(mcmc = CNimbleModelMCMC, niter = Iter, nburnin = Burnin, nchains = 1, thin = Thin, setSeed = Seeds[X], WAIC = EnableWAIC)

  # return samples
  return(list(Samples = Samples, Data = Data, Code = Code))
}

get_fit_data <- function(Surveys, GridFrac, CovConsSurv, CovTempSurv, DateIntervals, GenPopLookup, Order, Lag, VarTrend, FirstDate, LastDate, StaticVars, DynamicVars) {

  # reclassify genetic population IDs to three "populations"
  # 1 (Noosa) -> 1
  # 2 (North Coast) -> 1
  # 3 (North West) -> 2
  # 4 (Toowoomba)-> 2
  # 5 (South West) -> 2
  # 6 (Koala Coast) -> 3
  # 7 (Minjerribah) -> 3
  # 8 (Gold Coast) -> 3
  GenPopLookup <- GenPopLookup %>% as_tibble() %>% mutate(GENPOP_ID = case_match(GENPOP_ID, 1~1, 2~1, 3~2, 4~2, 5~2, 6~3, 7~3, 8~3) %>% as.integer())

  # get first and last date index values
  FirstDateID <- filter(DateIntervals, start_date <= FirstDate & end_date >= FirstDate) %>% pull(TimePeriodID)
  LastDateID <- filter(DateIntervals, start_date <= LastDate & end_date >= LastDate) %>% pull(TimePeriodID)

  # add a season variable to the temporally variable covariates
  # 1 = breading season (summer - October to March) and 0 = non-breeding season (winter - April to September)
  Seas <- ifelse(month(DateIntervals$start_date) == 10, 1, 0)
  Seas <- t(matrix(rep(Seas, dim(CovTempSurv)[1]), ncol = dim(CovTempSurv)[1], nrow = dim(CovTempSurv)[3]))
  CovTempSurv <- abind(CovTempSurv, Seas, along = 2)
  dimnames(CovTempSurv)[[2]][(dim(CovTempSurv)[[2]])] <- "hseas"

  # add a time variable to the temporally variable covariates
  Time <- 1:dim(CovTempSurv)[3]
  Time <- t(matrix(rep(Time, dim(CovTempSurv)[1]), ncol = dim(CovTempSurv)[1], nrow = dim(CovTempSurv)[3]))
  CovTempSurv <- abind(CovTempSurv, Time, along = 2)
  dimnames(CovTempSurv)[[2]][(dim(CovTempSurv)[[2]])] <- "htime"

  # Subset data for the specified date range
  # survey data - selecting only the transects for the time period selected
  Surveys$line_transect <- Surveys$line_transect %>% filter(Date >= FirstDate & Date <= LastDate)
  Surveys$strip_transect <- Surveys$strip_transect %>% filter(Date >= FirstDate & Date <= LastDate)
  Surveys$uaoa <- Surveys$uaoa %>% filter(Date >= FirstDate & Date <= LastDate)
  Surveys$perp_distance <- Surveys$perp_distance %>% filter(Date >= FirstDate & Date <= LastDate)
  # grid fractions - selecting only the grid fractions in the transects for the time period selected
  GridFrac <- GridFrac[which((GridFrac$TransectID %in% Surveys$line_transect$TransectID) | (GridFrac$TransectID %in% Surveys$strip_transect$TransectID) | (GridFrac$TransectID %in% Surveys$uaoa$TransectID)), ]
  # covariates  - selecting only the grids in transects for the time period selected
  # and selecting only the time periods needed for the temporally variable covariates
  CovConsSurv <- CovConsSurv[which((CovConsSurv$GridID %in% GridFrac$GridID)),]
  CovTempSurv <- CovTempSurv[which((CovTempSurv[,"GridID", 1] %in% GridFrac$GridID)), , (FirstDateID - Lag):LastDateID]

  # set up data for small grids (only for those within the surveyed areas) and the large grid and genetic populations
  NSGrids <- length(CovConsSurv$GridID)
  GenPopID <- left_join(as.data.frame(CovConsSurv$GridID), as.data.frame(GenPopLookup), by = c("CovConsSurv$GridID" = "GridID"))$GENPOP_ID
  NGPops <- length(unique(GenPopID))

  # static predictors

  # set up list for scale parameters
  ScaleParamsX <- list()

  # ground water dependent ecosystems
  # swap categories 1 and 2 around so that the reference category is the most common
  # this means that 1 = intermittent and freshwater, and
  # 2 = brackish, saline or fluctuating salinity with intermittent connectivity,
  X_hhgde <- as.vector(CovConsSurv[,"hhgde"]) %>% as.character() %>% case_match("1"~2,"2"~1,"3"~3,"4"~4, "5"~5, "6"~6) %>% as.numeric() %>% as.factor()

  # elevation
  X_htele <- as.vector(CovConsSurv[,"htele"]) %>% scale()
  ScaleParamsX$htele <- c(X_htele %>% attr("scaled:center"), X_htele %>% attr("scaled:scale"))
  X_htele <- X_htele %>% as.vector()

  # slope
  X_htslo <- as.vector(CovConsSurv[,"htslo"]) %>% scale()
  ScaleParamsX$htslo <- c(X_htslo %>% attr("scaled:center"), X_htslo %>% attr("scaled:scale"))
  X_htslo <- X_htslo %>% as.vector()

  # terrain ruggedness index
  X_htrug <- as.vector(CovConsSurv[,"htrug"]) %>% scale()
  ScaleParamsX$htrug <- c(X_htrug %>% attr("scaled:center"), X_htrug %>% attr("scaled:scale"))
  X_htrug <- X_htrug %>% as.vector()

  # soil - do a PCA and use first two components to reduce domensions
  X_hscec <- as.vector(CovConsSurv[,"hscec"]) # cation exchange
  X_hspho <- as.vector(CovConsSurv[,"hspho"]) # phosphorous
  X_hsnit <- as.vector(CovConsSurv[,"hsnit"]) # nitrogen
  X_hswat <- as.vector(CovConsSurv[,"hswat"]) # available water capacity
  Soil <- tibble(hscec = X_hscec, hspho = X_hspho, hsnit = X_hsnit, hswat = X_hswat)
  # run PCA
  Soil_PCA <- prcomp(~hscec + hspho + hsnit + hswat, data = Soil, scale = TRUE)
  # individual PCA scores for first two components
  IndScores <- as_tibble(1:nrow(Soil)) %>% left_join(bind_cols(as_tibble(as.numeric(rownames(get_pca_ind(Soil_PCA)$coord))), as_tibble(get_pca_ind(Soil_PCA)$coord[,1:2])), by = c("value" = "value"))

  # create PC component values for each small grid
  X_hspc1 <- IndScores$Dim.1 %>% scale() # PC1
  ScaleParamsX$hspc1 <- c(X_hspc1 %>% attr("scaled:center"), X_hspc1 %>% attr("scaled:scale"))
  X_hspc1 <- X_hspc1 %>% as.vector()
  X_hspc2 <- IndScores$Dim.2 %>% scale() # PC2
  ScaleParamsX$hspc2 <- c(X_hspc2 %>% attr("scaled:center"), X_hspc2 %>% attr("scaled:scale"))
  X_hspc2 <- X_hspc2 %>% as.vector()

  # long-term precipitation
  X_hcltp <- as.vector(CovConsSurv[,"hcltp"]) %>% scale()
  ScaleParamsX$hcltp <- c(X_hcltp %>% attr("scaled:center"), X_hcltp %>% attr("scaled:scale"))
  X_hcltp <- X_hcltp %>% as.vector()

  # long-term temperature
  X_hcltt <- as.vector(CovConsSurv[,"hcltt"]) %>% scale()
  ScaleParamsX$hcltt <- c(X_hcltt %>% attr("scaled:center"), X_hcltt %>% attr("scaled:scale"))
  X_hcltt <- X_hcltt %>% as.vector()

  # compile into a tibble
  X <- tibble(hhgde = X_hhgde, htele = X_htele, htslo = X_htslo, htrug = X_htrug, hspc1 = X_hspc1, hspc2 = X_hspc2, hcltp = X_hcltp, hcltt = X_hcltt)
  # impute any missing values
  if (any(is.na(X))) {
      X <- complete(mice(X, m = 1)) %>% as_tibble()
  }

  # dynamic predictors

  # set up list for scale parameters
  ScaleParamsY <- list()

  # time
  Y_htime <- as.vector(CovTempSurv[,"htime",])

  # season
  Y_hseas <- as.vector(CovTempSurv[,"hseas",])

  # persistent green
  Y_hhpgr <- as.vector(CovTempSurv[,"hhpgr",]) %>% scale()
  ScaleParamsY$hhpgr <- c(Y_hhpgr %>% attr("scaled:center"), Y_hhpgr %>% attr("scaled:scale"))
  Y_hhpgr <- Y_hhpgr %>% as.vector()

  # persistent green buffer (2km)
  Y_hhpgr2km <- as.vector(CovTempSurv[,"hhpgr2km",]) %>% scale()
  ScaleParamsY$hhpgr2km <- c(Y_hhpgr2km %>% attr("scaled:center"), Y_hhpgr2km %>% attr("scaled:scale"))
  Y_hhpgr2km <- Y_hhpgr2km %>% as.vector()

  # koala habitat
  Y_hhkha <- as.vector(CovTempSurv[,"hhkha",]) %>% as.factor()

  # woody cover - using proportion of each cell that is woody + sparse woody
  Y_hhfwc <- as.vector(CovTempSurv[,"hhfwc_1",] + CovTempSurv[,"hhfwc_2",]) %>% scale()
  ScaleParamsY$hhfwc <- c(Y_hhfwc %>% attr("scaled:center"), Y_hhfwc %>% attr("scaled:scale"))
  Y_hhfwc <- Y_hhfwc %>% as.vector()

  # precipitation - centre on mean long term precipitation
  Y_hcpre <- as.vector(CovTempSurv[,"hcpre",]) - (rep(as.vector(CovConsSurv[,"hcltp"]), length(as.vector(CovTempSurv[,"hcpre",]))/length(as.vector(CovConsSurv[,"hcltp"]))) / 2)
  # remove season effects
  Y_hcpre[which(as.vector(CovTempSurv[,"hseas",]) == 1)] <- Y_hcpre[which(as.vector(CovTempSurv[,"hseas",]) == 1)] - mean(Y_hcpre[which(as.vector(CovTempSurv[,"hseas",]) == 1)])
  Y_hcpre[which(as.vector(CovTempSurv[,"hseas",]) == 0)] <- Y_hcpre[which(as.vector(CovTempSurv[,"hseas",]) == 0)] - mean(Y_hcpre[which(as.vector(CovTempSurv[,"hseas",]) == 0)])
  # centre and scale
  Y_hcpre <- Y_hcpre %>% scale()
  ScaleParamsY$hcpre <- c(Y_hcpre %>% attr("scaled:center"), Y_hcpre %>% attr("scaled:scale"))
  Y_hcpre <- Y_hcpre %>% as.vector()

  # mean temperature - centre on long-term mean mean temperature
  Y_hctmn <- as.vector(CovTempSurv[,"hctmn",]) - (rep(as.vector(CovConsSurv[,"hcltt"]), length(as.vector(CovTempSurv[,"hctmn",]))/length(as.vector(CovConsSurv[,"hcltt"]))))
  # remove season effects
  Y_hctmn[which(as.vector(CovTempSurv[,"hseas",]) == 1)] <- Y_hctmn[which(as.vector(CovTempSurv[,"hseas",]) == 1)] - mean(Y_hctmn[which(as.vector(CovTempSurv[,"hseas",]) == 1)])
  Y_hctmn[which(as.vector(CovTempSurv[,"hseas",]) == 0)] <- Y_hctmn[which(as.vector(CovTempSurv[,"hseas",]) == 0)] - mean(Y_hctmn[which(as.vector(CovTempSurv[,"hseas",]) == 0)])
  # centre and scale
  Y_hctmn <- Y_hctmn %>% scale()
  ScaleParamsY$hctmn <- c(Y_hctmn %>% attr("scaled:center"), Y_hctmn %>% attr("scaled:scale"))
  Y_hctmn <- Y_hctmn %>% as.vector()

  # max temperature
  Y_hctma <- as.vector(CovTempSurv[,"hctma",]) - (rep(as.vector(CovConsSurv[,"hcltt"]), length(as.vector(CovTempSurv[,"hctma",]))/length(as.vector(CovConsSurv[,"hcltt"]))))
  # remove season effects
  Y_hctma[which(as.vector(CovTempSurv[,"hseas",]) == 1)] <- Y_hctma[which(as.vector(CovTempSurv[,"hseas",]) == 1)] - mean(Y_hctma[which(as.vector(CovTempSurv[,"hseas",]) == 1)])
  Y_hctma[which(as.vector(CovTempSurv[,"hseas",]) == 0)] <- Y_hctma[which(as.vector(CovTempSurv[,"hseas",]) == 0)] - mean(Y_hctma[which(as.vector(CovTempSurv[,"hseas",]) == 0)])
  # centre and scale
  Y_hctma <- Y_hctma %>% scale()
  ScaleParamsY$hctma <- c(Y_hctma %>% attr("scaled:center"), Y_hctma %>% attr("scaled:scale"))
  Y_hctma <- Y_hctma %>% as.vector()

  # land-use
  # reclassfy so that: 1 = natural, 2 = production from natural, 3 = intensive, 4 = other
  Y_htlus <- as.vector(CovTempSurv[,"htlus",]) %>% as.character() %>% case_match("1"~1,"2"~2,"3"~4,"4"~4, "5"~3, "6"~4, "7"~4) %>% as.numeric() %>% as.factor()

  # intensive land-use buffer (2km)
  Y_htilu2km <- as.vector(CovTempSurv[,"htilu2km",]) %>% scale()
  ScaleParamsY$htilu2km <- c(Y_htilu2km %>% attr("scaled:center"), Y_htilu2km %>% attr("scaled:scale"))
  Y_htilu2km <- Y_htilu2km %>% as.vector()

  # lot size buffer (2km)
  Y_htpls2km <- as.vector(CovTempSurv[,"htpls2km",]) %>% scale()
  ScaleParamsY$htpls2km <- c(Y_htpls2km %>% attr("scaled:center"), Y_htpls2km %>% attr("scaled:scale"))
  Y_htpls2km <- Y_htpls2km %>% as.vector()

  # compile into a tibble
  Y_temp <- tibble(htime = Y_htime, hseas = Y_hseas, hhpgr = Y_hhpgr, hhpgr2km = Y_hhpgr2km, hhkha = Y_hhkha, hhfwc = Y_hhfwc, hcpre = Y_hcpre, hctmn = Y_hctmn, hctma = Y_hctma, htlus = Y_htlus, htilu2km = Y_htilu2km, htpls2km = Y_htpls2km)

  # impute any missing values
  if (any(is.na(Y_temp))) {
      Y_temp <- complete(mice(Y_temp, m = 1)) %>% as_tibble()
  }

  # combine X and Y to check for collinearity
  XYData <- bind_cols(do.call(rbind, replicate(dim(CovTempSurv)[3], X, simplify = FALSE)), Y_temp)

  # check for collinearity among continuous predictors
  CorrXY <- cor(XYData %>% dplyr::select(-hhgde, -htime, -hseas, -hhkha, -htlus), use = "complete.obs", method = "spearman")

  # get the design matrix and remove collinear variables for X
  X <- model.matrix(as.formula(paste0("~ ", paste(StaticVars, collapse= " + "))), model.frame(as.formula(paste0("~ ", paste(StaticVars, collapse= " + "))), as.data.frame(X), na.action = "na.pass")) %>% as.data.frame()

  # remove the intercept term
  X <- X[, 2:ncol(X)] %>% as.data.frame()

  # get number of variables in X
  NX <- ncol(X)

  # get the design matrix and remove collinear variables for X
  # remove persistent green, max temperature, and intensive land-use buffer (based on correlations > 0.7)
  # persistent green correlated with persistent green buffer (r = 0.72)
  # max temperature correlated with mean temperature (r = 0.91)
  # intensive land-use buffer correlated with lot size (r = -0.73)
  Y_temp <- model.matrix(as.formula(paste0("~ ", paste(DynamicVars, collapse= " + "))), model.frame(as.formula(paste0("~ ", paste(DynamicVars, collapse= " + "))), as.data.frame(Y_temp), na.action = "na.pass"))

  # remove the intercept term
  Y_temp <- Y_temp[, 2:ncol(Y_temp)] %>% as.data.frame()

  # get number of variables in Y
  NY <- ncol(Y_temp)

  # create 3D matrix of covariates for Y
  Y <- array(NA, dim = c(dim(CovTempSurv)[1], dim(CovTempSurv)[3], NY))
  for (i in 1:NY) {
      Y[,,i] <- matrix(Y_temp[,i], nrow = dim(CovTempSurv)[1], ncol = dim(CovTempSurv)[3])
  }

  # predictors for detectability (note at present these are static spatial variables)

  # canopy height
  Z_hhcht <- as.vector(CovConsSurv[,"hhcht"]) %>% scale() %>% as.vector()

  # understorey fraction (up to 5m)
  Z_hhunf <- as.vector(CovConsSurv[,"hhunf"]) %>% scale() %>% as.vector()

  # interaction between canopy height and understorey fraction
  Z_hhchtunf <- (as.vector(CovConsSurv[,"hhcht"]) * as.vector(CovConsSurv[,"hhunf"])) %>% scale() %>% as.vector()

  # get lookup table to convert original grid IDs to IDs ordered from 1:number of small grids in the survey areas (to match the order in CovConsSurv and CovTempSurv)
  SGridIDsLookUp <- CovConsSurv$GridID %>% as_tibble() %>% mutate(GridID = value, OrderedID = seq_len(nrow(.))) %>% select(GridID, OrderedID)

  # strip transect data

  # join grid fractions data to grid fractions data for strip transects
  StripNestedGroupFrac <- GridFrac %>% filter(transect == "strip_transect") %>% as_tibble() %>% group_by(TransectID) %>% nest()
  StripJoinGroupFrac <- as_tibble(Surveys$strip_transect) %>% left_join(StripNestedGroupFrac, by = c("TransectID" = "TransectID"), keep = TRUE)
  # remove any strip transects with no grid fraction data
  StripJoinGroupFrac <- StripJoinGroupFrac %>% filter(!is.na(TransectID.y))

  # number of strip transects
  NStrips <- nrow(StripJoinGroupFrac)

  # get small grid IDs for each strip transect
  SGridIDsStrip <- StripJoinGroupFrac$data %>% map(~ .x$GridID) %>% unlist() %>% as_tibble() %>% left_join(SGridIDsLookUp, by = c("value" = "GridID"))
  SGridIDsStrip <- SGridIDsStrip$OrderedID %>% as.vector()

  # get small grid fractions for each strip transect
  SGridFracsStrip <- StripJoinGroupFrac$data %>% map(~ .x$fraction) %>% unlist() %>% as.vector()

  # get the start and end indexes for the small grids for each strip transect
  StripNumSGrids <- StripJoinGroupFrac$data %>% map(~ nrow(.x)) %>% unlist() %>% as.vector()
  SGridsStartStrip <- cumsum(StripNumSGrids) - StripNumSGrids + 1
  SGridsEndStrip <-cumsum(StripNumSGrids)

  # create matrices for small grid IDs and fractions for each strip transect
  IDStrip <- matrix(0, nrow = NStrips, ncol = max(StripNumSGrids))
  FracStrip <- matrix(0, nrow = NStrips, ncol = max(StripNumSGrids))
  for (i in 1:NStrips) {
          for (j in SGridsStartStrip[i]:SGridsEndStrip[i]) {
              IDStrip[i, j - SGridsStartStrip[i] + 1] <- SGridIDsStrip[j]
              FracStrip[i, j - SGridsStartStrip[i] + 1] <- SGridFracsStrip[j]
          }
  }

  # get maximum number of small grids in any strip transect
  NMaxSGridsStrip <- max(StripNumSGrids)

  # get the number of observers for each strip transect
  NumObsStrip <- StripJoinGroupFrac$Number_Observers %>% as.vector()

  # get the area of each strip transect
  AreaStrip <- StripJoinGroupFrac$TArea %>% as.vector()
  AreaStrip[which(AreaStrip == 0)] <- 4 #EDIT OUT WHEN FIXED

  # get the time step ID for each strip transect
  TimeIDStrip <- StripJoinGroupFrac$TimePeriodID %>% as.vector()
  TimeIDStrip <- TimeIDStrip - FirstDateID + 1 + Lag

  # get the koala count for each strip transect
  CntStrip <- StripJoinGroupFrac$Number_Sightings %>% as.vector()

  # get the detectability covariates for each strip
  Z_Strip_hhcht_Mat <- matrix(0, nrow = NStrips, ncol = max(StripNumSGrids))
  Z_Strip_hhunf_Mat <- matrix(0, nrow = NStrips, ncol = max(StripNumSGrids))
  Z_Strip_hhchtunf_Mat <- matrix(0, nrow = NStrips, ncol = max(StripNumSGrids))
  Z_Strip_hhcht <- rep(NA, NStrips)
  Z_Strip_hhunf <- rep(NA, NStrips)
  Z_Strip_hhchtunf <- rep(NA, NStrips)
  for (i in 1:NStrips) {
          for (j in SGridsStartStrip[i]:SGridsEndStrip[i]) {
              Z_Strip_hhcht_Mat[i, j - SGridsStartStrip[i] + 1] <- Z_hhcht[SGridIDsStrip[j]]
              Z_Strip_hhunf_Mat[i, j - SGridsStartStrip[i] + 1] <- Z_hhunf[SGridIDsStrip[j]]
              Z_Strip_hhchtunf_Mat[i, j - SGridsStartStrip[i] + 1] <- Z_hhchtunf[SGridIDsStrip[j]]
          }
          Z_Strip_hhcht[i] <- sum(Z_Strip_hhcht_Mat[i, ] * FracStrip[i, ], na.rm = TRUE)
          Z_Strip_hhunf[i] <- sum(Z_Strip_hhunf_Mat[i, ] * FracStrip[i, ], na.rm = TRUE)
          Z_Strip_hhchtunf[i] <- sum(Z_Strip_hhchtunf_Mat[i, ] * FracStrip[i, ], na.rm = TRUE)
  }

  # compile into a tibble
  Z_Strip <- tibble(hhcht = Z_Strip_hhcht, hhunf = Z_Strip_hhunf, hhchtunf = Z_Strip_hhchtunf)
  # impute any missing values
  if (any(is.na(Z_Strip))) {
      Z_Strip <- complete(mice(Z_Strip, m = 1)) %>% as_tibble()
  }

  # get the design matrix
  Z_Strip <- model.matrix(~ hhcht + hhunf + hhchtunf, model.frame(~ hhcht + hhunf + hhchtunf, as.data.frame(Z_Strip), na.action = "na.pass"))

  # all of area search data

  # join grid fractions data to grid fractions data for all of areas searches
  AoANestedGroupFrac <- GridFrac %>% filter(transect == "uaoa") %>% as_tibble() %>% group_by(TransectID) %>% nest()
  AoAJoinGroupFrac <- as_tibble(Surveys$uaoa) %>% left_join(AoANestedGroupFrac, by = c("TransectID" = "TransectID"), keep = TRUE)
  # remove any all of areas searches with no grid fraction data
  AoAJoinGroupFrac <- AoAJoinGroupFrac %>% filter(!is.na(TransectID.y))

  # number of all of area searches
  NAoAs <- nrow(AoAJoinGroupFrac)

  # get small grid IDs for each all of area search
  SGridIDsAoA <- AoAJoinGroupFrac$data %>% map(~ .x$GridID) %>% unlist() %>% as_tibble() %>% left_join(SGridIDsLookUp, by = c("value" = "GridID"))
  SGridIDsAoA <- SGridIDsAoA$OrderedID %>% as.vector()

  # get small grid fractions for each all of area search
  SGridFracsAoA <- AoAJoinGroupFrac$data %>% map(~ .x$fraction) %>% unlist() %>% as.vector()

  # get the start and end indexes for the small grids for each all of area search
  AoANumSGrids <- AoAJoinGroupFrac$data %>% map(~ nrow(.x)) %>% unlist() %>% as.vector()
  SGridsStartAoA <- cumsum(AoANumSGrids) - AoANumSGrids + 1
  SGridsEndAoA <-cumsum(AoANumSGrids)

  # create matrices for small grid IDs and fractions for each strip transect
  IDAoA <- matrix(0, nrow = NAoAs, ncol = max(AoANumSGrids))
  FracAoA <- matrix(0, nrow = NAoAs, ncol = max(AoANumSGrids))
  for (i in 1:NAoAs) {
          for (j in SGridsStartAoA[i]:SGridsEndAoA[i]) {
              IDAoA[i, j - SGridsStartAoA[i] + 1] <- SGridIDsAoA[j]
              FracAoA[i, j - SGridsStartAoA[i] + 1] <- SGridFracsAoA[j]
          }
  }

  # get maximum number of small grids in any all of area search
  NMaxSGridsAoA <- max(AoANumSGrids)

  # get the number of observers for each all of area search
  NumObsAoA <- AoAJoinGroupFrac$Number_Observers %>% as.vector()

  # get the area of each all of area search
  AreaAoA <- AoAJoinGroupFrac$TArea %>% as.vector()

  # get the time step ID for each all of area search
  TimeIDAoA <- AoAJoinGroupFrac$TimePeriodID %>% as.vector()
  TimeIDAoA <- TimeIDAoA - FirstDateID + 1 + Lag

  # get the koala count for each all of areas search
  CntAoA <- AoAJoinGroupFrac$Number_Sightings %>% as.vector()

  # get the detectability covariates for each AoA
  Z_AoA_hhcht_Mat <- matrix(0, nrow = NAoAs, ncol = max(AoANumSGrids))
  Z_AoA_hhunf_Mat <- matrix(0, nrow = NAoAs, ncol = max(AoANumSGrids))
  Z_AoA_hhchtunf_Mat <- matrix(0, nrow = NAoAs, ncol = max(AoANumSGrids))
  Z_AoA_hhcht <- rep(NA, NAoAs)
  Z_AoA_hhunf <- rep(NA, NAoAs)
  Z_AoA_hhchtunf <- rep(NA, NAoAs)
  for (i in 1:NAoAs) {
          for (j in SGridsStartAoA[i]:SGridsEndAoA[i]) {
              Z_AoA_hhcht_Mat[i, j - SGridsStartAoA[i] + 1] <- Z_hhcht[SGridIDsAoA[j]]
              Z_AoA_hhunf_Mat[i, j - SGridsStartAoA[i] + 1] <- Z_hhunf[SGridIDsAoA[j]]
              Z_AoA_hhchtunf_Mat[i, j - SGridsStartAoA[i] + 1] <- Z_hhchtunf[SGridIDsAoA[j]]
          }
          Z_AoA_hhcht[i] <- sum(Z_AoA_hhcht_Mat[i, ] * FracAoA[i, ], na.rm = TRUE)
          Z_AoA_hhunf[i] <- sum(Z_AoA_hhunf_Mat[i, ] * FracAoA[i, ], na.rm = TRUE)
          Z_AoA_hhchtunf[i] <- sum(Z_AoA_hhchtunf_Mat[i, ] * FracAoA[i, ], na.rm = TRUE)
  }

  # compile into a tibble
  Z_AoA <- tibble(hhcht = Z_AoA_hhcht, hhunf = Z_AoA_hhunf, hhchtunf = Z_AoA_hhchtunf)
  # impute any missing values
  if (any(is.na(Z_AoA))) {
      Z_AoA <- complete(mice(Z_AoA, m = 1)) %>% as_tibble()
  }

  # get the design matrix
  Z_AoA <- model.matrix(~ hhcht + hhunf + hhchtunf, model.frame(~ hhcht + hhunf + hhchtunf, as.data.frame(Z_AoA), na.action = "na.pass"))

  # line transect data

  # join grid fractions data to grid fractions data for line transects
  LineNestedGroupFrac <- GridFrac %>% filter(transect == "line_transect") %>% as_tibble() %>% group_by(TransectID) %>% nest()
  LineJoinGroupFrac <- as_tibble(Surveys$line_transect) %>% left_join(LineNestedGroupFrac, by = c("TransectID" = "TransectID"), keep = TRUE)
  # remove any line transects with no grid fraction data
  LineJoinGroupFrac <- LineJoinGroupFrac %>% filter(!is.na(TransectID.y))

  # number of line transects
  NLines <- nrow(LineJoinGroupFrac)

  # get small grid IDs for each line transect
  SGridIDsLine <- LineJoinGroupFrac$data %>% map(~ .x$GridID) %>% unlist() %>% as_tibble() %>% left_join(SGridIDsLookUp, by = c("value" = "GridID"))
  SGridIDsLine <- SGridIDsLine$OrderedID %>% as.vector()

  # get small grid fractions for each line transect
  SGridFracsLine <- LineJoinGroupFrac$data %>% map(~ .x$fraction) %>% unlist() %>% as.vector()

  # get the start and end indices for the small grids for each line transect
  LineNumSGrids <- LineJoinGroupFrac$data %>% map(~ nrow(.x)) %>% unlist() %>% as.vector()
  SGridsStartLine <- cumsum(LineNumSGrids) - LineNumSGrids + 1
  SGridsEndLine <-cumsum(LineNumSGrids)

  # create matrices for small grid IDs and fractions for each strip transect
  IDLine <- matrix(0, nrow = NLines, ncol = max(LineNumSGrids))
  FracLine <- matrix(0, nrow = NLines, ncol = max(LineNumSGrids))
  for (i in 1:NLines) {
          for (j in SGridsStartLine[i]:SGridsEndLine[i]) {
              IDLine[i, j - SGridsStartLine[i] + 1] <- SGridIDsLine[j]
              FracLine[i, j - SGridsStartLine[i] + 1] <- SGridFracsLine[j]
          }
  }

  # get maximum number of small grids in any line transect
  NMaxSGridsLine <- max(LineNumSGrids)

  # get the length of each line transect
  LengthLine <- LineJoinGroupFrac$Tlength %>% as.vector()

  # get the time step ID for each line transect
  TimeIDLine <- LineJoinGroupFrac$TimePeriodID %>% as.vector()
  TimeIDLine <- TimeIDLine - FirstDateID + 1 + Lag

  # get the koala count for each line transect
  CntLine <- LineJoinGroupFrac$Number_Sightings %>% as.vector()

  # get the detectability covariates for each Line
  Z_Line_hhcht_Mat <- matrix(0, nrow = NLines, ncol = max(LineNumSGrids))
  Z_Line_hhunf_Mat <- matrix(0, nrow = NLines, ncol = max(LineNumSGrids))
  Z_Line_hhchtunf_Mat <- matrix(0, nrow = NLines, ncol = max(LineNumSGrids))
  Z_Line_hhcht <- rep(NA, NLines)
  Z_Line_hhunf <- rep(NA, NLines)
  Z_Line_hhchtunf <- rep(NA, NLines)
  for (i in 1:NLines) {
          for (j in SGridsStartLine[i]:SGridsEndLine[i]) {
              Z_Line_hhcht_Mat[i, j - SGridsStartLine[i] + 1] <- Z_hhcht[SGridIDsLine[j]]
              Z_Line_hhunf_Mat[i, j - SGridsStartLine[i] + 1] <- Z_hhunf[SGridIDsLine[j]]
              Z_Line_hhchtunf_Mat[i, j - SGridsStartLine[i] + 1] <- Z_hhchtunf[SGridIDsLine[j]]
          }
          Z_Line_hhcht[i] <- sum(Z_Line_hhcht_Mat[i, ] * FracLine[i, ], na.rm = TRUE)
          Z_Line_hhunf[i] <- sum(Z_Line_hhunf_Mat[i, ] * FracLine[i, ], na.rm = TRUE)
          Z_Line_hhchtunf[i] <- sum(Z_Line_hhchtunf_Mat[i, ] * FracLine[i, ], na.rm = TRUE)
  }

  # compile into a tibble
  Z_Line <- tibble(hhcht = Z_Line_hhcht, hhunf = Z_Line_hhunf, hhchtunf = Z_Line_hhchtunf)
  # impute any missing values
  if (any(is.na(Z_Line))) {
      Z_Line <- complete(mice(Z_Line, m = 1)) %>% as_tibble()
  }

  # get the design matrix
  Z_Line <- model.matrix(~ hhcht + hhunf + hhchtunf, model.frame(~ hhcht + hhunf + hhchtunf, as.data.frame(Z_Line), na.action = "na.pass"))

  # get number of variables in Z
  NZ <- ncol(Z_Line)

  # get transect IDs for each perpendicular distance
  PDLineIDs <- match(Surveys$perp_distance$TransectID, LineJoinGroupFrac$TransectID.x)

  # get the number of perependicular distances - remove PDs that don't have a corresponding transect
  NPDists <- nrow(Surveys$perp_distance[which(!is.na(PDLineIDs)), ])

  # get perpendicular distances
  PDists <- Surveys$perp_distance$Perp_Dist[which(!is.na(PDLineIDs))] %>% as.vector()

  # remove missing Line IDs
  PDLineIDs <- PDLineIDs[which(!is.na(PDLineIDs))]

  # reset first and last date IDs
  FirstDateID_Orig <- FirstDateID
  LastDateID_Orig <- LastDateID
  LastDateID <- LastDateID - FirstDateID + 1 + Lag
  FirstDateID <- Lag + 1

  # set up data for spatial CAR process (1st order with equal weights) for genetic populations
  # assumes the following first-order spatial connections among genetic populations:
  # 1 (Noosa) -> 2
  # 2 (North Coast) -> 1, 3, 5, 6
  # 3 (North West) -> 2, 5, 4
  # 4 (Toowoomba)-> 3, 5
  # 5 (South West) -> 2, 3, 4, 6, 8
  # 6 (Koala Coast) -> 2, 5, 8
  # 7 (Minjerribah) -> disconnected
  # 8 (Gold Coast) -> 5, 6
  AdjS <- c()
  WeightsAdjS <- c()
  NumAdjS <- c()
  NumAdjS[1] <- 1
  NumAdjS[2] <- 4
  NumAdjS[3] <- 3
  NumAdjS[4] <- 2
  NumAdjS[5] <- 5
  NumAdjS[6] <- 3
  NumAdjS[7] <- 0
  NumAdjS[8] <- 2
  AdjS[1] <- 2
  AdjS[2] <- 1
  AdjS[3] <- 3
  AdjS[4] <- 5
  AdjS[5] <- 6
  AdjS[6] <- 2
  AdjS[7] <- 5
  AdjS[8] <- 4
  AdjS[9] <- 3
  AdjS[10] <- 5
  AdjS[11] <- 2
  AdjS[12] <- 3
  AdjS[13] <- 4
  AdjS[14] <- 6
  AdjS[15] <- 8
  AdjS[16] <- 2
  AdjS[17] <- 5
  AdjS[18] <- 8
  AdjS[19] <- 5
  AdjS[20] <- 6
  for (i in 1:20) {
      WeightsAdjS[i] <- 1
  }
  NGPopsAdjs <- length(AdjS)

  # set up data for temporal CAR process (1st order with equal weights)
  # annual time steps
  AdjT <- c()
  WeightsAdjT <- c()
  NumAdjT <- c()
  WeightsAdjT[1] <- 1
  AdjT[1] <- 2
  NumAdjT[1] <- 1
  for(i in 2:(floor(((LastDateID - FirstDateID + 1 - 1) / 2)))) {
      WeightsAdjT[2 + (i - 2) * 2] <- 1
      AdjT[2 + (i - 2) * 2] <- i - 1
      WeightsAdjT[3 + (i - 2) * 2] <- 1
      AdjT[3 + (i - 2) * 2] <- i + 1
      NumAdjT[i] <- 2
  }
  for(i in (floor(((LastDateID - FirstDateID + 1 - 1) / 2)) + 1):(floor(((LastDateID - FirstDateID + 1 - 1) / 2)) + 1)) {
      WeightsAdjT[(i - 2) * 2 + 2] <- 1
      AdjT[(i - 2) * 2 + 2] <- i - 1
      NumAdjT[i] <- 1
  }
  NTime <- length(NumAdjT)
  NTimeAdjs <- length(AdjT)

  # set up data for temporal CAR process (2nd order with weights acording to Breslow & Clayton (1993) - thin plate spline)
  # annual time steps
  AdjT2 <- c()
  WeightsAdjT2 <- c()
  NumAdjT2 <- c()
  WeightsAdjT2[1] <- 2
  WeightsAdjT2[2] <- -1
  AdjT2[1] <- 2
  AdjT2[2] <- 3
  NumAdjT2[1] <- 2
  WeightsAdjT2[3] <- 2
  WeightsAdjT2[4] <- 4
  WeightsAdjT2[5] <- -1
  AdjT2[3] <- 1
  AdjT2[4] <- 3
  AdjT2[5] <- 4
  NumAdjT2[2] <- 3
  for(i in 3:(floor(((LastDateID - FirstDateID + 1 - 1) / 2)) - 1)) {
      WeightsAdjT2[6 + (i - 3) * 4] <- -1
      WeightsAdjT2[7 + (i - 3) * 4] <- 4
      WeightsAdjT2[8 + (i - 3) * 4] <- 4
      WeightsAdjT2[9 + (i - 3) * 4] <- -1
      AdjT2[6 + (i - 3) * 4] <- i - 2
      AdjT2[7 + (i - 3) * 4] <- i - 1
      AdjT2[8 + (i - 3) * 4] <- i + 1
      AdjT2[9 + (i - 3) * 4] <- i + 2
      NumAdjT2[i] <- 4
  }
  for(i in (floor(((LastDateID - FirstDateID + 1 - 1) / 2))):(floor(((LastDateID - FirstDateID + 1 - 1) / 2)))) {  	WeightsAdjT2[6 + (i - 3) * 4] <- -1
      WeightsAdjT2[7 + (i - 3) * 4] <- 4
      WeightsAdjT2[8 + (i - 3) * 4] <- 2
      AdjT2[6 + (i - 3) * 4] <- i - 2
      AdjT2[7 + (i - 3) * 4] <- i - 1
      AdjT2[8 + (i - 3) * 4] <- i + 1
      NumAdjT2[i] <- 3
  }
  for(i in (floor(((LastDateID - FirstDateID + 1 - 1) / 2)) + 1):(floor(((LastDateID - FirstDateID + 1 - 1) / 2)) + 1)) {
      WeightsAdjT2[6 + (i - 3 - 1) * 4 + 3] <- -1
      WeightsAdjT2[7 + (i - 3 - 1) * 4 + 3] <- 2
      AdjT2[6 + (i - 3 - 1) * 4 + 3] <- i - 2
      AdjT2[7 + (i - 3 - 1) * 4 + 3] <- i - 1
      NumAdjT2[i] <- 2
  }
  NTime2 <- length(NumAdjT2)
  NTimeAdjs2 <- length(AdjT2)

  # set up nimble constants and data inputs for model with temporal CAR process order = 1
  NimbleConsts <- list(NGPops = NGPops, NGPopsAdjs = NGPopsAdjs, AdjS = AdjS, WeightsAdjS = WeightsAdjS, NumAdjS = NumAdjS, NTime = NTime, NTimeAdjs = NTimeAdjs, AdjT = AdjT, WeightsAdjT = WeightsAdjT, NumAdjT = NumAdjT, NTime2 = NTime2, NTimeAdjs2 = NTimeAdjs2, AdjT2 = AdjT2, WeightsAdjT2 = WeightsAdjT2, NumAdjT2 = NumAdjT2, NSGrids = NSGrids, GenPopID = GenPopID, FirstDateID = FirstDateID, LastDateID = LastDateID, NX = NX, NY = NY, NStrips = NStrips, SGridsStartStrip = SGridsStartStrip, SGridsEndStrip = SGridsEndStrip, SGridIDsStrip = SGridIDsStrip, SGridFracsStrip = SGridFracsStrip, AreaStrip = AreaStrip, TimeIDStrip = TimeIDStrip, NAoAs = NAoAs, SGridsStartAoA = SGridsStartAoA, SGridsEndAoA = SGridsEndAoA, SGridIDsAoA = SGridIDsAoA, SGridFracsAoA = SGridFracsAoA, AreaAoA = AreaAoA, TimeIDAoA = TimeIDAoA, NLines = NLines, SGridsStartLine = SGridsStartLine, SGridsEndLine = SGridsEndLine, SGridIDsLine = SGridIDsLine, SGridFracsLine = SGridFracsLine, LengthLine = LengthLine, TimeIDLine = TimeIDLine, PI = pi, NMaxSGridsAoA = NMaxSGridsAoA, NMaxSGridsStrip = NMaxSGridsStrip, NMaxSGridsLine = NMaxSGridsLine, NZ = NZ, NPDists = NPDists, PDLineIDs = PDLineIDs, Order = Order, Lag = Lag, VarTrend = VarTrend)
  NimbleData <- list(X = X, Y = Y, Z_Strip = Z_Strip, Z_AoA = Z_AoA, Z_Line = Z_Line, CntStrip = CntStrip, CntAoA = CntAoA, PDists = PDists, CntLine = CntLine)

  return(list(Constants = NimbleConsts, Data = NimbleData, NamesX = names(X), NamesY = names(Y_temp), ScalingX = ScaleParamsX, ScalingY = ScaleParamsY, FirstDate = FirstDate, LastDate = LastDate, FirstDateID_Orig = FirstDateID_Orig, LastDateID_Orig = LastDateID_Orig, SoilPCA = Soil_PCA, CorrXY = CorrXY, Order = Order, Lag = Lag, VarTrend = VarTrend, StaticVars = StaticVars, DynamicVars = DynamicVars))
}

# function to get data for a particular year for a model to make predictions
get_prediction_data <- function(Year, NimbleData, PredData, RainForestMask = TRUE) {

  # reclassify genetic population IDs to three "populations"
  # 1 (Noosa) -> 1
  # 2 (North Coast) -> 1
  # 3 (North West) -> 2
  # 4 (Toowoomba)-> 2
  # 5 (South West) -> 2
  # 6 (Koala Coast) -> 3
  # 7 (Minjerribah) -> 3
  # 8 (Gold Coast) -> 3
  PredData$GenPopLookup <- PredData$GenPopLookup %>% as_tibble() %>% mutate(GENPOP_ID = case_match(GENPOP_ID, 1~1, 2~1, 3~2, 4~2, 5~2, 6~3, 7~3, 8~3) %>% as.integer())

  # get first and last date IDs
  # check if the first date is in the data, use it, otherwise use the second date - to avoid issues at the start or end
  if (date(paste0(Year, "-04-01")) %in% PredData$DateIntervals$start_date) {
    FirstDateID <- filter(PredData$DateIntervals, start_date == date(paste0(Year, "-04-01"))) %>% pull(TimePeriodID)
  } else {
    FirstDateID <- filter(PredData$DateIntervals, start_date == date(paste0(Year, "-10-01"))) %>% pull(TimePeriodID)
  }
  # check if the second date is in the data, use it, otherwise use the first date - to avoid issues at the start or end
  if (date(paste0(Year, "-10-01")) %in% PredData$DateIntervals$start_date) {
    LastDateID <- filter(PredData$DateIntervals, start_date == date(paste0(Year, "-10-01"))) %>% pull(TimePeriodID)
  } else {
    LastDateID <- filter(PredData$DateIntervals, start_date == date(paste0(Year, "-04-01"))) %>% pull(TimePeriodID)
  }

  # add a season variable to the temporally variable covariates
  # 1 = breeding season (summer - October to March) and 0 = non-breeding season (winter - April to September)
  Seas <- ifelse(month(PredData$DateIntervals$start_date) == 10, 1, 0)
  Seas <- t(matrix(rep(Seas, dim(PredData$CovTemp)[1]), ncol = dim(PredData$CovTemp)[1], nrow = dim(PredData$CovTemp)[3]))
  PredData$CovTemp <- abind(PredData$CovTemp, Seas, along = 2)
  dimnames(PredData$CovTemp)[[2]][(dim(PredData$CovTemp)[[2]])] <- "hseas"

  # add a time step variable to the temporally variable covariates
  Time <- 1:dim(PredData$CovTemp)[3]
  Time <- t(matrix(rep(Time, dim(PredData$CovTemp)[1]), ncol = dim(PredData$CovTemp)[1], nrow = dim(PredData$CovTemp)[3]))
  PredData$CovTemp <- abind(PredData$CovTemp, Time, along = 2)
  dimnames(PredData$CovTemp)[[2]][(dim(PredData$CovTemp)[[2]])] <- "htime"

  # Subset temporally variable data for the specified year
  PredData$CovTemp <- PredData$CovTemp[ , , (FirstDateID - NimbleData$Constants$Lag):LastDateID]

  # set up small grid IDs
  SGridID = PredData$CovCons[,"GridID"]
  GenPopID <- left_join(as.data.frame(PredData$CovCons$GridID), as.data.frame(PredData$GenPopLookup), by = c("PredData$CovCons$GridID" = "GridID"))$GENPOP_ID

  # static predictors

  # ground water dependent ecosystems
  # swap categories 1 and 2 around so that the reference category is the most common
  # this means that 1 = intermittent and freshwater, and
  # 2 = brackish, saline or fluctuating salinity with intermittent connectivity,
  X_hhgde <- as.vector(PredData$CovCons[,"hhgde"]) %>% as.character() %>% case_match("1"~2,"2"~1,"3"~3,"4"~4, "5"~5, "6"~6) %>% as.numeric() %>% as.factor()

  # elevation
  X_htele <- as.vector(PredData$CovCons[,"htele"]) %>% scale(center = NimbleData$ScalingX$htele[1], scale = NimbleData$ScalingX$htele[2])
  X_htele <- X_htele %>% as.vector()

  # slope
  X_htslo <- as.vector(PredData$CovCons[,"htslo"]) %>% scale(center = NimbleData$ScalingX$htslo[1], scale = NimbleData$ScalingX$htslo[2])
  X_htslo <- X_htslo %>% as.vector()

  # terrain ruggedness index
  X_htrug <- as.vector(PredData$CovCons[,"htrug"]) %>% scale(center = NimbleData$ScalingX$htrug[1], scale = NimbleData$ScalingX$htrug[2])
  X_htrug <- X_htrug %>% as.vector()

  # soil - do a PCA and use first two components to reduce domensions
  X_hscec <- as.vector(PredData$CovCons[,"hscec"]) # cation exchange
  X_hspho <- as.vector(PredData$CovCons[,"hspho"]) # phosphorous
  X_hsnit <- as.vector(PredData$CovCons[,"hsnit"]) # nitrogen
  X_hswat <- as.vector(PredData$CovCons[,"hswat"]) # available water capacity
  Soil <- tibble(hscec = X_hscec, hspho = X_hspho, hsnit = X_hsnit, hswat = X_hswat)
  # predict from soil PCA
  PredPCA <- as_tibble(predict(NimbleData$SoilPCA, newdata = Soil))

  # create PC component values for each small grid
  X_hspc1 <- PredPCA$PC1 %>% scale(center = NimbleData$ScalingX$hspc1[1], scale = NimbleData$ScalingX$hspc1[2]) # PC1
  X_hspc1 <- X_hspc1 %>% as.vector()
  X_hspc2 <- PredPCA$PC2 %>% scale(center = NimbleData$ScalingX$hspc2[1], scale = NimbleData$ScalingX$hspc2[2]) # PC2
  X_hspc2 <- X_hspc2 %>% as.vector()

  # long-term precipitation
  X_hcltp <- as.vector(PredData$CovCons[,"hcltp"]) %>% scale(center = NimbleData$ScalingX$hcltp[1], scale = NimbleData$ScalingX$hcltp[2])
  X_hcltp <- X_hcltp %>% as.vector()

  # long-term temperature
  X_hcltt <- as.vector(PredData$CovCons[,"hcltt"]) %>% scale(center = NimbleData$ScalingX$hcltt[1], scale = NimbleData$ScalingX$hcltt[2])
  X_hcltt <- X_hcltt %>% as.vector()

  # compile into a tibble
  X <- tibble(hhgde = X_hhgde, htele = X_htele, htslo = X_htslo, htrug = X_htrug, hspc1 = X_hspc1, hspc2 = X_hspc2, hcltp = X_hcltp, hcltt = X_hcltt)

  # dynamic predictors

  # time
  Y_htime <- as.vector(PredData$CovTemp[,"htime",])

  # season
  Y_hseas <- as.vector(PredData$CovTemp[,"hseas",])

  # persistent green
  Y_hhpgr <- as.vector(PredData$CovTemp[,"hhpgr",]) %>% scale(center = NimbleData$ScalingY$hhpgr[1], scale = NimbleData$ScalingY$hhpgr[2])
  Y_hhpgr <- Y_hhpgr %>% as.vector()

  # persistent green buffer (2km)
  Y_hhpgr2km <- as.vector(PredData$CovTemp[,"hhpgr2km",]) %>% scale(center = NimbleData$ScalingY$hhpgr2km[1], scale = NimbleData$ScalingY$hhpgr2km[2])
  Y_hhpgr2km <- Y_hhpgr2km %>% as.vector()

  # koala habitat
  Y_hhkha <- as.vector(PredData$CovTemp[,"hhkha",]) %>% as.factor()

  # woody cover - using proportion of each cell that is woody + sparse woody
  Y_hhfwc <- as.vector(PredData$CovTemp[,"hhfwc_1",] + PredData$CovTemp[,"hhfwc_2",])
  Y_hhfwc <- Y_hhfwc %>% scale(center = NimbleData$ScalingY$hhfwc[1], scale = NimbleData$ScalingY$hhfwc[2])
  Y_hhfwc <- Y_hhfwc %>% as.vector()

  # precipitation - centre on mean long term precipitation
  Y_hcpre <- as.vector(PredData$CovTemp[,"hcpre",]) - (rep(as.vector(PredData$CovCons[,"hcltp"]), length(as.vector(PredData$CovTemp[,"hcpre",]))/length(as.vector(PredData$CovCons[,"hcltp"]))) / 2)
  # remove season effects
  Y_hcpre[which(as.vector(PredData$CovTemp[,"hseas",]) == 1)] <- Y_hcpre[which(as.vector(PredData$CovTemp[,"hseas",]) == 1)] - mean(Y_hcpre[which(as.vector(PredData$CovTemp[,"hseas",]) == 1)])
  Y_hcpre[which(as.vector(PredData$CovTemp[,"hseas",]) == 0)] <- Y_hcpre[which(as.vector(PredData$CovTemp[,"hseas",]) == 0)] - mean(Y_hcpre[which(as.vector(PredData$CovTemp[,"hseas",]) == 0)])
  # centre and scale
  Y_hcpre <- Y_hcpre %>% scale(center = NimbleData$ScalingY$hcpre[1], scale = NimbleData$ScalingY$hcpre[2])
  Y_hcpre <- Y_hcpre %>% as.vector()

  # mean temperature - centre on long-term mean mean temperature
  Y_hctmn <- as.vector(PredData$CovTemp[,"hctmn",]) - (rep(as.vector(PredData$CovCons[,"hcltt"]), length(as.vector(PredData$CovTemp[,"hctmn",]))/length(as.vector(PredData$CovCons[,"hcltt"]))))
  # remove season effects
  Y_hctmn[which(as.vector(PredData$CovTemp[,"hseas",]) == 1)] <- Y_hctmn[which(as.vector(PredData$CovTemp[,"hseas",]) == 1)] - mean(Y_hctmn[which(as.vector(PredData$CovTemp[,"hseas",]) == 1)])
  Y_hctmn[which(as.vector(PredData$CovTemp[,"hseas",]) == 0)] <- Y_hctmn[which(as.vector(PredData$CovTemp[,"hseas",]) == 0)] - mean(Y_hctmn[which(as.vector(PredData$CovTemp[,"hseas",]) == 0)])
  # centre and scale
  Y_hctmn <- Y_hctmn %>% scale(center = NimbleData$ScalingY$hctmn[1], scale = NimbleData$ScalingY$hctmn[2])
  Y_hctmn <- Y_hctmn %>% as.vector()

  # max temperature
  Y_hctma <- as.vector(PredData$CovTemp[,"hctma",]) - (rep(as.vector(PredData$CovCons[,"hcltt"]), length(as.vector(PredData$CovTemp[,"hctma",]))/length(as.vector(PredData$CovCons[,"hcltt"]))))
  # remove season effects
  Y_hctma[which(as.vector(PredData$CovTemp[,"hseas",]) == 1)] <- Y_hctma[which(as.vector(PredData$CovTemp[,"hseas",]) == 1)] - mean(Y_hctma[which(as.vector(PredData$CovTemp[,"hseas",]) == 1)])
  Y_hctma[which(as.vector(PredData$CovTemp[,"hseas",]) == 0)] <- Y_hctma[which(as.vector(PredData$CovTemp[,"hseas",]) == 0)] - mean(Y_hctma[which(as.vector(PredData$CovTemp[,"hseas",]) == 0)])
  # centre and scale
  Y_hctma <- Y_hctma %>% scale(center = NimbleData$ScalingY$hctma[1], scale = NimbleData$ScalingY$hctma[2])
  Y_hctma <- Y_hctma %>% as.vector()

  # land-use
  # reclassfy so that: 1 = natural, 2 = production from natural, 3 = intensive, 4 = other
  Y_htlus <- as.vector(PredData$CovTemp[,"htlus",]) %>% as.character() %>% case_match("1"~1,"2"~2,"3"~4,"4"~4, "5"~3, "6"~4, "7"~4) %>% as.numeric() %>% as.factor()

  # intensive land-use buffer (2km)
  Y_htilu2km <- as.vector(PredData$CovTemp[,"htilu2km",]) %>% scale(center = NimbleData$ScalingY$htilu2km[1], scale = NimbleData$ScalingY$htilu2km[2])
  Y_htilu2km <- Y_htilu2km %>% as.vector()

  # lot size buffer (2km)
  Y_htpls2km <- as.vector(PredData$CovTemp[,"htpls2km",]) %>% scale(center = NimbleData$ScalingY$htpls2km[1], scale = NimbleData$ScalingY$htpls2km[2])
  Y_htpls2km <- Y_htpls2km %>% as.vector()

  # variables to mask out areas that are not habitat

  # woody cover - using proportion of each cell that is woody + sparse woody
  # used to adjust abundance estimates based on proportion of grids cells that are woody
  Y_propwoody <- as.vector(PredData$CovTemp[,"hhfwc_1",] + PredData$CovTemp[,"hhfwc_2",])

  # koala habitat
  # used to mask out grid cells that are rainforest and non-habitat so that: 0 = rainforest and non-habitat, 1 = everything else
  Y_khabitat <- as.vector(PredData$CovTemp[,"hhkha",]) %>% as.character() %>% case_match("1"~1,"2"~1,"3"~1,"4"~0, "5"~1, "6"~1, "7"~1, "8"~0) %>% as.numeric()

  # get the proportion habitat by multiplying the above variables together
  if (RainForestMask) {
      Y_prophabitat <- Y_propwoody * Y_khabitat
  } else {
      Y_prophabitat <- Y_propwoody
  }

  # compile into a tibble
  Y_temp <- tibble(htime = Y_htime, hseas = Y_hseas, hhpgr = Y_hhpgr, hhpgr2km = Y_hhpgr2km, hhkha = Y_hhkha, hhfwc = Y_hhfwc, hcpre = Y_hcpre, hctmn = Y_hctmn, hctma = Y_hctma, htlus = Y_htlus, htilu2km = Y_htilu2km, htpls2km = Y_htpls2km, prophabitat = Y_prophabitat)

  # get the design matrix and remove collinear variables for X
  X <- model.matrix(as.formula(paste0("~ ", paste(NimbleData$StaticVars, collapse = " + "))), model.frame(as.formula(paste0("~ ", paste(NimbleData$StaticVars, collapse = " + "))), as.data.frame(X), na.action = "na.pass")) %>% as.data.frame()

  # remove the intercept term
  X <- X[, 2:ncol(X)] %>% as.data.frame()

  # get the design matrix and remove collinear variables for X
  Y_temp <- model.matrix(as.formula(paste0("~ ", paste(c(NimbleData$DynamicVars, "prophabitat"), collapse = " + "))), model.frame(as.formula(paste0("~ ", paste(c(NimbleData$DynamicVars, "prophabitat"), collapse = " + "))), as.data.frame(Y_temp), na.action = "na.pass"))

  # remove the intercept term
  Y_temp <- Y_temp[, 2:ncol(Y_temp)] %>% as.data.frame()

  # get number of variables in Y
  NY <- ncol(Y_temp)

  # create 3D matrix of covariates for Y
  Y <- array(NA, dim = c(dim(PredData$CovTemp)[1], dim(PredData$CovTemp)[3], NY))
  for (i in 1:NY) {
      Y[,,i] <- matrix(Y_temp[,i], nrow = dim(PredData$CovTemp)[1], ncol = dim(PredData$CovTemp)[3])
  }

  # reset First and Last Date IDs
  FirstDateID <- FirstDateID - (NimbleData$FirstDateID_Orig - NimbleData$Constants$FirstDateID)
  LastDateID <- LastDateID - (NimbleData$FirstDateID_Orig - NimbleData$Constants$FirstDateID)

  return(list(Year = Year, SGridID =  SGridID, GenPopID = GenPopID, X = X, Y = Y, FirstDateID_Orig = NimbleData$FirstDateID_Orig, LastDateID_Orig = NimbleData$LastDateID_Orig, FirstDateID = FirstDateID, LastDateID = LastDateID, Order = NimbleData$Constants$Order, Lag = NimbleData$Constants$Lag, VarTrend = NimbleData$Constants$VarTrend, StaticVars = NimbleData$StaticVars, DynamicVars = NimbleData$DynamicVars, NamesX = NimbleData$NamesX, NamesY = NimbleData$NamesY))
}

# function to get predictions from a model
# Data and output generated from get_prediction_data()
get_predictions <- function(MCMC, Data) {

  # get the spatio-temporal random effects if needed
  if (Data$VarTrend == 1) {
    STre <- select(as_tibble(MCMC),contains("std["))
  } else {
    # get temporal random effects
    Tre <- select(as_tibble(MCMC),contains("td["))
  }

  # get regression coefficients
  Betas <- select(as_tibble(MCMC),contains("beta_d["))

  # remove grids with no genetic population ID
  Data$X <- Data$X %>% filter(!is.na(Data$GenPopID))
  Data$Y <- Data$Y[which(!is.na(Data$GenPopID)), , ]
  Data$SGridID <- Data$SGridID[which(!is.na(Data$GenPopID))]
  Data$GenPopID <- Data$GenPopID[which(!is.na(Data$GenPopID))]

  # get fixed linear predictors
  LPFixed <- apply(Betas[,1:dim(Data$X)[2]], MARGIN = 1, function(y) {as.matrix(Data$X) %*% as.matrix(y)})

  # loop through time steps
  for (i in Data$FirstDateID:Data$LastDateID) {

    # get time variable linear predictors
    if (Data$VarTrend == 1) {
      LPTime <- LPFixed + apply(Betas[,(dim(Data$X)[2] + 1):(dim(Data$X)[2] + dim(Data$Y)[3] - 1)], MARGIN = 1, function(z) {as.matrix(Data$Y[, (i - Data$FirstDateID + 1), 1:(dim(Data$Y)[3] - 1)]) %*% as.matrix(z)}) + t(apply(as.matrix(Data$GenPopID), MARGIN = 1, function(x) {t(STre[, (3 * (floor((i - Data$FirstDateID_Orig) / 2) + 1 - 1) + x)])}))
    } else {
      LPTime <- LPFixed + apply(Betas[,(dim(Data$X)[2] + 1):(dim(Data$X)[2] + dim(Data$Y)[3] - 1)], MARGIN = 1, function(z) {as.matrix(Data$Y[, (i - Data$FirstDateID + 1), 1:(dim(Data$Y)[3] - 1)]) %*% as.matrix(z)}) + t(apply(as.matrix(Data$GenPopID), MARGIN = 1, function(x) {t(Tre[, floor((i - Data$FirstDateID_Orig) / 2) + 1])}))
    }

    # multiply density estimates by the habitat availability to account for habitat masked out (last covariate in Y)
    if (i == Data$FirstDateID) {
      Density <- sweep(exp(LPTime), MARGIN = 1, as.matrix(Data$Y[ , (i - Data$FirstDateID + 1 + Data$Lag), dim(Data$Y)[3]]), `*`)
    } else {
      Density <- Density + sweep(exp(LPTime), MARGIN = 1, as.matrix(Data$Y[ , (i - Data$FirstDateID + 1 + Data$Lag), dim(Data$Y)[3]]), `*`)
    }

    Density <- Density / (Data$LastDateID - Data$FirstDateID + 1)
  }

  Mean <- apply(Density, MARGIN = 1, function(x) {mean(x, na.rm = TRUE)})
  Lower <- apply(Density, MARGIN = 1, function(x) {quantile(x, 0.025, na.rm = TRUE)})
  Upper <- apply(Density, MARGIN = 1, function(x) {quantile(x, 0.975, na.rm = TRUE)})
  SD <- apply(Density, MARGIN = 1, function(x) {sd(x, na.rm = TRUE)})
  Spatial <- tibble(GridID = Data$SGridID, Expected = Mean, LowerCI = Lower, UpperCI = Upper, SD = SD) %>% mutate(Expected = ifelse(is.na(Expected), NA, Expected), LowerCI = ifelse(is.na(LowerCI), NA, LowerCI), UpperCI = ifelse(is.na(UpperCI), NA, UpperCI), SD = ifelse(is.na(SD), NA, SD))
  # remove grids with missing data
  Spatial <- Spatial %>% filter(!is.na(Expected))
  TotalSumMean <- mean(apply(Density * 25, MARGIN = 2, function(x) {sum(x, na.rm = TRUE)}))
  TotalSumLower <- quantile(apply(Density * 25, MARGIN = 2, function(x) {sum(x, na.rm = TRUE)}), 0.025, na.rm = TRUE)
  TotalSumUpper <- quantile(apply(Density * 25, MARGIN = 2, function(x) {sum(x, na.rm = TRUE)}), 0.975, na.rm = TRUE)
  TotalSD <- sd(apply(Density * 25, MARGIN = 2, function(x) {sum(x, na.rm = TRUE)}))

  Output <- list(Spatial = Spatial, TotalSumMean = TotalSumMean, TotalSumLowerCI = TotalSumLower, TotalSumUpperCI = TotalSumUpper, TotalSD = TotalSD)

  return(Output)
}


fcn_update_db <- function(db = "KoalaSurveyData2020_cur.accdb",
                          db.integrated = "Integrated_SEQKoalaDatabase.accdb",
                          survey_sites = "KoalaSurveySites_231108.shp",
                          update.spatial = T,
                          spatial_transects = "transects_v240517.shp",
                          line_transect_width = 28.69,
                          site.info = "site_info.csv",
                          monitoring.units = "Population_boundaries_v2.shp"){
  
  
  # load function to create line transects as line feature
  source("code/line-transect.r")
  
  cat(
    "------------------------------------------------------------","\n",
    "IT MAY TAKE AROUND 10 MINS TO RUN DEPENDING ON YOUR COMPUTER","\n",
    "------------------------------------------------------------", "\n\n")
  
  # look for the database (db) within the working directory
  path <- list.files(pattern = db, recursive = T, ignore.case = T, full.names = F, include.dirs = T)
  if (length(path) == 0) {
    stop(sprintf("No % s file found. Please, make sure you have a file named % s in this working directory",
                 db, db))} # close {} of error message
  if(length(path) > 1){
    message("Multiple files found:")
    for (i in seq_along(path)) {
      cat(i, ":", path[i], "\n")}
    choice <- as.integer(readline("Enter the number corresponding to the file you want to use: "))
    if (is.na(choice) || choice < 1 || choice > length(path)) {
      stop("Invalid selection. Process terminated. Please, rerun this function and choose one of the provided options.")}
    path <- path[choice]}
  
  # create the db path for the odbc connection
  db_path <- paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};",
                    "DBQ=",
                    getwd(), "/", path)
  
  # establish odbc connection
  channel <- RODBC::odbcDriverConnect(db_path)
  
  # retrieve survey data table only
  dat <- RODBC::sqlFetch(channel, "tblKoalaSurveyData2020_cur")
  
  # close odbc connection
  RODBC::odbcClose(channel)
  
  # include number of sightings
  
  
  # adding the 2018 Daisy Hill Line Transects
  # create a temporary dataframe for merging
  tmp <- dat |> 
    dplyr::mutate(XID = as.character(ID),
                  Start_Time = format(Start_Time, "%H:%M"),
                  End_Time = format(End_Time, "%H:%M"),
                  Time = format(Time, "%H:%M"),
                  Year = as.numeric(format(Date, "%Y")))
  
  # standardize column names
  SOL_2018 <- read.csv("input/line_transects_2018/V2_MASTER_Daisy Hill_Line_tidy.csv") |>  
    dplyr::mutate(Date = as.POSIXct(Date),
                  Weather = as.character(Weather),
                  DBH = as.character(DBH),
                  Wind = as.character(Wind),
                  Cloud_Cover = as.character(Cloud_Cover),
                  Canopy_Cover = as.character(Canopy_Cover)) |> 
    dplyr::select(-Concern, -Fatal_problem)
  
  # standardize the 2018 line transects dataframe 
  dat2 <- dplyr::bind_rows(tmp, SOL_2018) |> 
    dplyr::relocate(Year, .before = Date) |> 
    dplyr::relocate(Number_Sightings, .before = Sighting_Number) |>
    dplyr::rename(NumSightings = Number_Sightings) |> 
    dplyr::mutate(ID = dplyr::row_number())
  
  # include number of sightings per transect
  trans1 <- dat2 |>
    dplyr::mutate(Sighting_Number = ifelse(is.na(Sighting_Number),
                                           0, Sighting_Number)) |> 
    dplyr::group_by(Transect_ID) |>
    dplyr::mutate(NumSightings =
                    if (dplyr::n() > 1) {
                      NumSightings = dplyr::n()
                    }
                  else if (dplyr::n() == 1 & Sighting_Number != 0) {
                    NumSightings = 1
                  }
                  else {
                    NumSightings = 0
                  }) |> 
    dplyr::ungroup() |> 
    dplyr::mutate(Number_Sightings = ifelse(is.na(NumSightings),
                                            NumSightings,
                                            NumSightings)) |> 
    dplyr::relocate(Number_Sightings, .before = Sighting_Number) |> 
    dplyr::select(-NumSightings)
  
  
  # add length and area if not recorded
  trans.sf <- trans1 |> 
    dplyr::filter(!Method %in% c("IS", "VT")) |> 
    dplyr::group_by(Transect_ID) |> 
    dplyr::reframe(
      Start_Eastings = unique(Start_Eastings), 
      Start_Northings = unique(Start_Northings),
      End_Eastings = unique(End_Eastings),
      End_Northings = unique(End_Northings),
      Year = unique(Year)
    ) |> 
    fcn_line_transect_sf() |> 
    dplyr::mutate(TLength2 = as.numeric(st_length(geometry))) |> 
    dplyr::mutate(index = dplyr::row_number()) |> 
    sf::st_buffer(line_transect_width/2, endCapStyle = "FLAT")
  
  trans2 <- trans1 |> 
    dplyr::left_join(sf::st_drop_geometry(trans.sf) |> dplyr::select(Transect_ID, TLength2),
                     by = "Transect_ID") |> 
    dplyr::mutate(T_Length = ifelse(is.na(T_Length),
                                    TLength2, T_Length)) |>
    dplyr::mutate(T_Area = ifelse(is.na(T_Area),
                                  T_Length * line_transect_width/1e4, T_Area),
                  T_Width = ifelse(is.na(T_Width),
                                   line_transect_width, T_Width)) |> 
    dplyr::select(-TLength2)
  
  
  #-----------------------------------------------------------
  # standardise Site ID
  # look for a file with the spatial representation of survey sites within the working directory
  path <- list.files(pattern = survey_sites, recursive = T, ignore.case = T, full.names = F, include.dirs = T)
  if (length(path) == 0) {
    stop(sprintf("No % s file found. Please, make sure you have a file named % s in this working directory",
                 survey_sites, survey_sites))} # close {} of error message
  if(length(path) > 1){
    message("Multiple files found:")
    for (i in seq_along(path)) {
      cat(i, ":", path[i], "\n")}
    choice <- as.integer(readline("Enter the number corresponding to the file you want to use: "))
    if (is.na(choice) || choice < 1 || choice > length(path)) {
      stop("Invalid selection. Process terminated. Please, rerun this function and choose one of the provided options.")}
    path <- path[choice]}
  
  surveys <- sf::st_read(path) 
  sf::st_crs(surveys) <- "EPSG:7856"
  
  siteIds <-
    suppressWarnings({
      lapply(split(trans.sf, trans.sf$Transect_ID), function(d) {
        cat(
          paste0(
            "Matching Site_ID (",
            round(d$index/max(trans.sf$index)*100, 0),
            "% processed)","\n", "Transect_ID: ", d$Transect_ID, "\n\n"))
        
        sub <- sf::st_intersection(surveys, d) |>
          dplyr::select(Survey_Yr, Site, NAME, Subsite)
        
        if (nrow(sub) == 0) {
          Site_ID <- NA
          Subsite <- NA
          Site_Name <- NA
          
        } else if (nrow(sub) == 1) {
          Site_ID <- unique(sub$Site)
          Subsite <- unique(sub$Subsite)
          Site_Name <- unique(sub$NAME)
          
        } else {
          if (length(unique(sub$Site)) > 1) {
            sub <- sub |> 
              dplyr::mutate(inters.area = as.numeric(sf::st_area(geometry)) /
                              1e4) |>
              dplyr::select(Survey_Yr, Site, NAME, Subsite, inters.area)
            
            max <- max(sub$inters.area, na.rm = T)
            sub <- sub[sub$inters.area == max,]
            
            year <-
              unique(unique(d$Year))
            survey_yr <- sub$Survey_Yr
            
            if (any(year %in% survey_yr)) {
              filter <- survey_yr == year
            } else {
              filter <-
                survey_yr == survey_yr[which.max(survey_yr[which(survey_yr <= year)])]
            }
          } else {
            year <-
              unique(unique(d$Year))
            survey_yr <- sort(unique(sub$Survey_Yr))
            
            if (any(year %in% survey_yr)) {
              filter <- survey_yr == year
            } else {
              filter <-
                survey_yr == survey_yr[which.max(survey_yr[which(survey_yr <= year)])]
            }
          }
          
          
          if (purrr::is_empty(filter)) {
            diff <- abs(year - survey_yr)
            filter <-
              survey_yr == survey_yr[which.min(diff)]
          }
          
          if (length(filter) != nrow(sub)) {
            sub <- sf::st_intersection(surveys, d) |>
              dplyr::mutate(inters.area = as.numeric(sf::st_area(geometry)) /
                              1e4) |>
              dplyr::select(Survey_Yr, Site, NAME, Subsite, inters.area)
            sub <- sub[-which.min(sub$inters.area),]
            
            survey_yr <- sort(unique(sub$Survey_Yr))
            
            if (any(year %in% survey_yr)) {
              filter <- survey_yr == year
            } else {
              filter <-
                survey_yr == survey_yr[which.max(survey_yr[which(survey_yr <= year)])]
            }
          }
          
          sub <- sub |>
            dplyr::distinct(Survey_Yr, Site, .keep_all = T)
          
          Site_ID <- sub[filter, "Site"]$Site
          Subsite <- sub[filter, "Subsite"]$Subsite
          Site_Name <- sub[filter, "NAME"]$NAME
        }
        
        df <- data.frame(
          Transect_ID = unique(d$Transect_ID),
          Site_ID = Site_ID,
          Subsite = Subsite,
          Site_Name = Site_Name
        )
        
        return(df)
      }) |> dplyr::bind_rows()
    })
  
  # Check if there are any duplicated Transect_IDs resulting from different Site_ID values for the same transect
  siteIds <- siteIds |>
    dplyr::mutate(dups = duplicated(Transect_ID))
  
  # Identify transects without a Site_ID
  siteNAs <- siteIds |>
    dplyr::filter(is.na(Site_ID)) |> 
    dplyr::pull(unique(Transect_ID))
  
  siteNAs <- dat |> 
    dplyr::filter(!Method %in% c("IS", "VT") & Transect_ID %in% siteNAs)
  
  siteNAs.n <- nrow(siteNAs)
  
  if (siteNAs.n > 0) {
    assign("siteNAs", siteNAs, envir = .GlobalEnv)
    
    warning(sprintf("There are %s transects without a matching site_ID in %s. Please ensure you are using the most up-to-date spatial representation of survey sites, as outdated information could make it hard to get accurate data for the model. To see which transects are missing a matching site_ID, check the siteNAs dataframe in the top right panel (i.e., your global environment)", siteNAs.n, survey_sites))
  } # close if (siteNAs.n > 0)
  
  # Assign the spatially retrieved Site_ID value to each transect of the data frame
  trans3 <- trans2 |> 
    dplyr::ungroup() |> 
    dplyr::left_join(siteIds, by = "Transect_ID") |> 
    
    # Resolve problems with duplicated column names
    dplyr::mutate(Site_Name = ifelse(is.na(Site_Name.y), Site_Name.x, Site_Name.y),
                  Site_ID = ifelse(is.na(Site_ID.y), Site_ID.x, Site_ID.y),
                  SiteID_Master = Site_ID.x,
                  Subsite = ifelse(is.na(Subsite.y), NA, Subsite.y),
                  Subsite_Master = Subsite.x) |> 
    dplyr::select(-Site_Name.x, -Site_Name.y, -Site_ID.x, 
                  -Site_ID.y, -Subsite.x, -Subsite.y) |> 
    dplyr::relocate(Site_ID, .after = LGA) |> 
    dplyr::relocate(Subsite, .after = Site_ID) |> 
    dplyr::relocate(Site_Name, .after = Subsite) |> 
    dplyr::relocate(Number_Observers, .before = Observer_Obstacles) |> 
    
    # Remove unwanted columns
    dplyr::select(-Number_Observers_Master, -dups)
  
  # Define the source of the 2018 transects"V2_MASTER_Daisy Hill_Line"
  trans3 <- trans3 |> 
    dplyr::mutate(Source = ifelse(Year == 2018, "V2_MASTER_Daisy Hill_Line", "KoalaSurveyData2020_cur"))
  
  # Manually assign Site_ID and Subsite of a problematic transect
  trans3[trans3$Transect_ID == "221.0_10_SOL.20220503",
         "Site_ID"] <- 221
  trans3[trans3$Transect_ID == "221.0_10_SOL.20220503",
         "Subsite"] <- 3
  
  
  #-----------------------------------------------------------
  # UPDATE SHAPEFILE
  # look for the most up-to-date spatial representation of transects within the working directory
  if(update.spatial){
    path <- list.files(pattern = spatial_transects, recursive = T, ignore.case = T, full.names = F, include.dirs = T)
    if (length(path) == 0) {
      stop(sprintf("No % s file found. Please, make sure you have a file named % s in this working directory",
                   spatial_transects, spatial_transects))} # close {} of error message
    if(length(path) > 1){
      message("Multiple files found:")
      for (i in seq_along(path)) {
        cat(i, ":", path[i], "\n")}
      choice <- as.integer(readline("Enter the number corresponding to the file you want to use: "))
      if (is.na(choice) || choice < 1 || choice > length(path)) {
        stop("Invalid selection. Process terminated. Please, rerun this function and choose one of the provided options.")}
      path <- path[choice]}
    
    transects.sf <- sf::st_read(path) |> 
      dplyr::filter(!Source %in% gsub(".accdb", "", db)) 
    sf::st_crs(transects.sf) <- "EPSG:7856"
    
    lineTransect <- trans3 |> 
      
      # Remove incidental sightings and vehicle transects
      dplyr::filter(Method %in% c("DOL", "SOL") & Year > 2018) |> 
      
      # Reframe to retain only columns relevant for modelling
      dplyr::group_by(Transect_ID) |> 
      dplyr::reframe(
        Date = as.character(unique(Date)),
        Start_Eastings = unique(Start_Eastings), 
        Start_Northings = unique(Start_Northings),
        End_Eastings = unique(End_Eastings),
        End_Northings = unique(End_Northings),
        Method = "SOL", 
        Source = gsub(".accdb", "", db)
      ) 
    # Create a simple feature object to be exported as a shapefile
    suppressWarnings({
      lineTransectSf <- fcn_line_transect_sf(lineTransect) |> 
        sf::st_buffer(line_transect_width, endCapStyle = "FLAT") |>
        dplyr::mutate(Date = format(as.POSIXct(Date,"%Y-%m-%d"), "%d/%m/%Y")) |> 
        dplyr::select(TrnscID = Transect_ID, Date, Method, Source)
    })
    
    # Append to the source spatial representation
    trans.final <- dplyr::bind_rows(transects.sf, lineTransectSf) |> 
      dplyr::mutate(Year = as.integer(substring(Date, nchar(Date) - 4 + 1)))
    
    # Export as shapefile
    if(!file.exists("output/transects_spatial_representation")){
      dir.create("output/transects_spatial_representation")}
    suppressWarnings({
      sf::st_write(trans.final, 
                   paste0("output/transects_spatial_representation/transects_v",
                          format(Sys.Date(), "%y%m%d"),
                          ".shp"), append = F)
    })
    
    files.from <- list.files(pattern = paste0("transects_v", format(Sys.Date(), "%y%m%d")),
                             path = "output/transects_spatial_representation",
                             full.names = T)
    files.to <- sub("/([^/]*)\\.", "/Integrated_SEQKoalaDatabase_Spatial_updated.", files.from)
    files.to <- sub("output", "input", files.to)
    
    file.copy(from = files.from, 
              to = files.to,
              overwrite = T)
    
  } # close conditional to update spatial representation of transects
  
  # Fill any missing value in number of observers (i.e., SOL = 1 and DOL = 2)
  trans3[is.na(trans3$Number_Observers) & trans3$Method == "DOL", "Number_Observers"] <- 2
  trans3[is.na(trans3$Number_Observers) & trans3$Method == "SOL", "Number_Observers"] <- 1
  
  #-----------------------------------------------------------
  # Survey site information and monitoring unit
  path <- list.files(pattern = site.info, recursive = T, ignore.case = T, full.names = F, include.dirs = T)
  if (length(path) == 0) {
    stop(sprintf("No % s file found. Please, make sure you have a file named % s in this working directory",
                 site.info, site.info))} # close {} of error message
  if(length(path) > 1){
    message("Multiple files found:")
    for (i in seq_along(path)) {
      cat(i, ":", path[i], "\n")}
    choice <- as.integer(readline("Enter the number corresponding to the file you want to use: "))
    if (is.na(choice) || choice < 1 || choice > length(path)) {
      stop("Invalid selection. Process terminated. Please, rerun this function and choose one of the provided options.")}
    path <- path[choice]}
  
  # Load file with habitat information
  info <- read.csv(path)
  
  trans4 <- trans3 |>
    dplyr::left_join(
      info |> dplyr::select(LGA, KPA, Site.Number,
                            Land.Type, Habitat.Type, Site.Name),
      dplyr::join_by(Site_ID == Site.Number)) |>
    dplyr:: mutate(
      LGA = ifelse(is.na(LGA.x), LGA.y, LGA.x),
      KPA = ifelse(is.na(KPA.x),LGA.y, KPA.x),
      Hab_Type = ifelse(is.na(Hab_Type), Habitat.Type, Hab_Type),
      Land_Type = ifelse(is.na(Land_Type), Land.Type, Land_Type),
      Site_Name_Master = ifelse(is.na(Site_Name_Master), Site.Name, Site_Name_Master),
      Site_Name = ifelse(is.na(Site_Name), Site.Name, Site_Name)) |>
    dplyr::select(-LGA.y,-LGA.x,-KPA.x, -KPA.y, -Land.Type, -Habitat.Type, -Site.Name)
  
  # Assign monitoring unit from shapefile
  # This code follows the same logic applied to extract Site ID
  path <- list.files(pattern = monitoring.units, recursive = T, ignore.case = T, full.names = F, include.dirs = T) 
  path <- path[(substring(path, nchar(path) - 3 + 1) == "shp")]
  
  if (length(path) == 0) {
    stop(sprintf("No % s file found. Please, make sure you have a file named % s in this working directory",
                 monitoring.units, monitoring.units))} # close {} of error message
  if(length(path) > 1){
    message("Multiple files found:")
    for (i in seq_along(path)) {
      cat(i, ":", path[i], "\n")}
    choice <- as.integer(readline("Enter the number corresponding to the file you want to use: "))
    if (is.na(choice) || choice < 1 || choice > length(path)) {
      stop("Invalid selection. Process terminated. Please, rerun this function and choose one of the provided options.")}
    path <- path[choice]}
  
  units <- sf::st_read(path)
  sf::st_crs(units) <- 7856
  
  tmp <- trans4 |> 
    # Remove incidental sightings and vehicle transects
    dplyr::filter(Method %in% c("DOL", "SOL")) |> 
    # Reframe to retain only columns relevant for modelling
    dplyr::group_by(Transect_ID) |> 
    dplyr::reframe(
      Date = as.character(unique(Date)),
      Start_Eastings = unique(Start_Eastings), 
      Start_Northings = unique(Start_Northings),
      End_Eastings = unique(End_Eastings),
      End_Northings = unique(End_Northings),
      Method = "SOL", 
      Source = gsub(".accdb", "", db)
    ) |> 
    fcn_line_transect_sf() |>  
    dplyr::mutate(index = dplyr::row_number()) |> 
    sf::st_buffer(line_transect_width/2, endCapStyle = "FLAT")
  
  M.unit <- suppressWarnings({
    lapply(split(tmp, tmp$Transect_ID), function(d){
      cat(
        paste0(
          "Matching monitoring unit (",
          round(d$index/max(tmp$index)*100, 0),
          "% processed)","\n", "Transect_ID: ", d$Transect_ID, "\n\n"))
      
      sub <- sf::st_intersection(units, d) |> 
        dplyr::select(GENPOP_NAM, Transect_ID)
      
      if(nrow(sub) == 0) {
        Monitoring_Unit <- NA
        Transect_ID <- d$Transect_ID
        
      } else if (nrow(sub) == 1) {
        Monitoring_Unit <- unique(sub$GENPOP_NAM)
        Transect_ID <- unique(sub$Transect_ID)
        
      } else if(length(unique(sub$GENPOP_NAM)) > 1){
        sub <- sf::st_intersection(units, d) |>
          dplyr::mutate(inters.area = as.numeric(st_area(geometry)) / 1e4) |>
          dplyr::select(GENPOP_NAM, Transect_ID, inters.area) 
        
        max <- max(sub$inters.area, na.rm = T)
        sub <- sub[sub$inters.area == max, ]
        
        Monitoring_Unit <- sub$GENPOP_NAM
        Transect_ID <- unique(sub$Transect_ID)
      }
      
      df <- data.frame(Transect_ID = Transect_ID,
                       Monitoring_Unit = Monitoring_Unit)
      
      return(df)
    })
  }) |> dplyr::bind_rows()
  
  # Check any missing data
  M.unit.NAs <- M.unit |> 
    dplyr::filter(is.na(Monitoring_Unit))
  if (nrow(M.unit.NAs) > 0) {
    assign("M.unit.NAs", M.unit.NAs, envir = .GlobalEnv)
    stop(sprintf("There are %s transects without a matching monitoring unit. Please, check transects in the dataframe named M.unit.NAs before continuing.", nrow(M.unit.NAs)))} # close {} of error message
  
  # Assign a corresponding monitoring unit to the dataframe
  trans5 <- trans4 |> 
    dplyr::left_join(M.unit, by = "Transect_ID") |> 
    dplyr::rename(Monitoring_Unit = Monitoring_Unit.y) |> 
    dplyr::select(-Monitoring_Unit.x)
  
  #-----------------------------------------------------------
  LTPerpDist <- trans5 |>
    dplyr::filter(Method %in% c("SOL", "DOL") & Number_Sightings > 0) |> 
    dplyr::mutate(Perp_Dist = ifelse(is.na(Perp_Dist), 
                                     (Distance_Observer_Koala * 
                                        sin(Angle_Koala * pi/180) / 1), # convert from degrees to radians
                                     Perp_Dist)) |> 
    dplyr::relocate(Sighting_ID, .before = Transect_ID)
  
  # Check for missing data
  LTPerpDist.na <- LTPerpDist |> 
    dplyr::filter(is.na(Perp_Dist)) |> 
    
    # Remove one potentially duplicate entry (ID = 1816)
    dplyr::filter(!is.na(Koala_Eastings)) |> 
    sf::st_as_sf(coords = c("Koala_Eastings", "Koala_Northings"),
                 crs = 7856) # In this case is impossible to calculate the perpendicular distance because the Angle_Koala is also missing data
  
  
  # Calculate perpendicular distances using spatial representations of the transects and of koala sightings
  # Create an spatial object with transects
  LTPerpDist.sf <- LTPerpDist |> 
    dplyr::group_by(Transect_ID) |> 
    dplyr::reframe(
      Site_ID = unique(Transect_ID), 
      Date = unique(Date),
      T_Length = unique(T_Length),
      Number_Sightings = unique(Number_Sightings),
      Number_Observers = unique(Number_Observers),
      Start_Eastings = unique(Start_Eastings), 
      Start_Northings = unique(Start_Northings),
      End_Eastings = unique(End_Eastings),
      End_Northings = unique(End_Northings),
      Year = as.integer(unique(Year)),
      Sighting_ID = unique(Sighting_ID)
    ) |> 
    fcn_line_transect_sf() |> 
    dplyr::filter(Sighting_ID %in% LTPerpDist.na$Sighting_ID)
  
  
  # Calculate distance between the koala sighting and the transect line based on the spatial representation
  dist <- lapply(split(LTPerpDist.na, LTPerpDist.na$Sighting_ID), function(x){
    
    dist <- x |>
      sf::st_distance(
        dplyr::filter(LTPerpDist.sf,
                      LTPerpDist.sf$Sighting_ID == x$Sighting_ID)
      ) |> 
      as.numeric()
    
    z <- data.frame(Sighting_ID = x$Sighting_ID,
                    Perp_Dist = dist)
    
    return(z)
  }) |> 
    dplyr::bind_rows()
  
  
  # Update the data frame to fill the gap
  trans6 <- trans5 |> 
    dplyr::left_join(dist, by = "Sighting_ID") |> 
    dplyr::mutate(Perp_Dist = ifelse(is.na(Perp_Dist.x),
                                     Perp_Dist.y,
                                     Perp_Dist.x)) |> 
    dplyr::select(-Perp_Dist.x, -Perp_Dist.y)
  
  if(sum(duplicated(LTPerpDist$Sighting_ID)) > 0) {
    warning(sprintf("%i duplicate Sighting_ID", sum(duplicated(lineTransect$Sighting_ID))))
  } else {cat("There is a unique Sighting_ID per koala sighting")}
  
  # Final formatting
  # Create a template
  templ <- trans1[1, ] |> 
    dplyr::mutate(across(everything(), ~ NA)) |> 
    # Rename original Gen_Pop as Monitoring_Unit_Master
    dplyr::rename(Monitoring_Unit_Master = Gen_Pop) |> 
    dplyr::select(-ID, -XID) |> 
    # Change column type for combining
    dplyr::mutate(Date = as.character(Date))
  
  
  # Keep incidental sightings and vehicle transects separate
  IS_VT <- trans6 |> 
    dplyr::filter(!Method %in% c("SOL", "DOL")) |> 
    dplyr::mutate(Date = as.POSIXct(Date, format = "%d/%m/%Y"),
                  Tran_Plot = ifelse(is.na(Tran_Plot), Tran_plot, Tran_Plot)) |> 
    dplyr::rename("\"Date\"" = Date, "\"Time\"" = Time) 
  
  
  trans7 <- trans6 |>
    # Tidy up date format
    dplyr::mutate(Date = format(Date, format = "%d/%m/%Y"))  |> 
    # Assign empty values ("") instead of missing data (NAs)
    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        ~ ifelse(is.na(.), "", .))) |> 
    # Rename original Gen_Pop as Monitoring_Unit_Master
    dplyr::rename(Monitoring_Unit_Master = Gen_Pop) |> 
    # Tidying 
    dplyr::relocate(Perp_Dist, .after = Angle_Koala) |> 
    dplyr::select(-Site_Name_Master, -ID, -Distance_bt_Observers) |> 
    dplyr::mutate(Site_Area = as.numeric(Site_Area),
                  Subsite = as.numeric(Subsite))
  
  suppressWarnings({
    final <- dplyr::bind_rows(templ, trans7) |> 
      dplyr::mutate(Date = as.POSIXct(Date, format = "%d/%m/%Y"),
                    Tran_Plot = ifelse(is.na(Tran_Plot), Tran_plot, Tran_Plot),
                    Site_ID = as.character(Site_ID),
                    Total_Field_Staff = as.integer(Total_Field_Staff),
                    Start_Eastings = as.integer(Start_Eastings),
                    End_Eastings = as.integer(End_Eastings),
                    Start_Northings = as.integer(Start_Northings),
                    End_Northings = as.integer(End_Northings),
                    Number_Observers = as.integer(Number_Observers),
                    Count = as.integer(Count),
                    Observer_Eastings = as.integer(Observer_Eastings),
                    Observer_Northings = as.integer(Observer_Northings),
                    Koala_Eastings = as.integer(Koala_Eastings),
                    Koala_Northings = as.integer(Koala_Northings),
                    Koala_Elevation = as.integer(Koala_Elevation),
                    Koala_Bearing = as.integer(Koala_Bearing),
                    Angle_Koala = as.integer(Angle_Koala),
                    Field_Bearing = as.integer(Field_Bearing),
                    T_Length = as.double(T_Length),
                    T_Area = as.double(T_Area),
                    Perp_Dist = as.double(Perp_Dist),
                    Tree_Height = as.double(Tree_Height),
                    Koala_Height = as.double(Koala_Height)) |> 
      dplyr::rename("\"Date\"" = Date, "\"Time\"" = Time) |> 
      dplyr::relocate(Monitoring_Unit_Master, .after = Study_Area) |> 
      dplyr::relocate(Monitoring_Unit, .before = LGA) |> 
      dplyr::select(-Tran_plot, -XID)
    final <- final[-1, ]
  })
  
  # update db
  path <- list.files(pattern = db.integrated, recursive = T, ignore.case = T, full.names = F, include.dirs = T)
  if (length(path) > 1) {
    stop(sprintf("There are multiple files named %s. Please, make sure to keep only the most up-to-date koala survey database in this working directory. Databases found in: %s", db.integrated))} # close {} of error message
  
  # create a copy the database to avoid modifying the original
  if(!file.exists("output/integrated_database")){
    dir.create("output/integrated_database")}
  
  file.copy(from = path, 
            to = "output/integrated_database/Integrated_SEQKoalaDatabase_updated.accdb")
  
  # update path
  path.upd <- "output/integrated_database/Integrated_SEQKoalaDatabase_updated.accdb"
  
  # create the db path for the odbc connection
  db_path <- paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};",
                    "DBQ=",
                    getwd(), "/", path.upd)
  
  # establish odbc connection
  channel <- RODBC::odbcDriverConnect(db_path)
  
  # remove old SOL table survey data table only
  RODBC::sqlDrop(channel, "SOL_Compiled")
  
  RODBC::sqlSave(
    channel,
    dat = final,
    tablename = "SOL_Compiled",
    fast = F, 
    safer = F, 
    rownames = F, 
    colnames = F,
    varTypes = c("\"Date\"" = "date")
  )
  
  # close odbc connection
  RODBC::odbcClose(channel)
  
  file.copy(from = "output/integrated_database/Integrated_SEQKoalaDatabase_updated.accdb", 
            to = "input/databases/Integrated_SEQKoalaDatabase_updated.accdb", overwrite = T)
  
} # finish function


### Survey data from MS Access to R ###
fcn_all_tables_detect <- function(db_name = "integrated",
                                  script_name = c("line_transect_detect",
                                                  "strip_transect_detect", 
                                                  "uaoa_detect",
                                                  "perp_distance_detect"),
                                  observers_file = "group_observers_lookup.csv",
                                  columns = c("Weather",
                                              "Wind",
                                              "Cloud_Cover",
                                              "Canopy_Cover",
                                              "Subcanopy_Cover"),
                                  updated.db = T)
{
  # custom function to extract the maximum value when multiple values are separated by "_". For instance: in "2_3" only 3 is retained  
  split_sep <- function(column){
    ifelse(grepl("_", column),
           max(as.numeric(unlist(strsplit(column, "_"))), na.rm = T),
           column)}
  
  # extract and query tables from MS Access
  conn <- fcn_db_connect_table(db_name)
  table <- list()
  for(i in 1:length(script_name)){
    path <- list.files(
      pattern = script_name[i],
      recursive = T,
      ignore.case = T,
      full.names = F,
      include.dirs = T)
    
    query <- readr::read_file(path)
    table[[i]] <- RODBC::sqlQuery(conn, query, errors = TRUE)
    names(table)[i] <- gsub("_detect", "", script_name[i])
  }
  
  # change column names of "perp_distance" to avoid issues with the "Date" column 
  colnames(table[["perp_distance"]]) <- c("TransectID", "SightingID", "Perp_Dist", "Date")
  
  # close ODBC connection
  if (is.character(table)) {
    fcn_db_disconnect(conn)
    stop(table)
  }
  else {
    fcn_db_disconnect(conn)
  }
  
  # read lookup for obervers
  max.year <- lapply(table[which(names(table) %in% c("line_transect", "strip_transect", "uaoa"))], function(df){
    df <- df |> 
      mutate(Date = as.Date(Date, format = "%d/%m/%Y"),
             Year = lubridate::year(Date)) |>
      pull(Year) |> 
      max()
  }) |> unlist()
  if(max(max.year, na.rm = T) > 2023){
    warning(sprintf("Please, update the % s file to include groups of observers from 2024 onward. Otherwise, they will be recorded as missing data.", observers_file))}
  
  path <- list.files(pattern = observers_file, recursive = T, ignore.case = T, full.names = F, include.dirs = T)
  if (length(path) == 0) {
    stop(sprintf("No % s file found. Please, make sure you have a file named % s in this working directory",
        observers_file, observers_file))} # close {} of error message
  if(length(path) > 1){
    message("Multiple files found:")
    for (i in seq_along(path)) {
      cat(i, ":", path[i], "\n")}
    choice <- as.integer(readline("Enter the number corresponding to the file you want to use: "))
    if (is.na(choice) || choice < 1 || choice > length(path)) {
      stop("Invalid selection. Process terminated. Please, rerun this function and choose one of the provided options.")}
    path <- path[choice]}
  
   # generate a random sequence for groups of observers to preserve privacy
    set.seed(123)
    obs <- data.table::fread(path)[, ":="(start_date = as.POSIXct(start_date, format = "%d/%m/%Y"),
                                          end_date = as.POSIXct(end_date, format = "%d/%m/%Y"),
                                          ObserverGroup = sample(nrow(.SD), replace = F))]
    data.table::setkey(obs, start_date, end_date)

  # Process tables to include variables for koala detectability
  tables <- lapply(table[which(names(table) %in% c("line_transect", "strip_transect", "uaoa"))], function(df){
    
    # check duplicates
    if(anyDuplicated(df$TransectID)) {
      stop(sprintf("%s duplicate values found in TransectID", sum(duplicated(df$TransectID))))
    }
    
    df2 <- df |> 
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ifelse(. %in% c("NA's", ""), NA, .)),
                    # retain maximum values only
                    dplyr::across(.cols = dplyr::all_of(columns), split_sep),
                    # transform from character to numeric
                    dplyr::across(.cols = dplyr::all_of(columns), as.numeric),
                    # transform from character to Date
                    Date = as.Date(Date, format = "%d/%m/%Y"),
                    # extract day of the year
                    DayYear = lubridate::yday(Date),
                    # extract month
                    Month = lubridate::month(Date),
                    # extract time of the day
                    TimeDay = lubridate::hour(
                      lubridate::round_date(
                        as.POSIXct(
                          paste(Date, Start_Time),
                          format = "%Y-%m-%d %H:%M"),
                        unit = "hour")),
      ) |> 
      dplyr::select(-Start_Time) |> 
      # add TimePeriodID
      SEQKoalaDataPipeline::fcn_add_date_interval() 
    
    # standardize values to DESI's metadata
    if("Cloud_Cover" %in% names(df2)){
      df2 <- df2 |> 
        dplyr::mutate(Cloud_Cover = dplyr::case_when(Cloud_Cover > 5 & Cloud_Cover <= 20 ~ 1,
                                                     Cloud_Cover > 20 & Cloud_Cover <= 40 ~ 2,
                                                     Cloud_Cover > 40 & Cloud_Cover <= 60 ~ 3,
                                                     Cloud_Cover > 60 & Cloud_Cover <= 80 ~ 4,
                                                     Cloud_Cover > 80 & Cloud_Cover <= 100 ~ 5,
                                                     .default = Cloud_Cover))
    }
    
    if("Canopy_Cover" %in% names(df2)){
      df2 <- df2 |> 
        dplyr::mutate(Canopy_Cover = dplyr::case_when(Canopy_Cover == 1.5 ~ 2,
                                                      Canopy_Cover == 2.5 ~ 3,
                                                      Canopy_Cover == 3.5 ~ 4,
                                                      .default = Canopy_Cover)) |> 
        dplyr::mutate(Canopy_Cover = dplyr::case_when(Canopy_Cover > 5 & Canopy_Cover <= 20 ~ 1,
                                                      Canopy_Cover > 20 & Canopy_Cover <= 40 ~ 2,
                                                      Canopy_Cover > 40 & Canopy_Cover <= 60 ~ 3,
                                                      Canopy_Cover > 60 & Canopy_Cover <= 80 ~ 4,
                                                      Canopy_Cover > 80 & Canopy_Cover <= 100 ~ 5,
                                                      .default = Canopy_Cover))
    }
    
    if("Subcanopy_Cover" %in% names(df2)){
      df2 <- df2 |> 
        dplyr::mutate(Subcanopy_Cover = dplyr::case_when(Subcanopy_Cover == 1.5 ~ 2,
                                                         Subcanopy_Cover == 2.5 ~ 3,
                                                         Subcanopy_Cover == 3.5 ~ 4,
                                                         Subcanopy_Cover == 4.5 ~ 5,
                                                         .default = Subcanopy_Cover)) |> 
        dplyr::mutate(Subcanopy_Cover = dplyr::case_when(Subcanopy_Cover > 5 & Subcanopy_Cover <= 20 ~ 1,
                                                         Subcanopy_Cover > 20 & Subcanopy_Cover <= 40 ~ 2,
                                                         Subcanopy_Cover > 40 & Subcanopy_Cover <= 60 ~ 3,
                                                         Subcanopy_Cover > 60 & Subcanopy_Cover <= 80 ~ 4,
                                                         Subcanopy_Cover > 80 & Subcanopy_Cover <= 100 ~ 5,
                                                         .default = Subcanopy_Cover))
    }
    return(df2)
  }) # Close the lapply function
  
  # add group of observers
  tables <- lapply(tables, function(dt){
    dt <- data.table::as.data.table(dt)[, ":="(start_date = as.POSIXct(Date),
                                               end_date = as.POSIXct(Date))]
    dt <- suppressWarnings({ 
      data.table::foverlaps(dt, obs, by.x = c("start_date", "end_date"), type = "within") |> 
        tidytable::select(-start_date, -start_date, -i.start_date, -i.end_date, -group_observers, -end_date) |> 
        tidytable::relocate(ObserverGroup, .before = Weather) |> 
        as.data.frame()
      })
        if("ObserverID" %in% names(dt)){
        dt <- dt |> dplyr::select(-ObserverID)
        return(dt)} else {return(dt)}
    })
  
  # process perpedincular distances
  table[["perp_distance"]] <- table[["perp_distance"]] |> 
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ifelse(. %in% c("NA's", ""), NA, .)),
                  Date = as.Date(Date, format = "%d/%m/%Y")) |> 
    SEQKoalaDataPipeline::fcn_add_date_interval()
  
  # combine survey and perpendicular distances into a list
  master <- append(tables, table["perp_distance"])
  
  return(master)
}


### Spatial representation of line transect data ###
fcn_line_transect_sf_all_detect <- function (
    cols = c(
    "TransectID",
    "SiteID",
    "Date",
    "TrSiteID",
    "Tlength",
    "Number_Sightings",
    "Number_Observers",
    "Start_Eastings",
    "Start_Northings",
    "End_Eastings",
    "End_Northings",
    "geometry"
  ),
  table.list = master,
  table.name = "line_transect")
  {
    state <- fcn_get_state()
    if (state$use_integrated_db) {
      # retrieve line transect table
      db_original <- table.list[[table.name]]
      # stop if not using integrated table
      if (!state$use_integrated_db) {
        stop("use_integrated_db must be set to TRUE to get integrated db in SF format")
      }
      # line transects to sf
      sp <-
        sf::st_read(file.path(state$home_dir, state$gdb_path$total_db),
                    quiet = T) |>
        dplyr::mutate(TransectID = ifelse(is.na(TrnscID), Trns_ID, TrnscID))  |>
        dplyr::select(TransectID)
      if (anyDuplicated(sp$TransectID)) {
        warning(
          sprintf(
            "Number of records with duplicate TransectID detected in %s: %s. Dropping duplicates",
            state$gdb_path$total_db,
            sum(duplicated(sp$TransectID))
          )
        )
        sp <- sp[!duplicated(sp$TransectID), ]
      }

      # join survey data with sf object
      db_joined <- sp |>
        dplyr::inner_join(db_original, by = "TransectID")
      db_not_joined <-
        dplyr::anti_join(db_original, sp, by = "TransectID")
      
      # Create spatial representation if the db is not joined
      db_coords_present <- db_not_joined |>
        dplyr::filter(
          !is.na(Start_Eastings) &
            !is.na(Start_Northings) &
            !is.na(End_Eastings) & !is.na(End_Northings)
        )
      
      db_coords_absent <-
        dplyr::filter(db_not_joined,!(TransectID %in% db_coords_present))
      
      line_transect <- fcn_line_transect_sf(db_coords_present)
      buffer_width <- fcn_get_line_transect_buffer()
      line_transect_buffer <-
        sf::st_buffer(line_transect, endCapStyle = "FLAT", dist = buffer_width)
      db <- rbind(db_joined, line_transect_buffer)
      
      # Recover unjoined records to join at the site level
      site_join_results <-
        fcn_join_koala_survey_sites(db_coords_absent)
      db <- rbind(db, site_join_results$joined_table)
      unjoined_table <- site_join_results$unjoined_table
      
      if (nrow(sp) != nrow(db_original)) {
        warning(
          sprintf(
            "Integrated DB line transects (total=%s): %s joined with total shp, %s created with coordinates, %s joined with site information, %s not joined. Final result: %s records",
            nrow(db_original),
            nrow(db_joined),
            nrow(db_coords_present),
            nrow(site_join_results$joined_table),
            nrow(db_original) - nrow(db),
            nrow(db)
          )
        )
      }
      
      return(db)
    }
  }


### Spatial representation strip transects ###
fcn_strip_transect_sf_all_detect <- function(table.list = master,
                                             table.name = "strip_transect",
                                             date.format = "%d/%m/%Y")
{
  state <- fcn_get_state()
  # retrieve strip transect table
  strip_transects <- table.list[[table.name]]
  
  if (state$use_integrated_db) {
    transect_sf <-
      sf::st_read(file.path(state$home_dir, state$gdb_path$total_db),
                  quiet = T) |>
      dplyr::mutate(TransectID = ifelse(is.na(TrnscID), Trns_ID, TrnscID))  |>
      dplyr::select(TransectID)
    if (anyDuplicated(transect_sf$TransectID)) {
      warning(
        sprintf(
          "Number of records with duplicate TransectID detected in %s: %s. Dropping duplicates",
          state$gdb_path$total_db,
          sum(duplicated(transect_sf$TransectID))
        )
      )
      transect_sf <-
        transect_sf[!duplicated(transect_sf$TransectID),]
    }

    # join survey data
    joined_table <- dplyr::inner_join(transect_sf, strip_transects,
                                      by = c("TransectID"))
    unjoined_table <-
      dplyr::anti_join(strip_transects, transect_sf,
                       by = c("TransectID"))
    site_join_results <- fcn_join_koala_survey_sites(unjoined_table)
    joined_table <-
      rbind(joined_table, site_join_results$joined_table)
    unjoined_table <- site_join_results$unjoined_table
    if (nrow(site_join_results$joined_table) > 0) {
      warning(
        sprintf(
          "Strip Transect: attribute join incomplete, %s transects joined at the site level",
          nrow(site_join_results$joined_table)
        )
      )
    }
    if (nrow(unjoined_table) > 0) {
      warning(
        sprintf(
          "Strip Transect: attribute join incomplete, with %s transects dropped.",
          nrow(unjoined_table)
        )
      )
    }
    return(joined_table)
  } # if (state$use_integrated_db)
}

# Spatial representation of all-of-area searches
fcn_all_of_area_sf_all_detect <- function (table.list = master,
                                           table.name = "uaoa",
                                           date.format = "%d/%m/%Y") 
{
  state <- fcn_get_state()
  # retrieve urban all of area searches table
  all_of_areas <- table.list[[table.name]]
  
  if (state$use_integrated_db) {
    site_sf <-  sf::st_read(file.path(state$home_dir, state$gdb_path$total_db),
                            quiet = T) |>
      dplyr::mutate(TransectID = ifelse(is.na(TrnscID), Trns_ID, TrnscID))  |>
      dplyr::select(TransectID)
    if (anyDuplicated(site_sf$TransectID)) {
      warning(
        sprintf(
          "Number of records with duplicate TransectID detected in %s: %s. Dropping duplicates",
          state$gdb_path$total_db,
          sum(duplicated(site_sf$TransectID))
        )
      )
      site_sf <-
        site_sf[!duplicated(site_sf$TransectID),]
    }

    joined_table <- dplyr::inner_join(site_sf, all_of_areas, 
                                      by = c("TransectID"))
    unjoined_table <- dplyr::anti_join(all_of_areas, site_sf, 
                                       by = c("TransectID"))
    site_join_results <- fcn_join_koala_survey_sites(unjoined_table)
    joined_table <- rbind(joined_table, site_join_results$joined_table)
    unjoined_table <- site_join_results$unjoined_table
    if (nrow(site_join_results$joined_table) > 0) {
      warning(sprintf("Strip Transect: attribute join incomplete, %s transects joined at the site level", 
                      nrow(site_join_results$joined_table)))
    }
  }
  
  if (nrow(unjoined_table) > 0) {
    warning(paste0("Number of rows not joined: ", nrow(unjoined_table)))
  }
  study_area <- fcn_get_study_area()
  joined_table <- sf::st_intersection(joined_table, study_area) |> 
    dplyr::select(-AREA_HA, SHAPE_Length, SHAPE_Area)
  return(joined_table)
}


### Combine sf objects into a list ###
fcn_all_tables_sf_detect <- function (table_names = c("line_transect", 
                                                      "strip_transect", 
                                                      "uaoa")) 
{
  master_sf <- list()
  if ("line_transect" %in% table_names) {
    master_sf$line_transect = fcn_line_transect_sf_all_detect()}
  if ("strip_transect" %in% table_names) {
    master_sf$strip_transect = fcn_strip_transect_sf_all_detect()}
  if ("uaoa" %in% table_names){
    master_sf$uaoa = fcn_all_of_area_sf_all_detect()}
  master_sf <- lapply(master_sf, fcn_add_date_interval)
  return(master_sf)
}

### Covariate data frame ###
fcn_covariate_layer_df_detect <- function(layer = NULL) {
  # Extract name from string containing dates
  fcn_get_covariate_name <- function(input) {
    pattern <- ".*\\d{6}.tif$"
    if (grepl(pattern, input)) {
      output <- gsub("\\d{6}.tif$", "", input)
    } else {
      output <- gsub(".tif$", "", input)
    }
    return(output)
  }
  
  state <- SEQKoalaDataPipeline::fcn_get_state()
  match_method <- state$covariate_time_match
    
  covariate_description <- paste0(SEQKoalaDataPipeline::fcn_get_raster_path()$covariates, "\\..\\covariate_descriptions.csv") |>
    readr::read_csv(show_col_types = FALSE) |>
    dplyr::mutate(name = substr(Code, 0, 5)) |>
    dplyr::rename(static_dynamic = `Static/Dynamic`, continuous_discrete = `Continuous/Discrete`) |>
    dplyr::select(name, static_dynamic, continuous_discrete, Proportions )
  
  covariate_description <- covariate_description[!duplicated(covariate_description$name),]
  
  constant_covariates <- data.frame(filename = SEQKoalaDataPipeline::fcn_list_covariate_layers_constant())
  temporal_covariates <- data.frame(filename = SEQKoalaDataPipeline::fcn_list_covariate_layers_temporal())
  temporal_covariates <- temporal_covariates |>
    dplyr::mutate(date = as.numeric(as.numeric(gsub(".*[^0-9]([0-9]{6})[^0-9]*", "\\1", filename)))) |>
    dplyr::mutate(date = ifelse(is.na(date), NA, paste0("X", date))) |>
    dplyr::mutate(fullname = sub('\\.tif$', '', filename))
  
  df <- dplyr::bind_rows(list(constant = constant_covariates,
                              temporal = temporal_covariates),
                         .id = 'type') |>
    dplyr::mutate(name = sapply(filename, fcn_get_covariate_name)) |>
    dplyr::mutate(match_method = match_method[name]) |>
    dplyr::mutate(fullname = sub('\\.tif$', '', filename))
  
  # Exclude files that have a processing error, less than 100 bytes
  file_sizes <- file.info(paste0(state$home_dir, "\\", state$raster_path,'\\', df$filename))$size
  low_filesize <- 100
  if (any(file_sizes < low_filesize)) {
    invalid_files <- df$filename[file_sizes < low_filesize]
    warning(sprintf("Raster(s) %s appears to be less than 100bytes and therefore removed.", paste(invalid_files, collapse=", ")))
    df <- df[file_sizes > low_filesize,]
  }
  
  if (!is.null(layer)) {
    df <- df |>
      filter(name == layer)
  }
  
  df <- df |>
    dplyr::mutate(join_name = substr(name, 0, 5)) |>
    dplyr::inner_join(covariate_description, by = dplyr::join_by('join_name' == 'name')) |>
    dplyr::select(-join_name)
  
  return(df)
}

###
fcn_covariate_interval_mean_detect <- function (date, cov_name, get_df = TRUE) 
{
  cov_layer_df <- fcn_covariate_layer_df_detect()
  cov_layer_df_name <- cov_layer_df[cov_layer_df$name == cov_name, 
  ]
  cov_categorical <- any(cov_layer_df[cov_layer_df$name == 
                                        cov_name, ] == "Categorical")
  cov_dates <- cov_layer_df_name$date
  cov_dates_num <- gsub(".*[^0-9]([0-9]{6})[^0-9]*", "\\1", 
                        cov_dates) |> as.numeric() |> lubridate::ym()
  within_interval <- lubridate::`%within%`(cov_dates_num, date$interval)
  if (sum(within_interval) > 0) {
    cov_layer_list <- cov_layer_df_name[within_interval, 
                                        "filename"]
  }
  else {
    date_diff <- abs(cov_dates_num - date$middle_date)
    if (cov_name %in% c("hctmn", "hctma", "hcpre")) {
      months_in_interval <- function(start_date, end_date) {
        lubridate::month(seq.Date(from = as.Date(start_date, 
                                                 tz = Sys.timezone()), to = as.Date(end_date, 
                                                                                    tz = Sys.timezone()), by = "day")) |> unique()
      }
      date_interval_months <- months_in_interval(date$start_date, 
                                                 date$end_date)
      date_diff[!(lubridate::month(cov_dates_num) %in% 
                    date_interval_months)] <- NA
    }
    nearest_date <- which.min(date_diff)
    cov_layer_list <- cov_layer_df_name[nearest_date, "filename"]
  }
  cov_temporal <- SEQKoalaDataPipeline::fcn_extract_covariate_grid(cov_layer_list)
  if (length(cov_layer_list) > 1) {
    if (cov_categorical) {
      cov_temporal_mean <- terra::modal(cov_temporal[[2:terra::nlyr(cov_temporal)]], 
                                        ties = "first")
      cov_temporal_mean <- terra::round(cov_temporal_mean)
    }
    else {
      cov_temporal_mean <- terra::mean(cov_temporal[[2:terra::nlyr(cov_temporal)]])
    }
    cov_temporal <- c(cov_temporal["GridID"], cov_temporal_mean)
  }
  names(cov_temporal) <- c("GridID", cov_name)
  calc_proportion <- any(cov_layer_df[cov_layer_df$name == 
                                        cov_name, "Proportions"] == "Yes")
  if (calc_proportion) {
    cov_temporal_prop <- SEQKoalaDataPipeline::fcn_extract_covariate_grid(cov_layer_list, 
                                                    proportion = TRUE)
    name_list <- names(cov_temporal_prop[[2:terra::nlyr(cov_temporal_prop)]])
    cov_lyr_names <- stringr::str_split_i(name_list, "_", 
                                          1)
    cov_values <- stringr::str_split_i(name_list, "_", 2)
    cov_temporal_mean <- lapply(unique(cov_values), function(x) terra::mean(cov_temporal_prop[name_list[x == 
                                                                                                          cov_values]])) |> terra::rast()
    names(cov_temporal_mean) <- paste(cov_name, unique(cov_values), 
                                      sep = "_")
    cov_temporal <- c(cov_temporal, cov_temporal_mean)
  }
  if (get_df) {
    cov_temporal_df <- SEQKoalaDataPipeline::fcn_cov_grid_df(cov_temporal)
    return(cov_temporal_df)
  }
  else {
    return(cov_temporal)
  }
}


### extract temporal covariates ###
fcn_extract_cov_date_detect <- function(d, time_lag = NULL) {
  cov_layer_df <- fcn_covariate_layer_df_detect()
  cov_temporal_names <- cov_layer_df[cov_layer_df$type == 'temporal','name'] |>
    unique()
  
  print(sprintf("Processing covariates for date starting %s", d$start_date))
  
  # Apply for each of the covariate names
  covs <- purrr::map(cov_temporal_names, \(x) fcn_covariate_interval_mean_detect(d, x))
  cov_comb <- purrr::reduce(covs, dplyr::full_join, by = "GridID")
  
  # Compute lagged variables and join back to the dataframe
  if (!is.null(time_lag)) {
    for (lag in time_lag) {
      d_lag <- SEQKoalaDataPipeline::fcn_shift_months(d, -lag) # date object lagged
      covs_lag <- purrr::map(cov_temporal_names, ~ fcn_covariate_interval_detect(d_lag, .))
      cov_lag_comb <- purrr::reduce(covs_lag, dplyr::inner_join, by = "GridID")
      colnames(cov_lag_comb)[2:ncol(cov_lag_comb)] <- paste0(colnames(cov_lag_comb)[2:ncol(cov_lag_comb)],
                                                             "_",
                                                             lag, "mths")
      cov_comb <- dplyr::inner_join(cov_comb, cov_lag_comb, by = "GridID")
    }
  }
  return(cov_comb)
}


### develop covariate array ###
fcn_cov_array_detect <- function (cov_type = "temporal", dates = d, time_lag = NULL, write_path = paste0(out_dir, "/cov_raster")) 
{
  compute_constant <- cov_type %in% c("constant", "both")
  compute_temporal <- cov_type %in% c("temporal", "both")
  if (is.null(dates)) {
    dates <- SEQKoalaDataPipeline::fcn_get_date_intervals()
  }
  output <- list()
  cov_layer_df <- fcn_covariate_layer_df_detect()
  if (compute_constant) {
    cov_constant_names <- cov_layer_df[cov_layer_df$type == 
                                         "constant", "filename"]
    cov_constant <- SEQKoalaDataPipeline::fcn_extract_covariate_grid(cov_constant_names)
    cov_constant_df <- SEQKoalaDataPipeline::fcn_cov_grid_df(cov_constant)
    output$cov_constant <- cov_constant_df
    if (!is.null(write_path)) {
      saveRDS(cov_constant_df, file = paste0(write_path, 
                                             "\\cov_constant_array.rds"))
    }
    rm("cov_constant_df")
    gc()
  }
  if (compute_temporal) {
    dates <- list(dates)
    if (!is.null(write_path)) {
      lapply(dates, function(d) {
        cov_temporal <- fcn_extract_cov_date_detect(d, time_lag = NULL)
        readr::write_rds(cov_temporal, file = paste0(write_path, 
                                                     "\\cov_temporal_", d$id, ".rds"))
        rm(cov_temporal)
      })
    }
    else {
      cov_temporal <- lapply(dates, function(d) {
        SEQKoalaDataPipeline::fcn_extract_cov_date(d, time_lag = NULL)
      })
      cov_temporal_array <- do.call(abind::abind, c(cov_temporal, 
                                                    list(along = 3)))
      output$cov_temporal <- cov_temporal_array
    }
  }
  return(output)
}

### grid fraction ###
fcn_all_transect_grid_fractions <- function (buffer = c(0), keep_all = FALSE, grid_id_vec = NULL) 
{
  fcn_extract_raster_buffer <- function(df, fishnet, buffer = 0) {
    if (buffer > 0) {
      df <- sf::st_buffer(df, dist = buffer, joinStyle = 'ROUND')
    }
    res <- fcn_mixed_extract_raster(fishnet, df)
    return(res)
  }
  
  fishnet <- fcn_get_grid()
  master <- fcn_all_tables_sf_detect()
  buffer <- sort(buffer)
  master_grid <- lapply(master, function(df) {
    res <- fcn_extract_raster_buffer(df, fishnet, 0)
    if (!is.null(grid_id_vec)) {
      res <- res %>% dplyr::filter(GridID %in% grid_id_vec) %>% 
        dplyr::group_by(TransectID) %>% dplyr::mutate(fraction = fraction/sum(fraction)) %>% 
        dplyr::ungroup()
    }
    if (length(buffer) == 1) {
      var_name <- c("fraction")
    }
    else {
      var_name <- paste0("fraction_", buffer)
    }
    res <- res %>% dplyr::rename(`:=`(!!var_name[1], fraction))
    if (length(buffer) > 1) {
      for (i in 2:length(buffer)) {
        res_buffer <- fcn_extract_raster_buffer(df, fishnet, 
                                                buffer[i]) %>% dplyr::rename(`:=`(!!var_name[i], 
                                                                                  fraction))
        res <- res %>% dplyr::select(TransectID, GridID, 
                                     dplyr::all_of(var_name[1:(i - 1)])) %>% dplyr::right_join(res_buffer, 
                                                                                               by = c("TransectID", "GridID"))
      }
      res <- res %>% dplyr::mutate_at(var_name, ~replace(., 
                                                         is.na(.), 0))
    }
    if (!keep_all) {
      res <- res %>% dplyr::select("TransectID", "GridID", 
                                   "Date", dplyr::all_of(var_name))
      res <- fcn_add_date_interval(res)
    }
    return(res)
  })
  return(master_grid)
}
