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


