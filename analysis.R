# load libraries
library(tidyverse)
library(nimble)
library(lubridate)
library(coda)
library(MCMCvis)
library(abind)

# read utility functions
source("functions.R")

# load input data
Surveys <- readRDS("input/survey_data/master.rds")
GridFrac <- readRDS("input/survey_data/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data/cov_temporal_array_surveylocations.rds")
Adj2km_Queen <- readRDS("input/survey_data/adj_data_queen.rds")
Adj2km_Rook <- readRDS("input/survey_data/adj_data_rook.rds")
DateIntervals <- read_csv("input/survey_data/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))

# set up data for nimble models

# set the date range for analysis
# get first date - either date of first survey or user specified date
#FirstDate <- date("2015-01-01")
FirstDate <- date("2013-01-01")
#FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
# get last date - either date of last survey or user specified date
#LastDate <- date("2020-12-31")
#LastDate <- date("2013-01-01")
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
# get first and last date index values
FirstDateID <- filter(DateIntervals, start_date <= FirstDate & end_date >= FirstDate) %>% pull(TimePeriodID)
LastDateID <- filter(DateIntervals, start_date <= LastDate & end_date >= LastDate) %>% pull(TimePeriodID)

# set the Lag for the effect of predictors on growth rate in terms of number of 6-monthly time steps (can be maximum of 2 currently)
Lag <- 0

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
GridFrac <- GridFrac[which((GridFrac$TransectID %in% Surveys$line_transect$TransectID) | (GridFrac$TransectID %in% Surveys$strip_transect$TransectID) | (GridFrac$TransectID %in% Surveys$uaoa$TransectID)),]
# covariates  - selecting only the grids in transects for the time period selected
# and selecting only the time periods needed for the temporally variable covariates
CovConsSurv <- CovConsSurv[which((CovConsSurv$GridID %in% GridFrac$GridID)),]
CovTempSurv <- CovTempSurv[which((CovTempSurv[,"GridID", 1] %in% GridFrac$GridID)),, (FirstDateID - Lag):LastDateID]

# set up data for small grids (only for those within the surveyed areas)
NSGrids <- length(CovConsSurv$GridID)
LGridID <- left_join(as.data.frame(CovConsSurv$GridID), Adj2km_Queen$grid_lookup, by = c("CovConsSurv$GridID" = "GridID"))$SpGridID

# predictors for initial density for exponential population growth model

# persistent green
X_hhpgr <- as.vector(CovTempSurv[,"hhpgr", Lag + 1]) %>% scale() %>% as.vector()
# koala habitat
X_hhkha <- as.vector(CovTempSurv[,"hhkha", Lag + 1]) %>% as.factor()

# collate initial density predictor variables

# compile into a tibble
X <- tibble(hhpgr = X_hhpgr, hhkha = X_hhkha)
# get the design matrix
X <- model.matrix(~ hhpgr + hhkha, model.frame(~ hhpgr + hhkha, as.data.frame(X), na.action = "na.pass"))
# remove the intercept term
X <- X[, 2:ncol(X)] %>% as.data.frame()

# get number of variables in X
NX <- ncol(X)

# predictors for growth rate (to the next time step ie., from t-1 to t) for the exponential population growth model
# or density for the spatio-temporal regression model

# time
Y_htime <- as.vector(CovTempSurv[,"htime",]) %>% scale() %>% as.vector()
# season
Y_hseas <- as.vector(CovTempSurv[,"hseas",]) %>% as.vector()
# persistent green
Y_hhpgr <- as.vector(CovTempSurv[,"hhpgr",]) %>% scale() %>% as.vector()
# koala habitat
Y_hhkha <- as.vector(CovTempSurv[,"hhkha",]) %>% as.factor()

# collate predictor variables
# NOTE - time must be the first variable in the design matrix

# compile into a tibble
Y_temp <- tibble(htime = Y_htime, hseas = Y_hseas, hhpgr = Y_hhpgr, hhkha = Y_hhkha)
# get the design matrix
Y_temp <- model.matrix(~ hseas + hhpgr, model.frame(~ hseas + hhpgr, as.data.frame(Y_temp), na.action = "na.pass"))
# remove the intercept term
Y_temp <- Y_temp[, 2:ncol(Y_temp)] %>% as.data.frame()

# get number of variables in Y
NY <- ncol(Y_temp)

# create 3D matrix of covariates
Y <- array(NA, dim = c(dim(CovTempSurv)[1], dim(CovTempSurv)[3], NY))
for (i in 1:NY) {
	Y[,,i] <- matrix(Y_temp[,i], nrow = dim(CovTempSurv)[1], ncol = dim(CovTempSurv)[3])
}

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

# line transect data

# get the number of perependicular distances
NPDists <- nrow(Surveys$perp_distance) %>% as.vector()

# get perpendicular distances
PDists <- Surveys$perp_distance$Perp_Dist %>% as.vector()
PDists[is.na(PDists)] <- 10 #EDIT OUT WHEN FIXED

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

# reset first and last date IDs
LastDateID <- LastDateID - FirstDateID + 1 + Lag
FirstDateID <- Lag + 1

# set up data for spatial CAR process
AdjS <- Adj2km_Queen$adjacencyList$adj
WeightsAdjS <- Adj2km_Queen$adjacencyList$weights
NumAdjS <- Adj2km_Queen$adjacencyList$num
NLGrids <- length(Adj2km_Queen$adjacencyList$num)
NLGridAdjs <- length(Adj2km_Queen$adjacencyList$adj)

# set up data for temporal CAR process (1st order with equal weights)
# annual time steps
AdjT <- c()
WeightsAdjT <- c()
NumAdjT <- c()
WeightsAdjT[1] <- 1
AdjT[1] <- 2
NumAdjT[1] <- 1
for(i in 2:(floor(((LastDateID - FirstDateID + 1) / 2) - 0.01))) {
  	WeightsAdjT[2 + (i - 2) * 2] <- 1
	AdjT[2 + (i - 2) * 2] <- i - 1
  	WeightsAdjT[3 + (i - 2) * 2] <- 1
	AdjT[3 + (i - 2) * 2] <- i + 1
	NumAdjT[i] <- 2
}
for(i in (floor(((LastDateID - FirstDateID + 1) / 2) - 0.01) + 1):(floor(((LastDateID - FirstDateID + 1) / 2) - 0.01) + 1)) {
  	WeightsAdjT[(i - 2) * 2 + 2] <- 1
	AdjT[(i - 2) * 2 + 2] <- i - 1
	NumAdjT[i] <- 1
}
NTime <- length(NumAdjT)
NTimeAdjs <- length(AdjT)

# set up data for temporal CAR process (2nd order with weights acording)
# annual time steps
AdjT2 <- c()
WeightsAdjT2 <- c()
NumAdjT2 <- c()
WeightsAdjT2[1] <- 1
AdjT2[1] <- 2
NumAdjT2[1] <- 1
for(i in 2:(floor(((LastDateID - FirstDateID + 1) / 2) - 0.01))) {
  	WeightsAdjT2[2 + (i - 2) * 2] <- 1
	Adj2T[2 + (i - 2) * 2] <- i - 1
  	WeightsAdjT2[3 + (i - 2) * 2] <- 1
	AdjT2[3 + (i - 2) * 2] <- i + 1
	NumAdjT2[i] <- 2
}
for(i in (floor(((LastDateID - FirstDateID + 1) / 2) - 0.01) + 1):(floor(((LastDateID - FirstDateID + 1) / 2) - 0.01) + 1)) {
  	WeightsAdjT2[(i - 2) * 2 + 2] <- 1
	AdjT2[(i - 2) * 2 + 2] <- i - 1
	NumAdjT2[i] <- 1
}
NTime2 <- length(NumAdjT2)
NTimeAdjs2 <- length(AdjT2)

# get nimble code for spatio-temporal regression model
source("nimble_spt.R")

# spatio-temporal regression model fitting

NimbleConstsSPT <- list(NLGrids = NLGrids, NLGridAdjs = NLGridAdjs, AdjS = AdjS, WeightsAdjS = WeightsAdjS, NumAdjS = NumAdjS, NTime = NTime, NTimeAdjs = NTimeAdjs, AdjT = AdjT, WeightsAdjT = WeightsAdjT, NumAdjT = NumAdjT, NTime2 = NTime2, NTimeAdjs2 = NTimeAdjs2, AdjT2 = AdjT2, WeightsAdjT2 = WeightsAdjT2, NumAdjT2 = NumAdjT2, NSGrids = NSGrids, LGridID = LGridID, FirstDateID = FirstDateID, LastDateID = LastDateID, Lag = Lag, NY = NY, NStrips = NStrips, SGridsStartStrip = SGridsStartStrip, SGridsEndStrip = SGridsEndStrip, SGridIDsStrip = SGridIDsStrip, SGridFracsStrip = SGridFracsStrip, NumObsStrip = NumObsStrip, AreaStrip = AreaStrip, TimeIDStrip = TimeIDStrip, NAoAs = NAoAs, SGridsStartAoA = SGridsStartAoA, SGridsEndAoA = SGridsEndAoA, SGridIDsAoA = SGridIDsAoA, SGridFracsAoA = SGridFracsAoA, NumObsAoA = NumObsAoA, AreaAoA = AreaAoA, TimeIDAoA = TimeIDAoA, NPDists = NPDists, NLines = NLines, SGridsStartLine = SGridsStartLine, SGridsEndLine = SGridsEndLine, SGridIDsLine = SGridIDsLine, SGridFracsLine = SGridFracsLine, LengthLine = LengthLine, TimeIDLine = TimeIDLine, PI = pi, NMaxSGridsAoA = NMaxSGridsAoA, NMaxSGridsStrip = NMaxSGridsStrip, NMaxSGridsLine = NMaxSGridsLine)

NimbleDataSPT <- list(Y = Y, CntStrip = CntStrip, CntAoA = CntAoA, PDists = PDists, CntLine = CntLine)

NimbleInitsSPT <- list(sigma_sd = 0.5, sd = rep(0, NLGrids), sigma_td = 0.5, td = rep(0, NTime), beta_d = rep(0, NY), sigma_d = 0.5, d = matrix(10, nrow = NSGrids, ncol = (LastDateID - FirstDateID + 1)), pStrip = 0.5, pAoA = 0.5, sigma_hn = 1)

N.iter <- 5000
N.burnin <- 0
N.chains <- 1

NimbleModel <- nimbleModel(code = NimbleCodeSPT, constants = NimbleConstsSPT, data = NimbleDataSPT, inits = NimbleInitsSPT, calculate = TRUE)

CNimbleModel <- compileNimble(NimbleModel)

NimbleModelConf <- configureMCMC(NimbleModel, monitors = c("beta_d", "pStrip", "pAoA", "sigma_hn", "mean_sd", "mean_d", "td"))

NimbleModelMCMC <- buildMCMC(NimbleModelConf)

CNimbleModelMCMC <- compileNimble(NimbleModelMCMC, project = NimbleModel)

Samples <- runMCMC(mcmc = CNimbleModelMCMC, niter = N.iter, nburnin = N.burnin, nchains = N.chains)

MCMCsummary(object = Samples, round = 2)

MCMCplot(object = Samples, params = c("beta_di", "beta_l"))

MCMCtrace(object = Samples, pdf = FALSE, ind = TRUE, params = c("mean_d"))

# get nimble code for exponential population growth model
source("nimble_exp.R")

# exponential population growth model fitting

NimbleConsts <- list(NLGrids = NLGrids, NLGridAdjs = NLGridAdjs, AdjS = AdjS, WeightsAdjS = WeightsAdjS, NumAdjS = NumAdjS, NSGrids = NSGrids, LGridID = LGridID, FirstDateID = FirstDateID, LastDateID = LastDateID, Lag = Lag, NX = NX, NY = NY, NStrips = NStrips, SGridsStartStrip = SGridsStartStrip, SGridsEndStrip = SGridsEndStrip, SGridIDsStrip = SGridIDsStrip, SGridFracsStrip = SGridFracsStrip, NumObsStrip = NumObsStrip, AreaStrip = AreaStrip, TimeIDStrip = TimeIDStrip, NAoAs = NAoAs, SGridsStartAoA = SGridsStartAoA, SGridsEndAoA = SGridsEndAoA, SGridIDsAoA = SGridIDsAoA, SGridFracsAoA = SGridFracsAoA, NumObsAoA = NumObsAoA, AreaAoA = AreaAoA, TimeIDAoA = TimeIDAoA, NPDists = NPDists, NLines = NLines, SGridsStartLine = SGridsStartLine, SGridsEndLine = SGridsEndLine, SGridIDsLine = SGridIDsLine, SGridFracsLine = SGridFracsLine, LengthLine = LengthLine, TimeIDLine = TimeIDLine, PI = pi, NMaxSGridsAoA = NMaxSGridsAoA, NMaxSGridsStrip = NMaxSGridsStrip, NMaxSGridsLine = NMaxSGridsLine)

NimbleData <- list(X = X, Y = Y, CntStrip = CntStrip, CntAoA = CntAoA, PDists = PDists, CntLine = CntLine)

NimbleInits <- list(sigma_sdi = 0.5, sdi = rep(0, NLGrids), sigma_sl = 0.5, sl = rep(0, NLGrids), sigma_tl = 0.5, tl = rep(0, LastDateID - FirstDateID), beta_di = rep(0, NX), beta_l = rep(0, NY), sigma_di = 0.5, sigma_l = 0.5, di = rep(10, NSGrids), lambda = matrix(1, nrow = NSGrids, ncol = (LastDateID - FirstDateID)), pStrip = 0.5, pAoA = 0.5, sigma_hn = 1)

N.iter <- 1000
N.burnin <- 0
N.chains <- 1

NimbleModel <- nimbleModel(code = NimbleCode, constants = NimbleConsts, data = NimbleData, inits = NimbleInits, calculate = TRUE)

CNimbleModel <- compileNimble(NimbleModel)

NimbleModelConf <- configureMCMC(NimbleModel)

NimbleModelMCMC <- buildMCMC(NimbleModelConf)

CNimbleModelMCMC <- compileNimble(NimbleModelMCMC, project = NimbleModel)

Samples <- runMCMC(mcmc = CNimbleModelMCMC, niter = N.iter, nburnin = N.burnin, nchains = N.chains)

MCMCsummary(object = Samples, round = 2)

MCMCplot(object = Samples, params = c("beta_di", "beta_l"))

MCMCtrace(object = Samples, pdf = FALSE, ind = TRUE, params = c("beta_l"))


















