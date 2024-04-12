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
#GridFrac <- readRDS("input/survey_data/grid_fractions_complete_cov.rds")
#CovConsSurv <- readRDS("input/survey_data/cov_constant_array_complete_cov_surveylocations.rds")
#CovTempSurv <- readRDS("input/survey_data/cov_temporal_array_complete_cov_surveylocations.rds")
GridFrac <- readRDS("input/survey_data/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data/cov_temporal_array_surveylocations.rds")
Adj2km_Queen <- readRDS("input/survey_data/adj_data_queen.rds")
Adj2km_Rook <- readRDS("input/survey_data/adj_data_rook.rds")
DateIntervals <- read_csv("input/survey_data/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))

# set up data for nimble models

# set the date range for analysis
# get first date - either date of first survey or user specified date
FirstDate <- date("2015-01-01")
#FirstDate <- date("2013-01-01")
#FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
# get last date - either date of last survey or user specified date
LastDate <- date("2020-12-31")
#LastDate <- date("2013-01-01")
#LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
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

# set up data for large grids for spatial autocorrelation
NLGrids <- length(Adj2km_Queen$adjacencyList$num)
NLGridAdjs <- length(Adj2km_Queen$adjacencyList$adj)
Adj <- Adj2km_Queen$adjacencyList$adj
WeightsAdj <- Adj2km_Queen$adjacencyList$weights
NumAdj <- Adj2km_Queen$adjacencyList$num

# set up data for small grids (only for those within the surveyed areas)
NSGrids <- length(CovConsSurv$GridID)
LGridID <- left_join(as.data.frame(CovConsSurv$GridID), Adj2km_Queen$grid_lookup, by = c("CovConsSurv$GridID" = "GridID"))$SpGridID

# predictors for initial density

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

# predictors for growth rate (to the next time step ie., from t-1 to t)

# time
Y_htime <- as.vector(CovTempSurv[,"htime",]) %>% scale() %>% as.vector()
# season
Y_hseas <- as.vector(CovTempSurv[,"hseas",]) %>% as.vector()
# persistent green
Y_hhpgr <- as.vector(CovTempSurv[,"hhpgr",]) %>% scale() %>% as.vector()
# koala habitat
Y_hhkha <- as.vector(CovTempSurv[,"hhkha",]) %>% as.factor()

# collate predictor variables

# compile into a tibble
Y_temp <- tibble(htime = Y_htime, hseas = Y_hseas, hhpgr = Y_hhpgr, hhkha = Y_hhkha)
# get the design matrix
Y_temp <- model.matrix(~ htime + hseas + hhpgr, model.frame(~ htime + hseas + hhpgr, as.data.frame(Y_temp), na.action = "na.pass"))
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

# nimble code

NimbleCode <- nimbleCode({
	# NLGrids = number of large grids
	# NLGridAdjs = number of large grid adjacent pairs
	# Adj = adjacency matrix for large grids of length NLGridAdjs
	# WeightsAdj = weights for the adjacency matrix for large grids (2km) of length NLGridAdjs
	# NumAdj = number of neigbouring locations for each large grid (2km) of length NLGrids
	# NSGrids = number of small grids
	# LGridID = the larget grid ID for each small grid
	# FirstDateID = ID for the time step of the first survey
	# LastDateID = ID for the time step of the last survey
	# Lag = number of time steps lag for the predictors of growth rate
	# X = static predictors for initial koala density
	# Y = time-varying predictors for growth rate
	# NX = number of predictors for initial koala density
	# NY = number of time-varying predictors for growth rate
	# NStrips = number of strip transects
	# SGridsStartStrip = start index for small grids for each strip transect
    # SGridsEndStrip = end index for small grids for each strip transect
    # SGridIDsStrip = small grid IDs for each strip transect
	# SGridFracsStrip = small grid fractions for each strip transect
	# IDStrip = small grid IDs for each strip transect
	# FracStrip = small grid fractions for each strip transect
	# NMaxSGridsStrip = maximum number of small grids in any strip transect
	# NumObsStrip = number of observers for each strip transect
	# AreaStrip = area of each strip transect (in hectares)
	# TimeIDStrip = time step ID for each strip transect
	# CntStrip = koala count for each strip transect
	# NAoAs = number of all of area searches
	# SGridsStartAoA = start index for small grids for each all of area search
    # SGridsEndAoA = end index for small grids for each all of area search
    # SGridIDsAoA = small grid IDs for each all of area search
	# SGridFracsAoA = small grid fractions for each all of area search
	# IDAoA = small grid IDs for each all of area search
	# FracAoA = small grid fractions for each all of area search
	# NMaxSGridsAoA = maximum number of small grids in any all of areas search
	# NumObsAoA = number of observers for each all of area search
	# AreaAoA = area of each all of area search (in hectares)
	# TimeIDAoA = time step ID for each all of area search
	# CntSAoA = koala count for each all of area search
	# NPDists = number of perpendicular distances
	# PDists = perpendicular distances
	# NLines = number of line transects
	# SGridsStartLine = start index for small grids for each line transect
    # SGridsEndLine = end index for small grids for each line transect
    # SGridIDsLine = small grid IDs for each line transect
	# SGridFracsLine = small grid fractions for each line transect
	# IDLine = small grid IDs for each line transect
	# FracLine = small grid fractions for each line transect
	# NMaxSGridsLine = maximum number of small grids in any line transect
	# LengthLine = length of each line transect (in metres)
	# TimeIDLine = time step ID for each line transect
	# CntLine = koala count for each line transect
	# PI = pi

	# process model

	# large grid spatially correlated random-effect for initial densities (based on intrisic CAR model)
	# note here weights are taken to be 1
	sdi[1:NLGrids] ~ dcar_normal(adj = Adj[1:NLGridAdjs], weights = WeightsAdj[1:NLGridAdjs], num = NumAdj[1:NLGrids], tau_sdi)

	# large grid spatially correlated random-effect for growth rate (based on intrisic CAR model)
	# note here weights are taken to be 1
	sl[1:NLGrids] ~ dcar_normal(adj = Adj[1:NLGridAdjs], weights = WeightsAdj[1:NLGridAdjs], num = NumAdj[1:NLGrids], tau_sl)

	# loop through small grids
	for (i in 1:NSGrids) {
		# predictors for initial density incorporating spatially structured (at the large grid scale) and unstructured (at the small grid scale) stochasticity
		mu_di[i] <- exp(sdi[LGridID[i]] + inprod(beta_di[1:NX], X[i, 1:NX]))
		di[i] ~ dgamma(mean = mu_di[i] , sd = sigma_di)

		d[i, 1] <- di[i]

		# loop through remaining time steps
		for (t in 2:(LastDateID - FirstDateID + 1)) {
			# predictor for density including time-dependent variables, unstructured (at the small grid scale) stochasticity
			# and spatially structured (at the large grid scale) stochasticity
			mu_l[i, t - 1] <- exp(sl[LGridID[i]] + inprod(beta_l[1:NY], Y[i, FirstDateID - 1 + t - 1 - Lag, 1:NY]))

			lambda[i, t - 1] ~ dgamma(mean = mu_l[i, t - 1], sd = sigma_l)

			d[i, t] <- d[i, t - 1] * lambda[i, t - 1]
		}
	}

	# observation models

	# likelihoods

	# strip transects
	for (i in 1:NStrips) {

		# get true density estimate for the strip transect
		for (j in SGridsStartStrip[i]:SGridsEndStrip[i]) {
			dFracStrip[i, j - SGridsStartStrip[i] + 1] <- d[SGridIDsStrip[j], TimeIDStrip[i] - FirstDateID + 1] * SGridFracsStrip[j]
		}
		for (j in ((SGridsEndStrip[i] - SGridsStartStrip[i] + 2):(NMaxSGridsStrip + 1))) {
			dFracStrip[i, j] <- 0
		}
		dStrip[i] <- sum(dFracStrip[i, 1:NMaxSGridsStrip])

		# get true abundance estimate for the strip transect
		aStrip[i] <- round(dStrip[i] * AreaStrip[i])

		# likelihood for counts
		CntStrip[i] ~ dbin(1 - ((1 - pStrip) ^ NumObsStrip[i]), aStrip[i])
	}

	# all of area searches
	for (i in 1:NAoAs) {

		# get true density estimate for the all of area search
		for (j in SGridsStartAoA[i]:SGridsEndAoA[i]) {
			dFracAoA[i, j - SGridsStartAoA[i] + 1] <- d[SGridIDsAoA[j], TimeIDAoA[i] - FirstDateID + 1] * SGridFracsAoA[j]
		}
		for (j in ((SGridsEndAoA[i] - SGridsStartAoA[i] + 2):(NMaxSGridsAoA + 1))) {
			dFracAoA[i, j] <- 0
		}
		dAoA[i] <- sum(dFracAoA[i, 1:NMaxSGridsAoA])

		# get true abundance estimate for the strip transect
		aAoA[i] <- round(dAoA[i] * AreaAoA[i])

		# likelihood for counts
		CntAoA[i] ~ dbin(1 - ((1 - pAoA) ^ NumObsAoA[i]), aAoA[i])
	}

	# line transects

	# perpendicular distances
	for (i in 1:NPDists) {
		# likelihood function for half-normal detection function
		PDists[i] ~ T(dnorm(0, sd = sigma_hn), 0, )
	}

	# get the closed form solution for f0
	# valid when detection function is half-normal, detection on the line is certain, distances are not truncated, and detections are not grouped
	f0 <- sqrt(2 / (PI * (sigma_hn ^ 2)))

	# loop through line transects
	for (i in 1:NLines) {
		# get true density estimate for the line transect
		for (j in SGridsStartLine[i]:SGridsEndLine[i]) {
			dFracLine[i, j - SGridsStartLine[i] + 1] <- d[SGridIDsLine[j], TimeIDLine[i] - FirstDateID + 1] * SGridFracsLine[j]
		}
		for (j in ((SGridsEndLine[i] - SGridsStartLine[i] + 2):(NMaxSGridsLine + 1))) {
			dFracLine[i, j] <- 0
		}
		dLine[i] <- sum(dFracLine[i, 1:NMaxSGridsLine])

		# likelihood for counts
		# density divided by 10,000 to ensure desities are per hectare rather than per square metre
		CntLine[i] ~ dpois((dLine[i] * 2 * LengthLine[i]) / (f0 * 10000))
	}

	# priors
	tau_sdi <- sigma_sdi ^ -2
	sigma_sdi ~ dunif(0, 10)
	tau_sl <- sigma_sl ^ -2
	sigma_sl ~ dunif(0, 10)
	for (i in 1:NX) {
		beta_di[i] ~ dnorm(0, sd = 100)
	}
	for (i in 1:NY) {
		beta_l[i] ~ dnorm(0, sd = 100)
	}
	sigma_di ~ dunif(0, 10)
	sigma_l ~ dunif(0, 10)
	pStrip ~ dunif(0,1)
	pAoA ~ dunif(0,1)
	sigma_hn ~ dunif(0, 100)
})

NimbleConsts <- list(NLGrids = NLGrids, NLGridAdjs = NLGridAdjs, Adj = Adj, WeightsAdj = WeightsAdj, NumAdj = NumAdj, NSGrids = NSGrids, LGridID = LGridID, FirstDateID = FirstDateID, LastDateID = LastDateID, Lag = Lag, NX = NX, NY = NY, NStrips = NStrips, SGridsStartStrip = SGridsStartStrip, SGridsEndStrip = SGridsEndStrip, SGridIDsStrip = SGridIDsStrip, SGridFracsStrip = SGridFracsStrip, NumObsStrip = NumObsStrip, AreaStrip = AreaStrip, TimeIDStrip = TimeIDStrip, NAoAs = NAoAs, SGridsStartAoA = SGridsStartAoA, SGridsEndAoA = SGridsEndAoA, SGridIDsAoA = SGridIDsAoA, SGridFracsAoA = SGridFracsAoA, NumObsAoA = NumObsAoA, AreaAoA = AreaAoA, TimeIDAoA = TimeIDAoA, NPDists = NPDists, NLines = NLines, SGridsStartLine = SGridsStartLine, SGridsEndLine = SGridsEndLine, SGridIDsLine = SGridIDsLine, SGridFracsLine = SGridFracsLine, LengthLine = LengthLine, TimeIDLine = TimeIDLine, PI = pi, NMaxSGridsAoA = NMaxSGridsAoA, NMaxSGridsStrip = NMaxSGridsStrip, NMaxSGridsLine = NMaxSGridsLine)

NimbleData <- list(X = X, Y = Y, CntStrip = CntStrip, CntAoA = CntAoA, PDists = PDists, CntLine = CntLine)

NimbleInits <- list(sigma_sdi = 0.5, sdi = rep(0, NLGrids), sigma_sl = 0.5, sl = rep(0, NLGrids), beta_di = rep(0, NX), beta_l = rep(0, NY), sigma_di = 0.5, sigma_l = 0.5, di = rep(10, NSGrids), lambda = matrix(1, nrow = NSGrids, ncol = (LastDateID - FirstDateID)), pStrip = 0.5, pAoA = 0.5, sigma_hn = 1)

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

MCMCtrace(object = Samples, pdf = FALSE, ind = TRUE, params = c("beta_di"))


















