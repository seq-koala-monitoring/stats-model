# load libraries
library(tidyverse)
library(nimble)
library(lubridate)

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

# get date of first survey
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
# get date of last survey
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
# first and last date index values
FirstDateID <- filter(DateIntervals, start_date <= FirstDate & end_date >= FirstDate) %>% pull(TimePeriodID)
LastDateID <- filter(DateIntervals, start_date <= LastDate & end_date >= LastDate) %>% pull(TimePeriodID)

# data for large grids for spatial autocorrelation
NLGrids <- length(Adj2km_Queen$adjacencyList$num)
NLGridAdjs <- length(Adj2km_Queen$adjacencyList$adj)
Adj <- Adj2km_Queen$adjacencyList$adj
WeightsAdj <- Adj2km_Queen$adjacencyList$weights
NumAdj <- Adj2km_Queen$adjacencyList$num

# data for small grids (only for those within the surveyed areas)
NSGrids <- length(CovConsSurv$GridID)
LGridID <- left_join(as.data.frame(CovConsSurv$GridID), Adj2km_Queen$grid_lookup, by = c("CovConsSurv$GridID" = "GridID"))$'CovConsSurv$GridID'

# predictors for initial density

# persistent green
X_hhpgr <- as.vector(CovTempSurv[,"hhpgr", FirstDateID]) %>% scale() %>% as.vector()
# koala habitat
X_hhkha <- as.vector(CovTempSurv[,"hhkha",FirstDateID]) %>% as.factor()

# collate predictor variables

# compile into a tibble
X <- tibble(hhpgr = X_hhpgr, hhkha = X_hhkha)
# remove data where there are NA values for any covariate
#X <- X %>% filter(!(is.na(hhpgr) | is.na(hhkha)))
# get the design matrix
X <- model.matrix(~ hhpgr + hhkha, model.frame(~ hhpgr + hhkha, as.data.frame(X), na.action = "na.pass"))
# remove the intercept term
X <- X[, 2:ncol(X)]

# get number of variables in X
NX <- ncol(X)

# predictors for growth rate (these are all time varying predictors related to habitat, threats, and climate)

# set the Lag for the effect of predictors on growth rate in terms of number of 6-monthly time steps
Lag <- 0

# persistent green
Y_hhpgr <- as.vector(CovTempSurv[,"hhpgr",]) %>% scale() %>% as.vector()
# koala habitat
Y_hhkha <- as.vector(CovTempSurv[,"hhkha",]) %>% as.factor()

# collate predictor variables

# compile into a tibble
Y_temp <- tibble(hhpgr = Y_hhpgr, hhkha = Y_hhkha)
# remove data where there are NA values for any covariate
#Y_temp <- Y_temp %>% filter(!(is.na(hhpgr) | is.na(hhkha)))
# get the design matrix
Y_temp <- model.matrix(~ hhpgr + hhkha, model.frame(~ hhpgr + hhkha, as.data.frame(Y_temp), na.action = "na.pass"))
# remove the intercept term
Y_temp <- Y_temp[, 2:ncol(Y_temp)]

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

# get the time step ID for each strip transect
TimeIDStrip <- StripJoinGroupFrac$TimePeriodID %>% as.vector()

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

# get the koala count for each all of areas search
CntAoA <- AoAJoinGroupFrac$Number_Sightings %>% as.vector()

# line transect data

# get the number of perependicular distances
NPDists <- nrow(Surveys$perp_distance) %>% as.vector()

# get perpendicular distances
PDists <- Surveys$perp_distance$Perp_Dist %>% as.vector()

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

# get the koala count for each line transect
CntLine <- LineJoinGroupFrac$Number_Sightings %>% as.vector()

# nimble model

NimbleCode <- nimbleCode({
	# NLGrids = number of large grids (2km)
	# NLGridAdjs = number of large grid (2km) adjacent pairs
	# Adj = adjacency matrix for large grids (2km) of length NLGridAdjs
	# WeightsAdj = weights for the adjacency matrix for large grids (2km) of length NLGridAdjs
	# NumAdj = number of neigbouring locations for each large grid (2km) of length NLGrids
	# NSGrids = number of small grids (100m)
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

	# process model

	# large grid (2km) spatially correlated random-effect for initial densities (based on intrisic CAR model)
	# note here weights are taken to be 1
	#si[1:NLGrids] ~ dcar_normal(adj = Adj[1:NLGridAdjs], weights = WeightsAdj[1:NLGridAdjs], num = NumAdj[1:NLGrids], tau_si)

	# large grid (2km) spatially correlated random-effect for growth rate (based on intrisic CAR model)
	# note here weights are taken to be 1
	#sr[1:NLGrids] ~ dcar_normal(adj = Adj[1:NLGrids], weights = WeightsAdj[1:NLGrids], num = NumAdj[1:NLGrids], tau_sr)

	# loop through small grids (100m)
	#for (i in 1:NSGrids) {
		# linear predictor for initial density incorporating spatially structured (at the large grid scale) and unstructured (at the small grid scale) stochasticity
		#log(d[i, 1]) ~ dnorm(inprod(beta_i[1:NX], X[i, 1:NX]) + si[LGridID[i]], sd = sigma_i)

	#	log(d[i, 1]) ~ dnorm(inprod(beta_i[1:NX], X[i, 1:NX]) + start_int, sd = sigma_i)

		# intercept for growth rate based on spatially strucutred (at the large grid scale) stochasticity
	#	lp_int[i] <- intercept #sr[LGridID[i]]

		# loop through remaining time steps
	#	for (t in 2:(LastDateID - FirstDateID + 1)) {
			# linear predictor for density including time-dependent variables and unstructured (at the small grid scale) stochasticity
			# note the "- ((sigma_r ^ 2) / 2)" adjustment to ensure E(d) = d for all sigma_r
			# also includes spatially structured (at the large grid scale) stochasticity through lp_int[i]
	#		log(d[i, t]) ~ dnorm(log(d[i, t - 1]) + intercept + inprod(beta_r[1:NY], Y[i, t - 1 - Lag + FirstDateID - 1, 1:NY]) - ((sigma_r ^ 2) / 2), sd = sigma_r)
	#	}
	#}

	# observation models

	# likelihoods

	# strip transects
	for (i in 1:NStrips) {

		# get true density estimate for the strip transect
		#for (j in SGridsStartStrip[i]:SGridsEndStrip[i]) {
		#	dFracStrip[i, j - SGridsStartStrip[i] + 1] <- d[SGridIDsStrip[j], TimeIDStrip[i] - FirstDateID + 1] * SGridFracsStrip[j]
		#}
		#dStrip[i] <- sum(dFracStrip[i,])
		dStrip[i] ~ dgamma(1, 1)

		# get true abundance estimate for the strip transect
		aStrip[i] <- round(dStrip[i] * AreaStrip[i])     #ifelse((round(dStrip[i] * AreaStrip[i]) - (dStrip[i] * AreaStrip[i])) <= 0.5, round(dStrip[i] * AreaStrip[i]), trunc(dStrip[i] * AreaStrip[i]))

		# likelihood for counts
		CntStrip[i] ~ dbin(1 - ((1 - pStrip) ^ NumObsStrip[i]), aStrip[i])
	}

	# all of area searches
	for (i in 1:NAoAs) {

		# get true density estimate for the all of area search
		#for (j in SGridsStartAoA[i]:SGridsEndAoA[i]) {
		#	dFracAoA[i, j - SGridsStartAoA[i] + 1] <- d[SGridIDsAoA[j], TimeIDAoA[i] - FirstDateID + 1] * SGridFracsAoA[j]
		#}
		dAoA[i] <- 100#sum(dFracAoA[i,])

		# get true abundance estimate for the strip transect
		aAoA[i] <- round(dAoA[i] * AreaAoA[i])  #ifelse((round(dAoA[i] * AreaAoA[i]) - (dAoA[i] * AreaAoA[i])) <= 0.5, round(dAoA[i] * AreaAoA[i]), trunc(dAoA[i] * AreaAoA[i]))

		# likelihood for counts
		CntAoA[i] ~ dbin(1 - ((1 - pAoA) ^ NumObsAoA[i]), aAoA[i])
	}

	# line transects

	# perpendicular distances
	for (i in 1:NPDists) {
		# likelihood function for half-normal detection function
		PDist[i] ~ T(dnorm(0, sd = sigma_hn), 0, )
	}

	# get the closed form solution for f0
	# valid when detection function is half-normal, detection of the line is certain, not truncated, and detections not grouped
	f0 <- sqrt(2 / (pi * (sigma_hn ^ 2)))

	# loop through line transects
	for (i in 1:NLines)
	{
		# get true density estimate for the line transect
		#for (j in SGridsStartLine[i]:SGridsEndLine[i]) {
		#	dFracLine[i, j - SGridsStartLine[i] + 1] <- d[SGridIDsLine[j], TimeIDLine[i] - FirstDateID + 1] * SGridFracsLine[j]
		#}
		dLine[i] <- 100#sum(dFracLine[i,])

		# likelihood for counts
		# density divided by 10,000 to ensure desities are per hectare rather than per square metre
		CntLine[i] ~ dpois((dLine[i] * 2 * LengthLine[i]) / (f0 * 10000))
	}

	# priors
	#tau_si <- sigma_si ^ -2
	#sigma_si ~ dunif(0, 10)
	#tau_sr <- sigma_sr ^ -2
	#sigma_sr ~ dunif(0, 10)
	sigma_r ~ dunif(0, 10)
	pStrip ~ dunif(0,1)
	pAoA ~ dunif(0,1)
	sigma_hn ~ dunif(0, 10)
	beta_i[1:NX] ~ dnorm(0,0.001)
	beta_r[1:NY] ~ dnorm(0,0.001)
	start_int ~ dnorm(0,0.001)
	intercept ~ dnorm(0,0.001)
})

NimbleConsts <- list(NLGrids = NLGrids, NLGridAdjs = NLGridAdjs, Adj = Adj, WeightsAdj = WeightsAdj, NumAdj = NumAdj, NSGrids = NSGrids, FirstDateID = FirstDateID, LastDateID = LastDateID, Lag = Lag, NX = NX, NY = NY, NStrips = NStrips, SGridsStartStrip = SGridsStartStrip, SGridsEndStrip = SGridsEndStrip, SGridIDsStrip = SGridIDsStrip, SGridFracsStrip = SGridFracsStrip, NumObsStrip = NumObsStrip, AreaStrip = AreaStrip, TimeIDStrip = TimeIDStrip, NAoAs = NAoAs, SGridsStartAoA = SGridsStartAoA, SGridsEndAoA = SGridsEndAoA, SGridIDsAoA = SGridIDsAoA, SGridFracsAoA = SGridFracsAoA, NumObsAoA = NumObsAoA, AreaAoA = AreaAoA, TimeIDAoA = TimeIDAoA, NPDists = NPDists, NLines = NLines, SGridsStartLine = SGridsStartLine, SGridsEndLine = SGridsEndLine, SGridIDsLine = SGridIDsLine, SGridFracsLine = SGridFracsLine, LengthLine = LengthLine, TimeIDLine = TimeIDLine)

NimbleData <- list(X = X, Y = Y, CntStrip = CntStrip, CntAoA = CntAoA, PDists = PDists, CntLine = CntLine)

NimbleInits <- list(sigma_si = 1, sigma_sr = 1, sigma_r = 1, pStrip = 0.5, pAoA = 0.5, sigma_hn = 1, beta_i = rep(0, NX), beta_r = rep(0, NY))

NimbleDims <- list(dFracStrip = c(NStrips, NMaxSGridsStrip), dFracAoA = c(NAoAs, NMaxSGridsAoA), dFracLine = c(NLines, NMaxSGridsLine))

NimbleModel <- nimbleModel(code = NimbleCode, name = "koala1", constants = NimbleConsts, dimensions = NimbleDims, data = NimbleData, inits = NimbleInits)


pumpConsts <- list(N = 10,
                   t = c(94.3, 15.7, 62.9, 126, 5.24,
                       31.4, 1.05, 1.05, 2.1, 10.5))

pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))

pumpInits <- list(alpha = 1, beta = 1,
                  theta = rep(0.1, pumpConsts$N))






