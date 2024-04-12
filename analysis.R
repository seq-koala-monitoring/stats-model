# load libraries
library(tidyverse)
library(nimble)
library(lubridate)
library(coda)
library(MCMCvis)

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

# adjust code here to set the range for the entire data set or define specific date ranges

# get date of first survey
FirstDate <- date("2015-01-01")
#FirstDate <- date("2013-01-01")
#FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
# get date of last survey
LastDate <- date("2020-12-31")
#LastDate <- date("2013-01-01")
#LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
# first and last date index values
FirstDateID <- filter(DateIntervals, start_date <= FirstDate & end_date >= FirstDate) %>% pull(TimePeriodID)
LastDateID <- filter(DateIntervals, start_date <= LastDate & end_date >= LastDate) %>% pull(TimePeriodID)

# Subset data for the specified date range

# set the Lag for the effect of predictors on growth rate in terms of number of 6-monthly time steps
Lag <- 0

# survey data - selecting only the transects for the time period selected
Surveys$line_transect <- Surveys$line_transect %>% filter(Date >= FirstDate & Date <= LastDate)
Surveys$strip_transect <- Surveys$strip_transect %>% filter(Date >= FirstDate & Date <= LastDate)
Surveys$uaoa <- Surveys$uaoa %>% filter(Date >= FirstDate & Date <= LastDate)
#Surveys$perp_distance <- Surveys$perp_distance %>% filter(Date >= FirstDate & Date <= LastDate)
# grid fractions
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
LGridID <- left_join(as.data.frame(CovConsSurv$GridID), Adj2km_Queen$grid_lookup, by = c("CovConsSurv$GridID" = "GridID"))$'CovConsSurv$GridID'

# predictors for initial density

# persistent green
X_hhpgr <- as.vector(CovTempSurv[,"hhpgr", Lag + 1]) %>% scale() %>% as.vector()
# koala habitat
X_hhkha <- as.vector(CovTempSurv[,"hhkha", Lag + 1]) %>% as.factor()

# collate predictor variables

# compile into a tibble
X <- tibble(hhpgr = X_hhpgr, hhkha = X_hhkha)
# remove data where there are NA values for any covariate
#X <- X %>% filter(!(is.na(hhpgr) | is.na(hhkha)))
# get the design matrix
#X <- model.matrix(~ hhpgr, model.frame(~ hhpgr, as.data.frame(X), na.action = "na.pass"))
X <- model.matrix(~ hhpgr + hhkha, model.frame(~ hhpgr + hhkha, as.data.frame(X), na.action = "na.pass"))
# remove the intercept term
X <- X[, 2:ncol(X)] %>% as.data.frame()

# get number of variables in X
NX <- ncol(X)

# predictors for growth rate (these are all time varying predictors related to habitat, threats, and climate)

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
#Y_temp <- model.matrix(~ hhpgr, model.frame(~ hhpgr, as.data.frame(Y_temp), na.action = "na.pass"))
Y_temp <- model.matrix(~ hhpgr + hhkha, model.frame(~ hhpgr + hhkha, as.data.frame(Y_temp), na.action = "na.pass"))
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
	# PI = pi

	# process model

	# large grid (2km) spatially correlated random-effect for initial densities (based on intrisic CAR model)
	# note here weights are taken to be 1
	#si[1:NLGrids] ~ dcar_normal(adj = Adj[1:NLGridAdjs], weights = WeightsAdj[1:NLGridAdjs], num = NumAdj[1:NLGrids], tau_si)

	# large grid (2km) spatially correlated random-effect for growth rate (based on intrisic CAR model)
	# note here weights are taken to be 1
	#sr[1:NLGrids] ~ dcar_normal(adj = Adj[1:NLGrids], weights = WeightsAdj[1:NLGrids], num = NumAdj[1:NLGrids], tau_sr)

	# loop through small grids
	for (i in 1:NSGrids) {
		# linear predictor for initial density incorporating spatially structured (at the large grid scale) and unstructured (at the small grid scale) stochasticity
		#log(d[i, 1]) ~ dnorm(inprod(beta_i[1:NX], X[i, 1:NX]) #+ si[LGridID[i]], sd = sigma_i)

		mu_di[i] <- exp(int_di + inprod(beta_di[1:NX], X[i, 1:NX]))

		di[i] ~ dgamma(mean = mu_di[i] , sd = sigma_di)

		d[i, 1] <- di[i]

		# intercept for growth rate based on spatially strucutred (at the large grid scale) stochasticity
		#lp_int[i] <- intercept #sr[LGridID[i]]

		# loop through remaining time steps
		for (t in 2:(LastDateID - FirstDateID + 1)) {
			# linear predictor for density including time-dependent variables and unstructured (at the small grid scale) stochasticity
			# note the "- ((sigma_r ^ 2) / 2)" adjustment to ensure E(d) = d for all sigma_r
			# also includes spatially structured (at the large grid scale) stochasticity through lp_int[i]
			mu_l[i, t - 1] <- exp(int_l + inprod(beta_l[1:NY], Y[i, FirstDateID - 1 + t - 1 - Lag, 1:NY]))

			lambda[i, t - 1] ~ dgamma(mean = mu_l[i, t - 1], sd = sigma_l)

			d[i, t] <- d[i, t - 1] * lambda[i, t - 1]
		}

		#for (t in 1:(LastDateID - FirstDateID + 1)) {
			# linear predictor for density including time-dependent variables and unstructured (at the small grid scale) stochasticity
			# note the "- ((sigma_r ^ 2) / 2)" adjustment to ensure E(d) = d for all sigma_r
			# also includes spatially structured (at the large grid scale) stochasticity through lp_int[i]
	#		log(d[i, t]) ~ dnorm(log(d[i, t - 1]) + intercept + inprod(beta_r[1:NY], Y[i, t - 1 - Lag + FirstDateID - 1, 1:NY]) - ((sigma_r ^ 2) / 2), sd = sigma_r)
		#	d[i,t] ~ dgamma(mean = mud, sd = sdd)

		#}
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

		#for (j in SGridsStartStrip[i]:SGridsEndStrip[i]) {
		#	dFracStrip[i, j - SGridsStartStrip[i] + 1] <- d[SGridIDsStrip[j], TimeIDStrip[i] - FirstDateID + 1] * SGridFracsStrip[j]
		#}
		#dStrip[i] <- sum(dFracStrip[i,])
		#dStrip[i] ~ dgamma(1, 1)
		#dStrip[i] ~ dgamma(mean = mudstr, sd = sddstr)

		# get true abundance estimate for the strip transect
		aStrip[i] <- round(dStrip[i] * AreaStrip[i]) #  <- round(dStrip[i] * AreaStrip[i])     #ifelse((round(dStrip[i] * AreaStrip[i]) - (dStrip[i] * AreaStrip[i])) <= 0.5, round(dStrip[i] * AreaStrip[i]), trunc(dStrip[i] * AreaStrip[i]))

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

		#for (j in SGridsStartAoA[i]:SGridsEndAoA[i]) {
		#	dFracAoA[i, j - SGridsStartAoA[i] + 1] <- d[SGridIDsAoA[j], TimeIDAoA[i] - FirstDateID + 1] * SGridFracsAoA[j]
		#}
		#dAoA[i] <- #sum(dFracAoA[i,])
		#dAoA[i] ~ dgamma(1, 1)
		#dAoA[i] ~ dgamma(mean = mudaoa, sd = sddaoa)

		# get true abundance estimate for the strip transect
		aAoA[i] <- round(dAoA[i] * AreaAoA[i]) #<- round(dAoA[i] * AreaAoA[i])  #ifelse((round(dAoA[i] * AreaAoA[i]) - (dAoA[i] * AreaAoA[i])) <= 0.5, round(dAoA[i] * AreaAoA[i]), trunc(dAoA[i] * AreaAoA[i]))

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
		#for (j in 1:NMaxSGridsLine) {
			#if (j <= (SGridsEndLine[i] - SGridsStartLine[i] + 1)) {
		#		dFracLine[i, j] <- d[SGridIDsLine[SGridsStartLine[i] + j - 1], TimeIDLine[i] - FirstDateID + 1] * SGridFracsLine[SGridsStartLine[i] + j - 1] #* max(0, min(1, SGridsEndLine[i] - SGridsStartLine[i] + 1 - j + 1))
			#} else {
			#	dFracLine[i, j] <- 0
			#}
		#}

		for (j in SGridsStartLine[i]:SGridsEndLine[i]) {
			dFracLine[i, j - SGridsStartLine[i] + 1] <- d[SGridIDsLine[j], TimeIDLine[i] - FirstDateID + 1] * SGridFracsLine[j]
		}
		for (j in ((SGridsEndLine[i] - SGridsStartLine[i] + 2):(NMaxSGridsLine + 1))) {
			dFracLine[i, j] <- 0
		}

		dLine[i] <- sum(dFracLine[i, 1:NMaxSGridsLine])

		#dLine[i] ~ dgamma(mean = mudlin, sd = sddlin)

		# likelihood for counts
		# density divided by 10,000 to ensure desities are per hectare rather than per square metre
		CntLine[i] ~ dpois((dLine[i] * 2 * LengthLine[i]) / (f0 * 10000))
	}

	# priors
	#tau_si <- sigma_si ^ -2
	#sigma_si ~ dunif(0, 10)
	#tau_sr <- sigma_sr ^ -2
	#sigma_sr ~ dunif(0, 10)
	#sigma_r ~ dunif(0, 10)
	int_di ~ dnorm(0, sd = 100)
	int_l ~ dnorm(0, sd = 100)
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

NimbleConsts <- list(NLGrids = NLGrids, NLGridAdjs = NLGridAdjs, Adj = Adj, WeightsAdj = WeightsAdj, NumAdj = NumAdj, NSGrids = NSGrids, FirstDateID = FirstDateID, LastDateID = LastDateID, Lag = Lag, NX = NX, NY = NY, NStrips = NStrips, SGridsStartStrip = SGridsStartStrip, SGridsEndStrip = SGridsEndStrip, SGridIDsStrip = SGridIDsStrip, SGridFracsStrip = SGridFracsStrip, NumObsStrip = NumObsStrip, AreaStrip = AreaStrip, TimeIDStrip = TimeIDStrip, NAoAs = NAoAs, SGridsStartAoA = SGridsStartAoA, SGridsEndAoA = SGridsEndAoA, SGridIDsAoA = SGridIDsAoA, SGridFracsAoA = SGridFracsAoA, NumObsAoA = NumObsAoA, AreaAoA = AreaAoA, TimeIDAoA = TimeIDAoA, NPDists = NPDists, NLines = NLines, SGridsStartLine = SGridsStartLine, SGridsEndLine = SGridsEndLine, SGridIDsLine = SGridIDsLine, SGridFracsLine = SGridFracsLine, LengthLine = LengthLine, TimeIDLine = TimeIDLine, PI = pi, NMaxSGridsAoA = NMaxSGridsAoA, NMaxSGridsStrip = NMaxSGridsStrip, NMaxSGridsLine = NMaxSGridsLine)

NimbleData <- list(X = X, Y = Y, CntStrip = CntStrip, CntAoA = CntAoA, PDists = PDists, CntLine = CntLine)

#NimbleInits <- list(pStrip = 0.5, pAoA = 0.5, mudstr = 1, sddstr = 1, mudaoa = 1, sddaoa = 1, mudlin = 1, sddlin = 1, mud = 1, sdd = 1, sigma_hn = 1, dStrip = get.dens.inits(AreaStrip, CntStrip), dAoA = get.dens.inits(AreaAoA, CntAoA), dLine = rep(0.5, NLines), d = matrix(1, nrow = NSGrids, ncol = LastDateID - FirstDateID + 1))

NimbleInits <- list(int_di = 0, int_l = 0, beta_di = rep(0, NX), beta_l = rep(0, NY), sigma_di = 0.5, sigma_l = 0.5, di = rep(10, NSGrids), lambda = matrix(1, nrow = NSGrids, ncol = (LastDateID - FirstDateID)), pStrip = 0.5, pAoA = 0.5, sigma_hn = 1)

#NimbleDims <- list(dFracStrip = c(NStrips, NMaxSGridsStrip), dFracAoA = c(NAoAs, NMaxSGridsAoA), dFracLine = c(NLines, NMaxSGridsLine))

N.iter <- 50000
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

MCMCtrace(object = Samples, pdf = FALSE, ind = TRUE, params = c("beta_l"), iter = 50000)


















