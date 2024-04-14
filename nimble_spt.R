
# nimble code for spatio-temporal regression model

NimbleCodeSPT <- nimbleCode({
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

	# large grid spatially structured random-effect for densities (based on first order intrisic CAR model)
	# note here weights are taken to be 1
	sd[1:NLGrids] ~ dcar_normal(adj = AdjS[1:NLGridAdjs], weights = WeightsAdjS[1:NLGridAdjs], num = NumAdjS[1:NLGrids], tau = tau_sd)

	# determine order of temporal CAR process
	if (Order == 1) {
		# time-step temporally structured random-effect for densities (based on first order intrisic CAR model)
		# note here weights are taken to be 1 and mean is constrained to be zero
		td[1:NTime] ~ dcar_normal(adj = AdjT[1:NTimeAdjs], weights = WeightsAdjT[1:NTimeAdjs], num = NumAdjT[1:NTime], tau = tau_td, zero_mean = 1)
	} else {
		# time-step temporally structured random-effect for densities (based on second order intrisic CAR model)
		# note here weights are taken from .... and mean is constrained to be zero
		td[1:NTime2] ~ dcar_normal(adj = AdjT2[1:NTimeAdjs2], weights = WeightsAdjT2[1:NTimeAdjs2], num = NumAdjT2[1:NTime2], tau = tau_td, c = 2, zero_mean = 1)
	}

	# loop through small grids
	for (i in 1:NSGrids) {
		# loop through time steps
		for (t in 1:(LastDateID - FirstDateID + 1)) {
			# predictors for density including time-dependent variables, temporally unstructured (for time step) stochasticity
			# and spatially structured (at the large grid scale) stochasticity
			#mu_d[i, t] <- exp(sd[LGridID[i]] + td[t] + inprod(beta_d[1:NY], Y[i, FirstDateID - 1 + t - Lag, 1:NY]))
			log(mu_d[i, t]) <- sd[LGridID[i]] + td[floor((t / 2) - 0.01) + 1] + inprod(beta_d[1:NY], Y[i, FirstDateID - 1 + t - Lag, 1:NY])

			# density
			d[i, t] ~ dgamma(mean = mu_d[i, t], sd = sigma_d)
		}
	}

	# mean spatial random-effect
	mean_sd <- mean(sd[1:NLGrids])

	# mean density
	for (t in 1:(LastDateID - FirstDateID + 1)) {
		mean_d[t] <- mean(d[1:NSGrids, t])
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
		CntStrip[i] ~ dbin(pStrip, aStrip[i])
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
		CntAoA[i] ~ dbin(pAoA, aAoA[i])
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
	tau_sd <- sigma_sd ^ -2
	sigma_sd ~ dunif(0, 10)
	tau_td <- sigma_td ^ -2
	sigma_td ~ dunif(0, 10)
	for (i in 1:NY) {
		beta_d[i] ~ dnorm(0, sd = 100)
	}
	sigma_d ~ dunif(0, 10)
	pStrip ~ dunif(0,1)
	pAoA ~ dunif(0,1)
	sigma_hn ~ dunif(0, 100)
})
