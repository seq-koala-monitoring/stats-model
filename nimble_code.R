# THIS CONTAINS NIMBLE CODE FOR THE DIFFERENT MODELS

# DATA DEFINITIONS

# NLGrids = number of large grids
# NGPops = number of genetic populations
# NLGridAdjs = number of large grid adjacent pairs
# Adj = adjacency matrix for large grids of length NLGridAdjs
# WeightsAdj = weights for the adjacency matrix for large grids (2km) of length NLGridAdjs
# NumAdj = number of neigbouring locations for each large grid (2km) of length NLGrids
# NSGrids = number of small grids
# LGridID = the large grid ID for each small grid
# GPopID = the genetic population ID for each small grid
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

# MODELS WITH NO VARIABLE SELECTION

# 1. model with independent spatial random effects and no spatial variation in trends

code_sat_spind_const_trend <- nimbleCode({

	# process model

	# spatial random-effects
	# independent spatial random-effects for each genetic population
	for (i in 1:NGPops) {
		sd[i] ~ dnorm(mu_sd, sd = sigma_sd)
	}

	# temporal random-effects
	if (Order == 2) {
		# time-step temporally structured random-effect for densities (based on second order intrisic CAR model)
		# note here weights are taken from Breslow & Clayton (1993) and mean is constrained to be zero
		# annual time step
		td[1:NTime2] ~ dcar_normal(adj = AdjT2[1:NTimeAdjs2], weights = WeightsAdjT2[1:NTimeAdjs2], num = NumAdjT2[1:NTime2], tau = tau_td, c = 2, zero_mean = 1)
	} else {
		# time-step temporally structured random-effect for densities (based on first order intrisic CAR model)
		# note here weights are taken to be 1 and mean is constrained to be zero
		# annual time step
		td[1:NTime] ~ dcar_normal(adj = AdjT[1:NTimeAdjs], weights = WeightsAdjT[1:NTimeAdjs], num = NumAdjT[1:NTime], tau = tau_td, zero_mean = 1)
	}

	# loop through small grids
	for (i in 1:NSGrids) {

		# linear predictor for time-constant predictors
		a[i] <- inprod(beta_d[1:NX], X[i, ])

		# loop through time steps
		for (t in 1:(LastDateID - FirstDateID + 1)) {
			# linear predictor including time-constant predictors, time-varying predictors, spatial random effect
			# and temporal raondom effect
			log(mu_d[i, t]) <- sd[GenPopID[i]] + td[floor((t - 1) / 2) + 1] + a[i] + inprod(beta_d[(NX + 1):(NX + NY)], Y[i, FirstDateID - 1 + t - Lag, 1:NY])

			# density
			d[i, t] ~ dgamma(mean = mu_d[i, t], sd = sigma_d)
		}
	}

	# mean expected density
	for (t in 1:(LastDateID - FirstDateID + 1)) {
		mean_mud[t] <- mean(mu_d[1:NSGrids, t])
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between then the
		# probability of detection is pStrip = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_strip[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Strip[i, 1:NZ])
		# get f0
		f0_strip[i] <- sqrt(2 / (PI * (sigma_hn_strip[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_strip[i] <- cdfhnorm(q = 15, sigma = sigma_hn_strip[i])
		pStrip[i] <- 1 - ((1 - (cdfHN_strip[i] / (f0_strip[i] * 15))) ^ 2)

		# likelihood for counts
		CntStrip[i] ~ dpois(dStrip[i] * AreaStrip[i] * pStrip[i])
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between them the
		# probability of detection is pAoA = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_aoa[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_AoA[i, 1:NZ])
		# get f0
		f0_aoa[i] <- sqrt(2 / (PI * (sigma_hn_aoa[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_aoa[i] <- cdfhnorm(q = 15, sigma = sigma_hn_aoa[i])
		pAoA[i] <- 1 - ((1 - (cdfHN_aoa[i] / (f0_aoa[i] * 15))) ^ 2)

		# likelihood for counts
		CntAoA[i] ~ dpois(dAoA[i] * AreaAoA[i] * pAoA[i])
	}

	# line transects

	# perpendicular distances
	for (i in 1:NPDists) {
		log(sigma_hn[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Line[PDLineIDs[i], 1:NZ])

		# likelihood function for half-normal detection function
		PDists[i] ~ T(dnorm(0, sd = sigma_hn[i]), 0, )
	}

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

		# get the closed form solution for f0
		# valid when detection function is half-normal, detection on the line is certain, distances are not truncated, and detections are not grouped
		# simulate sd for the half-normal detection function for this transect
		log(sigma_hn_line[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Line[i, 1:NZ])
		# get f0
		f0_line[i] <- sqrt(2 / (PI * (sigma_hn_line[i] ^ 2)))

		# likelihood for counts
		# density divided by 10,000 to ensure desities are per hectare rather than per square metre
		CntLine[i] ~ dpois((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000))
	}

	# priors
	mu_sd ~ dnorm(0, sd = 100)
	sigma_sd ~ dunif(0, 10)
	tau_td <- sigma_td ^ -2
	sigma_td ~ dunif(0, 10)
	for (i in 1:(NX + NY)) {
		beta_d[i] ~ dnorm(0, sd = 100)
	}
	sigma_d ~ dunif(0, 10)
	int_shn ~ dnorm(0, sd = 100)
	for (i in 1:NZ) {
		beta_shn[i] ~ dnorm(0, sd = 100)
	}
})

# 2. model with independent spatial random effects and spatial variation in trends

code_sat_spind_var_trend <- nimbleCode({

	# process model

	# spatial random-effects
	# independent spatial random-effects for each genetic population
	for (i in 1:NGPops) {
		sd[i] ~ dnorm(mu_sd, sd = sigma_sd)
	}

	# spatio-temporal random effects
	if (Order == 2) {
		for (i in 1:NGPops){
    		# time-step temporally structured random-effect for densities (based on second order intrisic CAR model)
			# note here weights are taken from Breslow & Clayton (1993) and mean is constrained to be zero
			# annual time step
			std[i, 1:NTime] ~ dcar_normal(adj = AdjT2[1:NTimeAdjs2], weights = WeightsAdjT2[1:NTimeAdjs2], num = NumAdjT2[1:NTime2], tau = tau_std, c = 2, zero_mean = 1)
  		}
	} else {
		for (i in 1:NGPops){
    		# time-step temporally structured random-effect for densities (based on first order intrisic CAR model)
			# note here weights are taken to be 1 and mean is constrained to be zero
			# annual time step
			std[i, 1:NTime] ~ dcar_normal(adj = AdjT[1:NTimeAdjs], weights = WeightsAdjT[1:NTimeAdjs], num = NumAdjT[1:NTime], tau = tau_std, zero_mean = 1)
  		}
	}

	# loop through small grids
	for (i in 1:NSGrids) {

		# linear predictor for time-constant predictors
		a[i] <- inprod(beta_d[1:NX], X[i, ])

		# loop through time steps
		for (t in 1:(LastDateID - FirstDateID + 1)) {
			# linear predictor including time-constant predictors, time-varying predictors, spatial random effect
			# and temporal raondom effect
			log(mu_d[i, t]) <- sd[GenPopID[i]] + std[GenPopID[i], floor((t - 1) / 2) + 1] + a[i] + inprod(beta_d[(NX + 1):(NX + NY)], Y[i, FirstDateID - 1 + t - Lag, 1:NY])

			# density
			d[i, t] ~ dgamma(mean = mu_d[i, t], sd = sigma_d)
		}
	}

	# mean expected density
	for (t in 1:(LastDateID - FirstDateID + 1)) {
		mean_mud[t] <- mean(mu_d[1:NSGrids, t])
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between then the
		# probability of detection is pStrip = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_strip[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Strip[i, 1:NZ])
		# get f0
		f0_strip[i] <- sqrt(2 / (PI * (sigma_hn_strip[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_strip[i] <- cdfhnorm(q = 15, sigma = sigma_hn_strip[i])
		pStrip[i] <- 1 - ((1 - (cdfHN_strip[i] / (f0_strip[i] * 15))) ^ 2)

		# likelihood for counts
		CntStrip[i] ~ dpois(dStrip[i] * AreaStrip[i] * pStrip[i])
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between them the
		# probability of detection is pAoA = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_aoa[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_AoA[i, 1:NZ])
		# get f0
		f0_aoa[i] <- sqrt(2 / (PI * (sigma_hn_aoa[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_aoa[i] <- cdfhnorm(q = 15, sigma = sigma_hn_aoa[i])
		pAoA[i] <- 1 - ((1 - (cdfHN_aoa[i] / (f0_aoa[i] * 15))) ^ 2)

		# likelihood for counts
		CntAoA[i] ~ dpois(dAoA[i] * AreaAoA[i] * pAoA[i])
	}

	# line transects

	# perpendicular distances
	for (i in 1:NPDists) {
		log(sigma_hn[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Line[PDLineIDs[i], 1:NZ])

		# likelihood function for half-normal detection function
		PDists[i] ~ T(dnorm(0, sd = sigma_hn[i]), 0, )
	}

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

		# get the closed form solution for f0
		# valid when detection function is half-normal, detection on the line is certain, distances are not truncated, and detections are not grouped
		# simulate sd for the half-normal detection function for this transect
		log(sigma_hn_line[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Line[i, 1:NZ])
		# get f0
		f0_line[i] <- sqrt(2 / (PI * (sigma_hn_line[i] ^ 2)))

		# likelihood for counts
		# density divided by 10,000 to ensure desities are per hectare rather than per square metre
		CntLine[i] ~ dpois((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000))
	}

	# priors
	mu_sd ~ dnorm(0, sd = 100)
	sigma_sd ~ dunif(0, 10)
	tau_std <- sigma_std ^ -2
	sigma_std ~ dunif(0, 10)
	for (i in 1:(NX + NY)) {
		beta_d[i] ~ dnorm(0, sd = 100)
	}
	sigma_d ~ dunif(0, 10)
	int_shn ~ dnorm(0, sd = 100)
	for (i in 1:NZ) {
		beta_shn[i] ~ dnorm(0, sd = 100)
	}
})

# 3. model with car process spatial random effects and no spatial variation in trends

code_sat_spcar_const_trend <- nimbleCode({

	# process model

	# spatial random-effects
	# large grid spatially structured random-effect for densities (based on first order intrisic CAR model)
	# note here weights are taken to be 1
	sd[1:NLGrids] ~ dcar_normal(adj = AdjS[1:NLGridAdjs], weights = WeightsAdjS[1:NLGridAdjs], num = NumAdjS[1:NLGrids], tau = tau_sd)

	# temporal random-effects
	if (Order == 2) {
		# time-step temporally structured random-effect for densities (based on second order intrisic CAR model)
		# note here weights are taken from Breslow & Clayton (1993) and mean is constrained to be zero
		# annual time step
		td[1:NTime2] ~ dcar_normal(adj = AdjT2[1:NTimeAdjs2], weights = WeightsAdjT2[1:NTimeAdjs2], num = NumAdjT2[1:NTime2], tau = tau_td, c = 2, zero_mean = 1)
	} else {
		# time-step temporally structured random-effect for densities (based on first order intrisic CAR model)
		# note here weights are taken to be 1 and mean is constrained to be zero
		# annual time step
		td[1:NTime] ~ dcar_normal(adj = AdjT[1:NTimeAdjs], weights = WeightsAdjT[1:NTimeAdjs], num = NumAdjT[1:NTime], tau = tau_td, zero_mean = 1)
	}

	# loop through small grids
	for (i in 1:NSGrids) {

		# linear predictor for time-constant predictors
		a[i] <- inprod(beta_d[1:NX], X[i, ])

		# loop through time steps
		for (t in 1:(LastDateID - FirstDateID + 1)) {
			# linear predictor including time-constant predictors, time-varying predictors, spatial random effect
			# and temporal raondom effect
			log(mu_d[i, t]) <- sd[LGridID[i]] + td[floor((t - 1) / 2) + 1] + a[i] + inprod(beta_d[(NX + 1):(NX + NY)], Y[i, FirstDateID - 1 + t - Lag, 1:NY])

			# density
			d[i, t] ~ dgamma(mean = mu_d[i, t], sd = sigma_d)
		}
	}

	# mean expected density
	for (t in 1:(LastDateID - FirstDateID + 1)) {
		mean_mud[t] <- mean(mu_d[1:NSGrids, t])
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between then the
		# probability of detection is pStrip = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_strip[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Strip[i, 1:NZ])
		# get f0
		f0_strip[i] <- sqrt(2 / (PI * (sigma_hn_strip[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_strip[i] <- cdfhnorm(q = 15, sigma = sigma_hn_strip[i])
		pStrip[i] <- 1 - ((1 - (cdfHN_strip[i] / (f0_strip[i] * 15))) ^ 2)

		# likelihood for counts
		CntStrip[i] ~ dpois(dStrip[i] * AreaStrip[i] * pStrip[i])
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between them the
		# probability of detection is pAoA = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_aoa[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_AoA[i, 1:NZ])
		# get f0
		f0_aoa[i] <- sqrt(2 / (PI * (sigma_hn_aoa[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_aoa[i] <- cdfhnorm(q = 15, sigma = sigma_hn_aoa[i])
		pAoA[i] <- 1 - ((1 - (cdfHN_aoa[i] / (f0_aoa[i] * 15))) ^ 2)

		# likelihood for counts
		CntAoA[i] ~ dpois(dAoA[i] * AreaAoA[i] * pAoA[i])
	}

	# line transects

	# perpendicular distances
	for (i in 1:NPDists) {
		log(sigma_hn[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Line[PDLineIDs[i], 1:NZ])

		# likelihood function for half-normal detection function
		PDists[i] ~ T(dnorm(0, sd = sigma_hn[i]), 0, )
	}

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

		# get the closed form solution for f0
		# valid when detection function is half-normal, detection on the line is certain, distances are not truncated, and detections are not grouped
		# simulate sd for the half-normal detection function for this transect
		log(sigma_hn_line[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Line[i, 1:NZ])
		# get f0
		f0_line[i] <- sqrt(2 / (PI * (sigma_hn_line[i] ^ 2)))

		# likelihood for counts
		# density divided by 10,000 to ensure desities are per hectare rather than per square metre
		CntLine[i] ~ dpois((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000))
	}

	# priors
	tau_sd <- sigma_sd ^ -2
	sigma_sd ~ dunif(0, 10)
	tau_td <- sigma_td ^ -2
	sigma_td ~ dunif(0, 10)
	for (i in 1:(NX + NY)) {
		beta_d[i] ~ dnorm(0, sd = 100)
	}
	sigma_d ~ dunif(0, 10)
	int_shn ~ dnorm(0, sd = 100)
	for (i in 1:NZ) {
		beta_shn[i] ~ dnorm(0, sd = 100)
	}
})


# 4. model with car process spatial random effects and spatial variation in trends

code_sat_spcar_var_trend <- nimbleCode({

	# process model

	# spatial random-effects
	# large grid spatially structured random-effect for densities (based on first order intrisic CAR model)
	# note here weights are taken to be 1
	sd[1:NLGrids] ~ dcar_normal(adj = AdjS[1:NLGridAdjs], weights = WeightsAdjS[1:NLGridAdjs], num = NumAdjS[1:NLGrids], tau = tau_sd)

	# spatio-temporal random effects
	if (Order == 2) {
		for (i in 1:NGPops){
    		# time-step temporally structured random-effect for densities (based on second order intrisic CAR model)
			# note here weights are taken from Breslow & Clayton (1993) and mean is constrained to be zero
			# annual time step
			std[i, 1:NTime] ~ dcar_normal(adj = AdjT2[1:NTimeAdjs2], weights = WeightsAdjT2[1:NTimeAdjs2], num = NumAdjT2[1:NTime2], tau = tau_std, c = 2, zero_mean = 1)
  		}
	} else {
		for (i in 1:NGPops){
    		# time-step temporally structured random-effect for densities (based on first order intrisic CAR model)
			# note here weights are taken to be 1 and mean is constrained to be zero
			# annual time step
			std[i, 1:NTime] ~ dcar_normal(adj = AdjT[1:NTimeAdjs], weights = WeightsAdjT[1:NTimeAdjs], num = NumAdjT[1:NTime], tau = tau_std, zero_mean = 1)
  		}
	}

	# loop through small grids
	for (i in 1:NSGrids) {

		# linear predictor for time-constant predictors
		a[i] <- inprod(beta_d[1:NX], X[i, ])

		# loop through time steps
		for (t in 1:(LastDateID - FirstDateID + 1)) {
			# linear predictor including time-constant predictors, time-varying predictors, spatial random effect
			# and temporal raondom effect
			log(mu_d[i, t]) <- sd[LGridID[i]] + std[GenPopID[i], floor((t - 1) / 2) + 1] + a[i] + inprod(beta_d[(NX + 1):(NX + NY)], Y[i, FirstDateID - 1 + t - Lag, 1:NY])

			# density
			d[i, t] ~ dgamma(mean = mu_d[i, t], sd = sigma_d)
		}
	}

	# mean expected density
	for (t in 1:(LastDateID - FirstDateID + 1)) {
		mean_mud[t] <- mean(mu_d[1:NSGrids, t])
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between then the
		# probability of detection is pStrip = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_strip[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Strip[i, 1:NZ])
		# get f0
		f0_strip[i] <- sqrt(2 / (PI * (sigma_hn_strip[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_strip[i] <- cdfhnorm(q = 15, sigma = sigma_hn_strip[i])
		pStrip[i] <- 1 - ((1 - (cdfHN_strip[i] / (f0_strip[i] * 15))) ^ 2)

		# likelihood for counts
		CntStrip[i] ~ dpois(dStrip[i] * AreaStrip[i] * pStrip[i])
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between them the
		# probability of detection is pAoA = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_aoa[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_AoA[i, 1:NZ])
		# get f0
		f0_aoa[i] <- sqrt(2 / (PI * (sigma_hn_aoa[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_aoa[i] <- cdfhnorm(q = 15, sigma = sigma_hn_aoa[i])
		pAoA[i] <- 1 - ((1 - (cdfHN_aoa[i] / (f0_aoa[i] * 15))) ^ 2)

		# likelihood for counts
		CntAoA[i] ~ dpois(dAoA[i] * AreaAoA[i] * pAoA[i])
	}

	# line transects

	# perpendicular distances
	for (i in 1:NPDists) {
		log(sigma_hn[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Line[PDLineIDs[i], 1:NZ])

		# likelihood function for half-normal detection function
		PDists[i] ~ T(dnorm(0, sd = sigma_hn[i]), 0, )
	}

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

		# get the closed form solution for f0
		# valid when detection function is half-normal, detection on the line is certain, distances are not truncated, and detections are not grouped
		# simulate sd for the half-normal detection function for this transect
		log(sigma_hn_line[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Line[i, 1:NZ])
		# get f0
		f0_line[i] <- sqrt(2 / (PI * (sigma_hn_line[i] ^ 2)))

		# likelihood for counts
		# density divided by 10,000 to ensure desities are per hectare rather than per square metre
		CntLine[i] ~ dpois((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000))
	}

	# priors
	tau_sd <- sigma_sd ^ -2
	sigma_sd ~ dunif(0, 10)
	tau_std <- sigma_std ^ -2
	sigma_std ~ dunif(0, 10)
	for (i in 1:(NX + NY)) {
		beta_d[i] ~ dnorm(0, sd = 100)
	}
	sigma_d ~ dunif(0, 10)
	int_shn ~ dnorm(0, sd = 100)
	for (i in 1:NZ) {
		beta_shn[i] ~ dnorm(0, sd = 100)
	}
})

# MODELS WITH VARIABLE SELECTION

# 5. model with independent spatial random effects and no spatial variation in trends

code_sel_spind_const_trend <- nimbleCode({

	# variable selection
	for(i in 1:(NX + NY)) {
    	z[i] ~ dbern(psi) # indicator variable for each coefficient
    	zbeta_d[i] <- z[i] * beta_d[i]  ## indicator * beta
  	}

	# process model

	# spatial random-effects
	# independent spatial random-effects for each genetic population
	for (i in 1:NGPops) {
		sd[i] ~ dnorm(mu_sd, sd = sigma_sd)
	}

	# temporal random-effects
	if (Order == 2) {
		# time-step temporally structured random-effect for densities (based on second order intrisic CAR model)
		# note here weights are taken from Breslow & Clayton (1993) and mean is constrained to be zero
		# annual time step
		td[1:NTime2] ~ dcar_normal(adj = AdjT2[1:NTimeAdjs2], weights = WeightsAdjT2[1:NTimeAdjs2], num = NumAdjT2[1:NTime2], tau = tau_td, c = 2, zero_mean = 1)
	} else {
		# time-step temporally structured random-effect for densities (based on first order intrisic CAR model)
		# note here weights are taken to be 1 and mean is constrained to be zero
		# annual time step
		td[1:NTime] ~ dcar_normal(adj = AdjT[1:NTimeAdjs], weights = WeightsAdjT[1:NTimeAdjs], num = NumAdjT[1:NTime], tau = tau_td, zero_mean = 1)
	}

	# loop through small grids
	for (i in 1:NSGrids) {

		# linear predictor for time-constant predictors
		a[i] <- inprod(zbeta_d[1:NX], X[i, ])

		# loop through time steps
		for (t in 1:(LastDateID - FirstDateID + 1)) {
			# linear predictor including time-constant predictors, time-varying predictors, spatial random effect
			# and temporal raondom effect
			log(mu_d[i, t]) <- sd[GenPopID[i]] + td[floor((t - 1) / 2) + 1] + a[i] + inprod(zbeta_d[(NX + 1):(NX + NY)], Y[i, FirstDateID - 1 + t - Lag, 1:NY])

			# density
			d[i, t] ~ dgamma(mean = mu_d[i, t], sd = sigma_d)
		}
	}

	# mean expected density
	for (t in 1:(LastDateID - FirstDateID + 1)) {
		mean_mud[t] <- mean(mu_d[1:NSGrids, t])
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between then the
		# probability of detection is pStrip = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_strip[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Strip[i, 1:NZ])
		# get f0
		f0_strip[i] <- sqrt(2 / (PI * (sigma_hn_strip[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_strip[i] <- cdfhnorm(q = 15, sigma = sigma_hn_strip[i])
		pStrip[i] <- 1 - ((1 - (cdfHN_strip[i] / (f0_strip[i] * 15))) ^ 2)

		# likelihood for counts
		CntStrip[i] ~ dpois(dStrip[i] * AreaStrip[i] * pStrip[i])
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between them the
		# probability of detection is pAoA = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_aoa[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_AoA[i, 1:NZ])
		# get f0
		f0_aoa[i] <- sqrt(2 / (PI * (sigma_hn_aoa[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_aoa[i] <- cdfhnorm(q = 15, sigma = sigma_hn_aoa[i])
		pAoA[i] <- 1 - ((1 - (cdfHN_aoa[i] / (f0_aoa[i] * 15))) ^ 2)

		# likelihood for counts
		CntAoA[i] ~ dpois(dAoA[i] * AreaAoA[i] * pAoA[i])
	}

	# line transects

	# perpendicular distances
	for (i in 1:NPDists) {
		log(sigma_hn[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Line[PDLineIDs[i], 1:NZ])

		# likelihood function for half-normal detection function
		PDists[i] ~ T(dnorm(0, sd = sigma_hn[i]), 0, )
	}

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

		# get the closed form solution for f0
		# valid when detection function is half-normal, detection on the line is certain, distances are not truncated, and detections are not grouped
		# simulate sd for the half-normal detection function for this transect
		log(sigma_hn_line[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Line[i, 1:NZ])
		# get f0
		f0_line[i] <- sqrt(2 / (PI * (sigma_hn_line[i] ^ 2)))

		# likelihood for counts
		# density divided by 10,000 to ensure desities are per hectare rather than per square metre
		CntLine[i] ~ dpois((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000))
	}

	# priors
	mu_sd ~ dnorm(0, sd = 100)
	sigma_sd ~ dunif(0, 10)
	tau_td <- sigma_td ^ -2
	sigma_td ~ dunif(0, 10)
	for (i in 1:(NX + NY)) {
		beta_d[i] ~ dnorm(0, sd = 100)
	}
	psi ~ dunif(0,1)
	sigma_d ~ dunif(0, 10)
	int_shn ~ dnorm(0, sd = 100)
	for (i in 1:NZ) {
		beta_shn[i] ~ dnorm(0, sd = 100)
	}
})





































# saturated model with spatial car process and no spatial variation in trends

code_sat_spcar_const_trend <- nimbleCode({
	# NLGrids = number of large grids
	# NGPops = number of genetic populations
	# NLGridAdjs = number of large grid adjacent pairs
	# Adj = adjacency matrix for large grids of length NLGridAdjs
	# WeightsAdj = weights for the adjacency matrix for large grids (2km) of length NLGridAdjs
	# NumAdj = number of neigbouring locations for each large grid (2km) of length NLGrids
	# NSGrids = number of small grids
	# LGridID = the large grid ID for each small grid
	# GPopID = the genetic population ID for each small grid
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

	# spatial random-effects
	# large grid spatially structured random-effect for densities (based on first order intrisic CAR model)
	# note here weights are taken to be 1
	sd[1:NLGrids] ~ dcar_normal(adj = AdjS[1:NLGridAdjs], weights = WeightsAdjS[1:NLGridAdjs], num = NumAdjS[1:NLGrids], tau = tau_sd)

	# temporal random-effects
	if (Order == 2) {
		# time-step temporally structured random-effect for densities (based on second order intrisic CAR model)
		# note here weights are taken from Breslow & Clayton (1993) and mean is constrained to be zero
		# annual time step
		td[1:NTime2] ~ dcar_normal(adj = AdjT2[1:NTimeAdjs2], weights = WeightsAdjT2[1:NTimeAdjs2], num = NumAdjT2[1:NTime2], tau = tau_td, c = 2, zero_mean = 1)
	} else {
		# time-step temporally structured random-effect for densities (based on first order intrisic CAR model)
		# note here weights are taken to be 1 and mean is constrained to be zero
		# annual time step
		td[1:NTime] ~ dcar_normal(adj = AdjT[1:NTimeAdjs], weights = WeightsAdjT[1:NTimeAdjs], num = NumAdjT[1:NTime], tau = tau_td, zero_mean = 1)
	}

	# loop through small grids
	for (i in 1:NSGrids) {

		# linear predictor for time-constant predictors
		a[i] <- inprod(beta_d[1:NX], X[i, ])

		# loop through time steps
		for (t in 1:(LastDateID - FirstDateID + 1)) {
			# linear predictor including time-constant predictors, time-varying predictors, spatial random effect
			# and temporal raondom effect
			log(mu_d[i, t]) <- sd[LGridID[i]] + td[floor((t - 1) / 2) + 1] + a[i] + inprod(beta_d[(NX + 1):(NX + NY)], Y[i, FirstDateID - 1 + t - Lag, 1:NY])

			# density
			d[i, t] ~ dgamma(mean = mu_d[i, t], sd = sigma_d)
		}
	}

	# mean spatial random-effect
	mean_sd <- mean(sd[1:NLGrids])

	# mean expected density
	for (t in 1:(LastDateID - FirstDateID + 1)) {
		mean_mud[t] <- mean(mu_d[1:NSGrids, t])
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between then the
		# probability of detection is pStrip = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_strip[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Strip[i, 1:NZ])
		# get f0
		f0_strip[i] <- sqrt(2 / (PI * (sigma_hn_strip[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_strip[i] <- cdfhnorm(q = 15, sigma = sigma_hn_strip[i])
		pStrip[i] <- 1 - ((1 - (cdfHN_strip[i] / (f0_strip[i] * 15))) ^ 2)

		# likelihood for counts
		CntStrip[i] ~ dpois(dStrip[i] * AreaStrip[i] * pStrip[i])
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between them the
		# probability of detection is pAoA = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_aoa[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_AoA[i, 1:NZ])
		# get f0
		f0_aoa[i] <- sqrt(2 / (PI * (sigma_hn_aoa[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_aoa[i] <- cdfhnorm(q = 15, sigma = sigma_hn_aoa[i])
		pAoA[i] <- 1 - ((1 - (cdfHN_aoa[i] / (f0_aoa[i] * 15))) ^ 2)

		# likelihood for counts
		CntAoA[i] ~ dpois(dAoA[i] * AreaAoA[i] * pAoA[i])
	}

	# line transects

	# perpendicular distances
	for (i in 1:NPDists) {
		log(sigma_hn[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Line[PDLineIDs[i], 1:NZ])

		# likelihood function for half-normal detection function
		PDists[i] ~ T(dnorm(0, sd = sigma_hn[i]), 0, )
	}

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

		# get the closed form solution for f0
		# valid when detection function is half-normal, detection on the line is certain, distances are not truncated, and detections are not grouped
		# simulate sd for the half-normal detection function for this transect
		log(sigma_hn_line[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Line[i, 1:NZ])
		# get f0
		f0_line[i] <- sqrt(2 / (PI * (sigma_hn_line[i] ^ 2)))

		# likelihood for counts
		# density divided by 10,000 to ensure desities are per hectare rather than per square metre
		CntLine[i] ~ dpois((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000))
	}

	# priors
	tau_sd <- sigma_sd ^ -2
	sigma_sd ~ dunif(0, 10)
	tau_td <- sigma_td ^ -2
	sigma_td ~ dunif(0, 10)
	for (i in 1:(NX + NY)) {
		beta_d[i] ~ dnorm(0, sd = 100)
	}
	sigma_d ~ dunif(0, 10)
	int_shn ~ dnorm(0, sd = 100)
	for (i in 1:NZ) {
		beta_shn[i] ~ dnorm(0, sd = 100)
	}
})

# saturated model with spatial car process and spatial variation in trends

code_sat_spcar_var_trend <- nimbleCode({
	# NLGrids = number of large grids
	# NGPops = number of genetic populations
	# NLGridAdjs = number of large grid adjacent pairs
	# Adj = adjacency matrix for large grids of length NLGridAdjs
	# WeightsAdj = weights for the adjacency matrix for large grids (2km) of length NLGridAdjs
	# NumAdj = number of neigbouring locations for each large grid (2km) of length NLGrids
	# NSGrids = number of small grids
	# LGridID = the large grid ID for each small grid
	# GPopID = the genetic population ID for each small grid
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

	# spatial random-effects
	# large grid spatially structured random-effect for densities (based on first order intrisic CAR model)
	# note here weights are taken to be 1
	sd[1:NLGrids] ~ dcar_normal(adj = AdjS[1:NLGridAdjs], weights = WeightsAdjS[1:NLGridAdjs], num = NumAdjS[1:NLGrids], tau = tau_sd)

	# spatio-temporal random effects
	if (Order == 2) {
		for (i in 1:NGPops){
    		std[i, 1:NTime] ~ dcar_normal(adj = AdjT2[1:NTimeAdjs2], weights = WeightsAdjT2[1:NTimeAdjs2], num = NumAdjT2[1:NTime2], tau = tau_std, c = 2, zero_mean = 1)
  		}
	} else {
		for (i in 1:NGPops){
    		std[i, 1:NTime] ~ dcar_normal(adj = AdjT[1:NTimeAdjs], weights = WeightsAdjT[1:NTimeAdjs], num = NumAdjT[1:NTime], tau = tau_std, zero_mean = 1)
  		}
	}

	# loop through small grids
	for (i in 1:NSGrids) {

		# linear predictor for time-constant predictors
		a[i] <- inprod(beta_d[1:NX], X[i, ])

		# loop through time steps
		for (t in 1:(LastDateID - FirstDateID + 1)) {
			# linear predictor including time-constant predictors, time-varying predictors, spatial random effect
			# and temporal raondom effect
			log(mu_d[i, t]) <- sd[LGridID[i]] + td[floor((t - 1) / 2) + 1] + std[GenPopID[i], floor((t - 1) / 2) + 1] + a[i] + inprod(beta_d[(NX + 1):(NX + NY)], Y[i, FirstDateID - 1 + t - Lag, 1:NY])

			# density
			d[i, t] ~ dgamma(mean = mu_d[i, t], sd = sigma_d)
		}
	}

	# mean spatial random-effect
	mean_sd <- mean(sd[1:NLGrids])

	# mean expected density
	for (t in 1:(LastDateID - FirstDateID + 1)) {
		mean_mud[t] <- mean(mu_d[1:NSGrids, t])
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between then the
		# probability of detection is pStrip = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_strip[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Strip[i, 1:NZ])
		# get f0
		f0_strip[i] <- sqrt(2 / (PI * (sigma_hn_strip[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_strip[i] <- cdfhnorm(q = 15, sigma = sigma_hn_strip[i])
		pStrip[i] <- 1 - ((1 - (cdfHN_strip[i] / (f0_strip[i] * 15))) ^ 2)

		# likelihood for counts
		CntStrip[i] ~ dpois(dStrip[i] * AreaStrip[i] * pStrip[i])
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between them the
		# probability of detection is pAoA = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_aoa[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_AoA[i, 1:NZ])
		# get f0
		f0_aoa[i] <- sqrt(2 / (PI * (sigma_hn_aoa[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_aoa[i] <- cdfhnorm(q = 15, sigma = sigma_hn_aoa[i])
		pAoA[i] <- 1 - ((1 - (cdfHN_aoa[i] / (f0_aoa[i] * 15))) ^ 2)

		# likelihood for counts
		CntAoA[i] ~ dpois(dAoA[i] * AreaAoA[i] * pAoA[i])
	}

	# line transects

	# perpendicular distances
	for (i in 1:NPDists) {
		log(sigma_hn[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Line[PDLineIDs[i], 1:NZ])

		# likelihood function for half-normal detection function
		PDists[i] ~ T(dnorm(0, sd = sigma_hn[i]), 0, )
	}

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

		# get the closed form solution for f0
		# valid when detection function is half-normal, detection on the line is certain, distances are not truncated, and detections are not grouped
		# simulate sd for the half-normal detection function for this transect
		log(sigma_hn_line[i]) <- int_shn + inprod(beta_shn[1:NZ], Z_Line[i, 1:NZ])
		# get f0
		f0_line[i] <- sqrt(2 / (PI * (sigma_hn_line[i] ^ 2)))

		# likelihood for counts
		# density divided by 10,000 to ensure desities are per hectare rather than per square metre
		CntLine[i] ~ dpois((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000))
	}

	# priors
	tau_sd <- sigma_sd ^ -2
	sigma_sd ~ dunif(0, 10)
	tau_std <- sigma_std ^ -2
	sigma_std ~ dunif(0, 10)
	for (i in 1:(NX + NY)) {
		beta_d[i] ~ dnorm(0, sd = 100)
	}
	sigma_d ~ dunif(0, 10)
	int_shn ~ dnorm(0, sd = 100)
	for (i in 1:NZ) {
		beta_shn[i] ~ dnorm(0, sd = 100)
	}
})
