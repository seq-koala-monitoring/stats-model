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

# MODEL WITH NO VARIABLE SELECTION - SATURATED MODEL (CAN ALSO BE USED TO DO VARIABLE SELECTION VIA RJMCMC)

nimble_sat_model <- nimbleCode({

	# process model
	
	# spatially correlated random-effects for genetic populations
	sd[1:NGPops] ~ dcar_normal(adj = AdjS[1:NGPopsAdjs], weights = WeightsAdjS[1:NGPopsAdjs], num = NumAdjS[1:NGPops], tau = tau_sd)

	# temporal random effects
	if (Order == 2) {
		# time-step temporally structured random-effect for densities (based on second order intrisic CAR model) and the same among genetic populations
		# note here weights are taken from Breslow & Clayton (1993)
		# annual time step

		# if spatially variable trend
		if (VarTrend == 1) {
			for (i in 1:NGPops) {
				std[i, 1:NTime2] ~ dcar_normal(adj = AdjT2[1:NTimeAdjs2], weights = WeightsAdjT2[1:NTimeAdjs2], num = NumAdjT2[1:NTime2], tau = tau_std, c = 2, zero_mean = 1)
			}
		} else {
			td[1:NTime2] ~ dcar_normal(adj = AdjT2[1:NTimeAdjs2], weights = WeightsAdjT2[1:NTimeAdjs2], num = NumAdjT2[1:NTime2], tau = tau_td, c = 2, zero_mean = 1)
		}
	} else {
		# time-step temporally structured random-effect for densities (based on first order intrisic CAR model) and the same among genetic populations
		# note here weights are taken to be 1 and mean is constrained to be zero
		# annual time step

		# if spatially variable trend
		if (VarTrend == 1) {			
			for (i in 1:NGPops) {
			std[i, 1:NTime] ~ dcar_normal(adj = AdjT[1:NTimeAdjs], weights = WeightsAdjT[1:NTimeAdjs], num = NumAdjT[1:NTime], tau = tau_std, zero_mean = 1)

			}
		} else {
			td[1:NTime] ~ dcar_normal(adj = AdjT[1:NTimeAdjs], weights = WeightsAdjT[1:NTimeAdjs], num = NumAdjT[1:NTime], tau = tau_td, zero_mean = 1)
		}
	}

	# loop through small grids
	for (i in 1:NSGrids) {

		# linear predictor for time-constant predictors
		a[i] <- inprod(beta_d[1:NX], X[i, 1:NX])

		# loop through time steps
		for (t in 1:(LastDateID - FirstDateID + 1)) {
			# linear predictor including time-constant predictors, time-varying predictors, spatial random effect
			# and temporal random effect
			if (VarTrend == 1) {
				log(d[i, t]) <- a[i] + sd[GenPopID[i]] + std[GenPopID[i], floor((t - 1) / 2) + 1] + inprod(beta_d[(NX + 1):(NX + NY)], Y[i, FirstDateID - 1 + t - Lag, 1:NY])
			} else {
				log(d[i, t]) <- a[i] + sd[GenPopID[i]] + td[floor((t - 1) / 2) + 1] + inprod(beta_d[(NX + 1):(NX + NY)], Y[i, FirstDateID - 1 + t - Lag, 1:NY])
			}
		}
	}

	# observation models

	# likelihoods

	# observation group random effects
	#for (i in 1:NObsGrps) {
	#	o[i] ~ dnorm(0, sd = sigma_o)
	#}

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
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / (f0 * w) where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between then the
		# probability of detection is pStrip = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_strip[i]) <- inprod(beta_shn[1:NZ], Z_Strip[i, 1:NZ]) #+ o[StripObsGrp[i]]
		# get f0
		f0_strip[i] <- sqrt(2 / (PI * (sigma_hn_strip[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_strip[i] <- cdfhnorm(q = 15, sigma = sigma_hn_strip[i])
		pStrip[i] <- 1 - ((1 - (cdfHN_strip[i] / (f0_strip[i] * 15))) ^ 2)

		# likelihood for counts
		# poisson version
		#CntStrip[i] ~ dpois(dStrip[i] * AreaStrip[i] * pStrip[i])
		# negative binomial version
		# using parameterisation where p = phi / (phi + mu) and r = phi where phi is the dispersion parameter
		# here varianvce = mu + (mu ^ 2) / phi and mean = mu
		CntStrip[i] ~ dnegbin(nbpStrip[i], nbrst)
		nbpStrip[i] <- nbrst / (nbrst + (dStrip[i] * AreaStrip[i] * pStrip[i]))
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
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / (f0 * w) where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between them the
		# probability of detection is pAoA = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_aoa[i]) <- inprod(beta_shn[1:NZ], Z_AoA[i, 1:NZ]) #+ o[AoAObsGrp[i]]
		# get f0
		f0_aoa[i] <- sqrt(2 / (PI * (sigma_hn_aoa[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_aoa[i] <- cdfhnorm(q = 15, sigma = sigma_hn_aoa[i])
		pAoA[i] <- 1 - ((1 - (cdfHN_aoa[i] / (f0_aoa[i] * 15))) ^ 2)

		# likelihood for counts
		# poisson version
		#CntAoA[i] ~ dpois(dAoA[i] * AreaAoA[i] * pAoA[i])
		# negative binomial version
		# using parameterisation where p = phi / (phi + mu) and r = phi where phi is the dispersion parameter
		# here varianvce = mu + (mu ^ 2) / phi and mean = mu
		CntAoA[i] ~ dnegbin(nbpAoA[i], nbraa)
		nbpAoA[i] <- nbraa / (nbraa + (dAoA[i] * AreaAoA[i] * pAoA[i]))
	}

	# line transects

	# perpendicular distances
	for (i in 1:NPDists) {
		log(sigma_hn[i]) <- inprod(beta_shn[1:NZ], Z_Line[PDLineIDs[i], 1:NZ]) #+ o[LineObsGrp[PDLineIDs[i]]]

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

		# simulate sd for the half-normal detection function for this transect
		log(sigma_hn_line[i]) <- inprod(beta_shn[1:NZ], Z_Line[i, 1:NZ]) #+ o[LineObsGrp[i]]
		# get f0
		f0_line[i] <- sqrt(2 / (PI * (sigma_hn_line[i] ^ 2)))

		# likelihood for counts
		# density divided by 10,000 to ensure desities are per hectare rather than per square metre
		# poisson version
		#CntLine[i] ~ dpois((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000))
		# negative binomial version
		# using parameterisation where p = phi / (phi + mu) and r = phi where phi is the dispersion parameter
		# here varianvce = mu + (mu ^ 2) / phi and mean = mu
		CntLine[i] ~ dnegbin(nbpLine[i], nbrli)
		nbpLine[i] <- nbrli / (nbrli + ((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000)))
	}

	# priors
	
	# process model
	tau_sd <- sigma_sd ^ -2
	sigma_sd ~ T(dnorm(0, sd = 0.78), 0, ) # mildly informative prior restricting random effects to roughly within an 30% population density difference between populations with probability ~0.66. A 30% change either side of 1 on the log scale = log(1.3) - log(0.7) = 0.62 so want to set the expected value for sigma to 0.62 (i.e., roughly within 1 sd of the mean for a 66% probability for a normal random effect). Since expected value of half-normal = sigma * srt(2) / sqrt(pi), then chosen expected value of sigma = 0.62*sqrt(pi)/sqrt(2) = 0.78
	if (VarTrend == 1) {
		tau_std <- sigma_std ^ -2
		sigma_std ~ T(dnorm(0, sd = 0.78), 0, ) # mildly informative prior restricting random effects to roughly within an 30% population density change per unit time with probability ~0.66. A 30% change either side of 1 on the log scale = log(1.3) - log(0.7) = 0.62 so want to set the expected value for sigma to 0.62 (i.e., roughly within 1 sd of the mean for a 66% probability for a normal random effect). Since expected value of half-normal = sigma * srt(2) / sqrt(pi), then chosen expected value of sigma = 0.62*sqrt(pi)/sqrt(2) = 0.78
	} else {
		tau_td <- sigma_td ^ -2
		sigma_td ~ T(dnorm(0, sd = 0.78), 0, ) # mildly informative prior restricting random effects to roughly within an 30% population density change per unit time with probability ~0.66. A 30% change either side of 1 on the log scale = log(1.3) - log(0.7) = 0.62 so want to set the expected value for sigma to 0.62 (i.e., roughly within 1 sd of the mean for a 66% probability for a normal random effect). Since expected value of half-normal = sigma * srt(2) / sqrt(pi), then chosen expected value of sigma = 0.62*sqrt(pi)/sqrt(2) = 0.78
	}
 	for (i in 1:(NX + NY)) {
		beta_d[i] ~ dnorm(0, sd = 100)
	}
	
	# observation model
	for (i in 1:NZ) {
		beta_shn[i] ~ dnorm(0, sd = 100)
	}
	#sigma_o ~ T(dnorm(0, sd = 1.38), 0, ) # mildly informative prior restricting random effects to roughly within an 50% difference in detectability from the mean for observer groups with probability ~0.66. A 50% change either side of 1 on the log scale = log(1.5) - log(0.5) = 1.10 so want to set the expected value for sigma to 1.10 (i.e., roughly within 1 sd of the mean for a 66% probability for a normal random effect). Since expected value of half-normal = sigma * srt(2) / sqrt(pi), then chosen expected value of sigma = 1.10*sqrt(pi)/sqrt(2) = 1.38
	nbrst ~ dunif(0, 50)
	nbraa ~ dunif(0, 50)
	nbrli ~ dunif(0, 50)
})

# MODEL FOR DOING POSTERIOR PREDICTIVE CHECKS ON A SATURATED MODEL

nimble_sat_model_check <- nimbleCode({

	# process model
	
	# spatially correlated random-effects for genetic populations
	sd[1:NGPops] ~ dcar_normal(adj = AdjS[1:NGPopsAdjs], weights = WeightsAdjS[1:NGPopsAdjs], num = NumAdjS[1:NGPops], tau = tau_sd)

	# temporal random effects
	if (Order == 2) {
		# time-step temporally structured random-effect for densities (based on second order intrisic CAR model) and the same among genetic populations
		# note here weights are taken from Breslow & Clayton (1993)
		# annual time step

		# if spatially variable trend
		if (VarTrend == 1) {
			for (i in 1:NGPops) {
				std[i, 1:NTime2] ~ dcar_normal(adj = AdjT2[1:NTimeAdjs2], weights = WeightsAdjT2[1:NTimeAdjs2], num = NumAdjT2[1:NTime2], tau = tau_td, c = 2, zero_mean = 1)
			}
		} else {
			td[1:NTime2] ~ dcar_normal(adj = AdjT2[1:NTimeAdjs2], weights = WeightsAdjT2[1:NTimeAdjs2], num = NumAdjT2[1:NTime2], tau = tau_td, c = 2, zero_mean = 1)
		}
	} else {
		# time-step temporally structured random-effect for densities (based on first order intrisic CAR model) and the same among genetic populations
		# note here weights are taken to be 1 and mean is constrained to be zero
		# annual time step

		# if spatially variable trend
		if (VarTrend == 1) {			
			for (i in 1:NGPops) {
			std[i, 1:NTime] ~ dcar_normal(adj = AdjT[1:NTimeAdjs], weights = WeightsAdjT[1:NTimeAdjs], num = NumAdjT[1:NTime], tau = tau_std, zero_mean = 1)

			}
		} else {
			td[1:NTime] ~ dcar_normal(adj = AdjT[1:NTimeAdjs], weights = WeightsAdjT[1:NTimeAdjs], num = NumAdjT[1:NTime], tau = tau_td, zero_mean = 1)
		}
	}

	# loop through small grids
	for (i in 1:NSGrids) {

		# linear predictor for time-constant predictors
		a[i] <- inprod(beta_d[1:NX], X[i, 1:NX])

		# loop through time steps
		for (t in 1:(LastDateID - FirstDateID + 1)) {
			# linear predictor including time-constant predictors, time-varying predictors, spatial random effect
			# and temporal random effect
			if (VarTrend == 1) {
				log(d[i, t]) <- a[i] + sd[GenPopID[i]] + std[GenPopID[i], floor((t - 1) / 2) + 1] + inprod(beta_d[(NX + 1):(NX + NY)], Y[i, FirstDateID - 1 + t - Lag, 1:NY])
			} else {
				log(d[i, t]) <- a[i] + sd[GenPopID[i]] + td[floor((t - 1) / 2) + 1] + inprod(beta_d[(NX + 1):(NX + NY)], Y[i, FirstDateID - 1 + t - Lag, 1:NY])
			}
		}
	}

	# observation models

	# likelihoods

	# observation group random effects
	#for (i in 1:NObsGrps) {
	#	o[i] ~ dnorm(0, sd = sigma_o)
	#}

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
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / (f0 * w) where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between then the
		# probability of detection is pStrip = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_strip[i]) <- inprod(beta_shn[1:NZ], Z_Strip[i, 1:NZ]) #+ o[StripObsGrp[i]]
		# get f0
		f0_strip[i] <- sqrt(2 / (PI * (sigma_hn_strip[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_strip[i] <- cdfhnorm(q = 15, sigma = sigma_hn_strip[i])
		pStrip[i] <- 1 - ((1 - (cdfHN_strip[i] / (f0_strip[i] * 15))) ^ 2)

		# likelihood for counts
		# poisson version
		#CntStrip[i] ~ dpois(dStrip[i] * AreaStrip[i] * pStrip[i])
		# negative binomial version
		# using parameterisation where p = phi / (phi + mu) and r = phi where phi is the dispersion parameter
		# here varianvce = mu + (mu ^ 2) / phi and mean = mu
		CntStrip[i] ~ dnegbin(nbpStrip[i], nbrst)
		nbpStrip[i] <- nbrst / (nbrst + (dStrip[i] * AreaStrip[i] * pStrip[i]))

		# goodness-of-fit
		# compute discrepancy statistics for observed data
		Strip_res_obs[i] <- CntStrip[i] - (dStrip[i] * AreaStrip[i] * pStrip[i])
		# compute discrepancy statistics for simulated data
		CntStrip_sim[i] ~ dnegbin(nbpStrip[i], nbrst)
		Strip_res_sim[i] <- CntStrip_sim[i] - (dStrip[i] * AreaStrip[i] * pStrip[i])
		# get the negative-binomial variance for calculating the chi-squared discrepancy statistics
		Strip_var[i] <- (dStrip[i] * AreaStrip[i] * pStrip[i]) + (((dStrip[i] * AreaStrip[i] * pStrip[i]) ^ 2) / nbrst)
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
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / (f0 * w) where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between them the
		# probability of detection is pAoA = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_aoa[i]) <- inprod(beta_shn[1:NZ], Z_AoA[i, 1:NZ]) #+ o[AoAObsGrp[i]]
		# get f0
		f0_aoa[i] <- sqrt(2 / (PI * (sigma_hn_aoa[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_aoa[i] <- cdfhnorm(q = 15, sigma = sigma_hn_aoa[i])
		pAoA[i] <- 1 - ((1 - (cdfHN_aoa[i] / (f0_aoa[i] * 15))) ^ 2)

		# likelihood for counts
		# poisson version
		#CntAoA[i] ~ dpois(dAoA[i] * AreaAoA[i] * pAoA[i])
		# negative binomial version
		# using parameterisation where p = phi / (phi + mu) and r = phi where phi is the dispersion parameter
		# here varianvce = mu + (mu ^ 2) / phi and mean = mu		
		CntAoA[i] ~ dnegbin(nbpAoA[i], nbraa)
		nbpAoA[i] <- nbraa / (nbraa + (dAoA[i] * AreaAoA[i] * pAoA[i]))

		# goodness-of-fit
		# compute discrepancy statistics for observed data
		AoA_res_obs[i] <- CntAoA[i] - (dAoA[i] * AreaAoA[i] * pAoA[i])
		# compute discrepancy statistics for simulated data
		CntAoA_sim[i] ~ dnegbin(nbpAoA[i], nbraa)
		AoA_res_sim[i] <- CntAoA_sim[i] - (dAoA[i] * AreaAoA[i] * pAoA[i])
		# get the negative-binomial variance for calculating the chi-squared discrepancy statistics
		AoA_var[i] <- (dAoA[i] * AreaAoA[i] * pAoA[i]) + (((dAoA[i] * AreaAoA[i] * pAoA[i]) ^ 2) / nbraa)
	}

	# line transects

	# perpendicular distances
	for (i in 1:NPDists) {
		log(sigma_hn[i]) <- inprod(beta_shn[1:NZ], Z_Line[PDLineIDs[i], 1:NZ]) #+ o[LineObsGrp[PDLineIDs[i]]]

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

		# simulate sd for the half-normal detection function for this transect
		log(sigma_hn_line[i]) <- inprod(beta_shn[1:NZ], Z_Line[i, 1:NZ]) #+ o[LineObsGrp[i]]
		# get f0
		f0_line[i] <- sqrt(2 / (PI * (sigma_hn_line[i] ^ 2)))

		# likelihood for counts
		# density divided by 10,000 to ensure desities are per hectare rather than per square metre
		# poisson version
		#CntLine[i] ~ dpois((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000))
		# negative binomial version
		# using parameterisation where p = phi / (phi + mu) and r = phi where phi is the dispersion parameter
		# here varianvce = mu + (mu ^ 2) / phi and mean = mu
		CntLine[i] ~ dnegbin(nbpLine[i], nbrli)
		nbpLine[i] <- nbrli / (nbrli + ((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000)))

		# goodness-of-fit
		# compute discrepancy statistics for observed data
		Line_res_obs[i] <- CntLine[i] - ((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000))
		# compute discrepancy statistics for simulated data
		CntLine_sim[i] ~ dnegbin(nbpLine[i], nbrli)
		Line_res_sim[i] <- CntLine_sim[i] - ((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000))
		# get the negative-binomial variance for calculating the chi-squared discrepancy statistics
		Line_var[i] <- ((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000)) + ((((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000)) ^ 2) / nbrli)
	}

	# priors
	
	# process model
	tau_sd <- sigma_sd ^ -2
	sigma_sd ~ T(dnorm(0, sd = 0.78), 0, ) # mildly informative prior restricting random effects to roughly within an 30% population density difference between populations with probability ~0.66. A 30% change either side of 1 on the log scale = log(1.3) - log(0.7) = 0.62 so want to set the expected value for sigma to 0.62 (i.e., roughly within 1 sd of the mean for a 66% probability for a normal random effect). Since expected value of half-normal = sigma * srt(2) / sqrt(pi), then chosen expected value of sigma = 0.62*sqrt(pi)/sqrt(2) = 0.78
	if (VarTrend == 1) {
		tau_std <- sigma_std ^ -2
		sigma_std ~ T(dnorm(0, sd = 0.78), 0, ) # mildly informative prior restricting random effects to roughly within an 30% population density change per unit time with probability ~0.66. A 30% change either side of 1 on the log scale = log(1.3) - log(0.7) = 0.62 so want to set the expected value for sigma to 0.62 (i.e., roughly within 1 sd of the mean for a 66% probability for a normal random effect). Since expected value of half-normal = sigma * srt(2) / sqrt(pi), then chosen expected value of sigma = 0.62*sqrt(pi)/sqrt(2) = 0.78
	} else {
		tau_td <- sigma_td ^ -2
		sigma_td ~ T(dnorm(0, sd = 0.78), 0, ) # mildly informative prior restricting random effects to roughly within an 30% population density change per unit time with probability ~0.66. A 30% change either side of 1 on the log scale = log(1.3) - log(0.7) = 0.62 so want to set the expected value for sigma to 0.62 (i.e., roughly within 1 sd of the mean for a 66% probability for a normal random effect). Since expected value of half-normal = sigma * srt(2) / sqrt(pi), then chosen expected value of sigma = 0.62*sqrt(pi)/sqrt(2) = 0.78
	}
 	for (i in 1:(NX + NY)) {
		beta_d[i] ~ dnorm(0, sd = 100)
	}
	
	# observation model
	for (i in 1:NZ) {
		beta_shn[i] ~ dnorm(0, sd = 100)
	}
	#sigma_o ~ T(dnorm(0, sd = 1.38), 0, ) # mildly informative prior restricting random effects to roughly within an 50% difference in detectability from the mean for observer groups with probability ~0.66. A 50% change either side of 1 on the log scale = log(1.5) - log(0.5) = 1.10 so want to set the expected value for sigma to 1.10 (i.e., roughly within 1 sd of the mean for a 66% probability for a normal random effect). Since expected value of half-normal = sigma * srt(2) / sqrt(pi), then chosen expected value of sigma = 1.10*sqrt(pi)/sqrt(2) = 1.38
	nbrst ~ dunif(0, 50)
	nbraa ~ dunif(0, 50)
	nbrli ~ dunif(0, 50)
})

nimble_sat_model_check_old <- nimbleCode({

	# process model
	if (Order == 2) {
		# time-step temporally structured random-effect for densities (based on second order intrisic CAR model) and the same among genetic populations
		# note here weights are taken from Breslow & Clayton (1993)
		# annual time step

		# if spatially variable trend
		if (VarTrend == 1) {
			for (i in 1:NGPops) {
				std[i, 1:NTime2] ~ dcar_normal(adj = AdjT2[1:NTimeAdjs2], weights = WeightsAdjT2[1:NTimeAdjs2], num = NumAdjT2[1:NTime2], tau = tau_td, c = 2)
			}
		} else {
			td[1:NTime2] ~ dcar_normal(adj = AdjT2[1:NTimeAdjs2], weights = WeightsAdjT2[1:NTimeAdjs2], num = NumAdjT2[1:NTime2], tau = tau_td, c = 2)
		}
	} else {
		# time-step temporally structured random-effect for densities (based on first order intrisic CAR model) and the same among genetic populations
		# note here weights are taken to be 1 and mean is constrained to be zero
		# annual time step

		# if spatially variable trend
		if (VarTrend == 1) {
			for (i in 1:NGPops) {
			std[i, 1:NTime] ~ dcar_normal(adj = AdjT[1:NTimeAdjs], weights = WeightsAdjT[1:NTimeAdjs], num = NumAdjT[1:NTime], tau = tau_std)

			}
		} else {
			td[1:NTime] ~ dcar_normal(adj = AdjT[1:NTimeAdjs], weights = WeightsAdjT[1:NTimeAdjs], num = NumAdjT[1:NTime], tau = tau_td)
		}
	}

	# loop through small grids
	for (i in 1:NSGrids) {

		# linear predictor for time-constant predictors
		a[i] <- inprod(beta_d[1:NX], X[i, 1:NX])

		# loop through time steps
		for (t in 1:(LastDateID - FirstDateID + 1)) {
			# linear predictor including time-constant predictors, time-varying predictors, spatial random effect
			# and temporal random effect
			if (VarTrend == 1) {
				log(d[i, t]) <- a[i] + std[GenPopID[i], floor((t - 1) / 2) + 1] + inprod(beta_d[(NX + 1):(NX + NY)], Y[i, FirstDateID - 1 + t - Lag, 1:NY])
			} else {
				log(d[i, t]) <- a[i] + td[floor((t - 1) / 2) + 1] + inprod(beta_d[(NX + 1):(NX + NY)], Y[i, FirstDateID - 1 + t - Lag, 1:NY])
			}
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

		# to calculate detection error assume observers are placed 15m apart and detection follows the same
		# process as for the line transects. In the case the detection probability for each observer
		# is P = integral(0, w, g(x)dx) / w where g(x) is the half-normal detection function
		# and w is the width (in this case 15m). Here integral(0, w, g(x)dx) = integral(0, w, f(x)dx) / f0 where
		# f(x) is the half-normal distribution. Since two observers are observing each 15m gap between then the
		# probability of detection is pStrip = 1 - ((1 - P) ^ 2)
		# first simulate sd for the half-normal detection function for this transect
		log(sigma_hn_strip[i]) <- inprod(beta_shn[1:(NZ - 2)], Z_Strip[i, 1:(NZ - 2)])
		# get f0
		f0_strip[i] <- sqrt(2 / (PI * (sigma_hn_strip[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_strip[i] <- cdfhnorm(q = 15, sigma = sigma_hn_strip[i])
		pStrip[i] <- 1 - ((1 - (cdfHN_strip[i] / (f0_strip[i] * 15))) ^ 2)

		# likelihood for counts
		# poisson version
		#CntStrip[i] ~ dpois(dStrip[i] * AreaStrip[i] * pStrip[i])
		# negative binomial version
		# using parameterisation where p = phi / (phi + mu) and r = phi where phi is the dispersion parameter
		# here varianvce = mu + (mu ^ 2) / phi and mean = mu
		CntStrip[i] ~ dnegbin(nbpStrip[i], nbrst)
		nbpStrip[i] <- nbrst / (nbrst + (dStrip[i] * AreaStrip[i] * pStrip[i]))

		# goodness-of-fit
		# compute discrepancy statistics for observed data
		Strip_res_obs[i] <- CntStrip[i] - (dStrip[i] * AreaStrip[i] * pStrip[i])
		# compute discrepancy statistics for simulated data
		CntStrip_sim[i] ~ dnegbin(nbpStrip[i], nbrst)
		Strip_res_sim[i] <- CntStrip_sim[i] - (dStrip[i] * AreaStrip[i] * pStrip[i])
		# get the negaive-binomial variance for calculating the chi-squared discrepancy statistics
		Strip_var[i] <- (dStrip[i] * AreaStrip[i] * pStrip[i]) + (((dStrip[i] * AreaStrip[i] * pStrip[i]) ^ 2) / nbrst)
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
		log(sigma_hn_aoa[i]) <- inprod(beta_shn[1:(NZ - 2)], Z_AoA[i, 1:(NZ - 2)])
		# get f0
		f0_aoa[i] <- sqrt(2 / (PI * (sigma_hn_aoa[i] ^ 2)))
		# get area under f(x) for half-normal detection function out to 15 m
		cdfHN_aoa[i] <- cdfhnorm(q = 15, sigma = sigma_hn_aoa[i])
		pAoA[i] <- 1 - ((1 - (cdfHN_aoa[i] / (f0_aoa[i] * 15))) ^ 2)

		# likelihood for counts
		# poisson version
		#CntAoA[i] ~ dpois(dAoA[i] * AreaAoA[i] * pAoA[i])
		# negative binomial version
		# using parameterisation where p = phi / (phi + mu) and r = phi where phi is the dispersion parameter
		# here varianvce = mu + (mu ^ 2) / phi and mean = mu
		CntAoA[i] ~ dnegbin(nbpAoA[i], nbraa)
		nbpAoA[i] <- nbraa / (nbraa + (dAoA[i] * AreaAoA[i] * pAoA[i]))

		# goodness-of-fit
		# compute discrepancy statistics for observed data
		AoA_res_obs[i] <- CntAoA[i] - (dAoA[i] * AreaAoA[i] * pAoA[i])
		# compute discrepancy statistics for simulated data
		CntAoA_sim[i] ~ dnegbin(nbpAoA[i], nbraa)
		AoA_res_sim[i] <- CntAoA_sim[i] - (dAoA[i] * AreaAoA[i] * pAoA[i])
		# get the negaive-binomial variance for calculating the chi-squared discrepancy statistics
		AoA_var[i] <- (dAoA[i] * AreaAoA[i] * pAoA[i]) + (((dAoA[i] * AreaAoA[i] * pAoA[i]) ^ 2) / nbraa)		
	}

	# line transects

	# perpendicular distances
	for (i in 1:NPDists) {
		log(sigma_hn[i]) <- inprod(beta_shn[1:(NZ - 2)], Z_Line[PDLineIDs[i], 1:(NZ - 2)])

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

		# simulate sd for the half-normal detection function for this transect
		log(sigma_hn_line[i]) <- inprod(beta_shn[1:(NZ - 2)], Z_Line[i, 1:(NZ - 2)])
		# get f0
		f0_line[i] <- sqrt(2 / (PI * (sigma_hn_line[i] ^ 2)))

		# likelihood for counts
		# density divided by 10,000 to ensure desities are per hectare rather than per square metre
		# poisson version
		#CntLine[i] ~ dpois((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000))
		# negative binomial version
		# using parameterisation where p = phi / (phi + mu) and r = phi where phi is the dispersion parameter
		# here varianvce = mu + (mu ^ 2) / phi and mean = mu
		CntLine[i] ~ dnegbin(nbpLine[i], nbrli)
		nbpLine[i] <- nbrli / (nbrli + ((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000)))

		# goodness-of-fit
		# compute discrepancy statistics for observed data
		Line_res_obs[i] <- CntLine[i] - ((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000))
		# compute discrepancy statistics for simulated data
		CntLine_sim[i] ~ dnegbin(nbpLine[i], nbrli)
		Line_res_sim[i] <- CntLine_sim[i] - ((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000))
		# get the negaive-binomial variance for calculating the chi-squared discrepancy statistics
		Line_var[i] <- ((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000)) + ((((dLine[i] * 2 * LengthLine[i]) / (f0_line[i] * 10000)) ^ 2) / nbrli)		
	}

	# priors
	if (VarTrend == 1) {
		tau_std <- sigma_std ^ -2
		sigma_std ~ T(dnorm(0, sd = 0.75), 0, ) # mildly informative prior restricting random effects to roughly within an 80% population density change with probability ~0.66. Chosen sd = 0.6*sqrt(pi)/sqrt(2)
	} else {
		tau_td <- sigma_td ^ -2
		sigma_td ~ T(dnorm(0, sd = 0.75), 0, ) # mildly informative prior restricting random effects to roughly within an 80% population density change with probability ~0.66. Chosen sd = 0.6*sqrt(pi)/sqrt(2)
	}
 	for (i in 1:(NX + NY)) {
		beta_d[i] ~ dnorm(0, sd = 100)
	}
	for (i in 1:(NZ - 2)) {
		beta_shn[i] ~ dnorm(0, sd = 100)
	}
	nbrst ~ dunif(0, 50)
	nbraa ~ dunif(0, 50)
	nbrli ~ dunif(0, 50)
})
