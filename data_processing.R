# load libraries
library(tidyverse)
library(nimble)
library(abind)
library(mice)
library(factoextra)
library(terra)

# read utility functions
source("functions.R")

# load input data
Surveys <- readRDS("input/survey_data_500/master.rds")
GridFrac <- readRDS("input/survey_data_500/grid_fractions.rds")
CovConsSurv <- readRDS("input/survey_data_500/cov_constant_array_surveylocations.rds")
CovTempSurv <- readRDS("input/survey_data_500/cov_temporal_array_surveylocations.rds")
AdjL_Queen <- readRDS("input/survey_data_500/adj_data_queen.rds")
DateIntervals <- read_csv("input/survey_data_500/date_interval_lookup.csv") %>% mutate(end_date = as.Date(end_date))
GenPopLookup <- readRDS("input/survey_data_500/gen_pop_lookup.rds")

# set up data for nimble models

# set order for temporal CAR process
Order <- 1
# set the Lag for the effect of predictors on density in terms of number of 6-monthly time steps (can be maximum of 2 currently)
Lag <- 0
# set whether to include CAR spatial model at 5 km grid scale, if not then assumed indepedent random effect based on management unit (0 = no, 1 = yes)
SpCAR <- 0
# set whether ot assume different trend in each management unit (0 = no, 1 = yes)
VarTrend <- 1

# set the date range for analysis
# get first date - either date of first survey or user specified date
#FirstDate <- date("2015-01-01")
FirstDate <- min(c(min(Surveys$line_transect$Date), min(Surveys$strip_transect$Date), min(Surveys$uaoa$Date)))
# get last date - either date of last survey or user specified date
#LastDate <- date("2013-01-01")
LastDate <- max(c(max(Surveys$line_transect$Date), max(Surveys$strip_transect$Date), max(Surveys$uaoa$Date)))
# get first and last date index values
FirstDateID <- filter(DateIntervals, start_date <= FirstDate & end_date >= FirstDate) %>% pull(TimePeriodID)
LastDateID <- filter(DateIntervals, start_date <= LastDate & end_date >= LastDate) %>% pull(TimePeriodID)

# generate array for masks
# load grids
grid.rast <- rast("input/survey_data_500/grid_raster.tif")
grid.vect <- vect("input/survey_data_500/grid_vec.shp")
# load and resample land uses
lu1999.res <- resample(rast("covariates/output/mask/lu1999_mask.tif"), grid.rast, method = "near")
lu2017.res <- resample(rast("covariates/output/mask/lu2017_mask.tif"), grid.rast, method = "near")
# load lot sizes and created a named list
files <- list.files("covariates/output/mask", pattern = "htpls", full.names = T)
names <- list.files("covariates/output/mask", pattern = "htpls", full.names = F)
names <- as.character(parse_number(names))
htpls <- lapply(files, rast)
names(htpls) <- names
# summarize raster values over grid polygons (it takes some minutes to run)
htpls.grid <- lapply(htpls, function(r) {
  # mean lot size per grid cell
  g <- exactextractr::exact_resample(r, grid.rast, "mean") 
  # rename it
  names(g) <- "mean_lot_size"
  # extract mean lot size per GridID
  gid <- g %>% 
    extract(grid.vect, bind = T)})
# add raster date to the previous list
for(i in 1:length(htpls.grid)){
  htpls.grid[[i]]$date <- names(htpls.grid)[i]}
# Load look up tables to mask cells by land use only
lookup1999 <- read_csv("covariates/output/mask/land_use_lookup1999_mask.csv")
lookup2017 <- read_csv("covariates/output/mask/land_use_lookup2017_mask.csv")
# use the median lot size (ha) in Brisbane CBD and inner suburbs (2021) as threshold 
# this value represents highly urbanized areas where koalas are unlikely to occur
median.lot.size.brisbane <- 0.63 
# create temporal mask dataframes (it takes some minutes to run)
htpmask <- lapply(htpls.grid, function(r){
  # set reference dates
  date.1999 <- as.Date("1999-06-01")
  date.2017 <- as.Date("2017-06-01")
  # set lot size date
  date.vect <- as.Date(paste0(str_sub(unique(r$date), end = -3), "-", str_sub(unique(r$date), start = 5), "-", "01"))
  # Use either land use raster from 1999 or 2017 based on the closest date
  if(date.vect <= date.1999 | 
     abs(as.numeric(date.1999 - date.vect)) <= abs(as.numeric(date.2017 - date.vect))){
    j <- extract(lu1999.res, r, bind = T) %>% as.data.frame() %>% left_join(lookup1999, join_by("Tertiary" == "code")) %>% mutate(Mask = case_when(Tertiary.y %in% c("Urban residential", "Public services", "Recreation and culture") & mean_lot_size <= median.lot.size.brisbane ~ 1, mask == 1 ~ 1, .default = 0)) %>% select(GridID, Mask, date)
  } else {
    j <- extract(lu2017.res, r, bind = T) %>% as.data.frame() %>% left_join(lookup2017, join_by("Tertiary" == "code")) %>% mutate(Mask = case_when(Tertiary.y %in% c("Urban residential", "Public services", "Recreation and culture") & mean_lot_size <= median.lot.size.brisbane ~ 1, mask == 1 ~ 1, .default = 0)) %>% select(GridID, Mask, date)
  }})
# assign a TimePeriodID to each temporal dataframe
# store TimePeriodIds
TimePeriod <- vector()
for(i in 1:length(htpmask)){
  name <- names(htpmask)[i]
  name <- as.Date(paste0(str_sub(parse_number(name), end = -3), "-", str_sub(parse_number(name), start = 5), "-", "01"))
  sub <- DateIntervals[which(DateIntervals$start_date <= name & DateIntervals$end_date >= name),]
  htpmask[[i]]$TimePeriodID <- sub$TimePeriodID
  TimePeriod[i] <- as.character(sub$TimePeriodID)}
# name the previous list
names(htpmask) <- TimePeriod
# create an empty list to store a mask dataframe per TimePeriodID
lu.mask.list <- vector(mode = "list", length = length(unique(DateIntervals$TimePeriodID)))
names(lu.mask.list) <- as.character(unique(DateIntervals$TimePeriodID))
# create a vector to represent each TimePeriodID (from 1 to 60) 
index <- as.numeric(names(htpmask))
# append each temporal mask dataframe to the corresponding slot in the empty list
# slots where the TimePeriodID does not have an exact match are filled with data from the closest matching TimePeriodID in the mask dataframe
for(i in 1:length(luMask.list)) {
  n <- which.min(abs(index - i))
  lu.mask.list[[i]] <- htpmask[[n]] %>% select(-TimePeriodID, -date)}
# create an array
lu.mask.array <- abind::abind(lu.mask.list, along = 3)
# remove temporary objects and free up memory
rm(i, n, index, files, names,lookup1999, lookup2017, htpls, lu1999.res, lu2017.res, htpls.grid, sub)
gc()

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
LGridID <- left_join(as.data.frame(CovConsSurv$GridID), AdjL_Queen$grid_lookup, by = c("CovConsSurv$GridID" = "GridID"))$SpGridID
NLGrids <- length(AdjL_Queen$adjacencyList$num)
GenPopID <- left_join(as.data.frame(CovConsSurv$GridID), as.data.frame(GenPopLookup), by = c("CovConsSurv$GridID" = "GridID"))$GENPOP_ID
NGPops <- length(unique(GenPopID))

# static predictors

# set up list for scale parameters
ScaleParamsX <- list()

# ground water dependent ecosystems
# reclassify so that: 1 = intermittent and freshwater, 2 = near-permanent and freshwater, 3 = exclusion and recharge zones, 4 = other
X_hhgde <- as.vector(CovConsSurv[,"hhgde"]) %>% as.character() %>% case_match("1"~4,"2"~1,"3"~2,"4"~4, "5"~4, "6"~3) %>% as.numeric() %>% as.factor()

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
# inspect scree plots
fviz_eig(Soil_PCA)
fviz_eig(Soil_PCA, choice = "eigenvalue")
# inspect bi-plot of PCA components
fviz_pca_var(Soil_PCA, axes = c(1, 2), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
# individual PCA scores for first two components
IndScores <- as_tibble(1:nrow(Soil)) %>% left_join(bind_cols(as_tibble(as.numeric(rownames(get_pca_ind(Soil_PCA)$coord))), as_tibble(get_pca_ind(Soil_PCA)$coord[,1:2])))

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
# reclassfy so that: 1 = remnant high suitability core, 2 = remnant medium suitability core, 3 = non-remnant high suitability core, 4 = other
Y_hhkha <- as.vector(CovTempSurv[,"hhkha",]) %>% as.character() %>% case_match("1"~1,"2"~2,"3"~4,"4"~4, "5"~3, "6"~4, "7"~4, "8"~4) %>% as.numeric() %>% as.factor()

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
write.csv(CorrXY, file="output/collinearity/corr_xy.csv")

# check for collinearity among continuous predictors in X
Corr_X <- cor(X %>% dplyr::select(-hhgde), use = "complete.obs", method = "spearman")
write.csv(Corr_X, file="output/collinearity/cor_x.csv")

# check for collinearity among continuous predictors in Y
Corr_Y <- cor(Y_temp %>% dplyr::select(-htime, -hseas, -hhkha, -htlus), use = "complete.obs", method = "spearman")
write.csv(Corr_Y, file="output/collinearity/cor_y.csv")

# get the design matrix and remove collinear variables for X
# remove elevation and terrain ruggedness index (based on correlations > 0.7)
# ruggedness correlated with slope (r = 0.97)
# elevation correlated with temperature (r = -0.79)
# also removed ground water dependent ecosystems here as this results in very poor mixing in the MCMC runs
# due to poor representation of the categories in the survey data
X <- model.matrix(~ htslo + hspc1 + hspc2 + hcltp + hcltt + hhgde, model.frame(~ htslo + hspc1 + hspc2 + hcltp + hcltt + hhgde, as.data.frame(X), na.action = "na.pass"))

# remove the intercept term
X <- X[, 2:ncol(X)] %>% as.data.frame()

# get number of variables in X
NX <- ncol(X)

# get the design matrix and remove collinear variables for X
# remove persistent green, max temperature, and intensive land-use buffer (based on correlations > 0.7)
# persistent green correlated with persistent green buffer (r = 0.72)
# max temperature correlated with mean temperature (r = 0.91)
# intensive land-use buffer correlated with lot size (r = -0.73)
Y_temp <- model.matrix(~ hhfwc + hhpgr2km + htpls2km + hcpre + hctmn + hseas + hhkha + htlus, model.frame(~ hseas + hhfwc + hhpgr2km + htpls2km + hcpre + hctmn + hseas + hhkha + htlus, as.data.frame(Y_temp), na.action = "na.pass"))

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

# set up data for spatial CAR process
AdjS <- AdjL_Queen$adjacencyList$adj
WeightsAdjS <- AdjL_Queen$adjacencyList$weights
NumAdjS <- AdjL_Queen$adjacencyList$num
NLGridAdjs <- length(AdjL_Queen$adjacencyList$adj)

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

# set up data for temporal CAR process (2nd order with weights acording to ....)
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
NimbleConsts <- list(NLGrids = NLGrids, NLGridAdjs = NLGridAdjs, AdjS = AdjS, WeightsAdjS = WeightsAdjS, NumAdjS = NumAdjS, NGPops = NGPops, NTime = NTime, NTimeAdjs = NTimeAdjs, AdjT = AdjT, WeightsAdjT = WeightsAdjT, NumAdjT = NumAdjT, NTime2 = NTime2, NTimeAdjs2 = NTimeAdjs2, AdjT2 = AdjT2, WeightsAdjT2 = WeightsAdjT2, NumAdjT2 = NumAdjT2, NSGrids = NSGrids, LGridID = LGridID, GenPopID = GenPopID, FirstDateID = FirstDateID, LastDateID = LastDateID, NX = NX, NY = NY, NStrips = NStrips, SGridsStartStrip = SGridsStartStrip, SGridsEndStrip = SGridsEndStrip, SGridIDsStrip = SGridIDsStrip, SGridFracsStrip = SGridFracsStrip, AreaStrip = AreaStrip, TimeIDStrip = TimeIDStrip, NAoAs = NAoAs, SGridsStartAoA = SGridsStartAoA, SGridsEndAoA = SGridsEndAoA, SGridIDsAoA = SGridIDsAoA, SGridFracsAoA = SGridFracsAoA, AreaAoA = AreaAoA, TimeIDAoA = TimeIDAoA, NLines = NLines, SGridsStartLine = SGridsStartLine, SGridsEndLine = SGridsEndLine, SGridIDsLine = SGridIDsLine, SGridFracsLine = SGridFracsLine, LengthLine = LengthLine, TimeIDLine = TimeIDLine, PI = pi, NMaxSGridsAoA = NMaxSGridsAoA, NMaxSGridsStrip = NMaxSGridsStrip, NMaxSGridsLine = NMaxSGridsLine, NZ = NZ, NPDists = NPDists, PDLineIDs = PDLineIDs, Order = Order, Lag = Lag, SpCAR = SpCAR, VarTrend = VarTrend)
NimbleData <- list(X = X, Y = Y, Z_Strip = Z_Strip, Z_AoA = Z_AoA, Z_Line = Z_Line, CntStrip = CntStrip, CntAoA = CntAoA, PDists = PDists, CntLine = CntLine)

# save nimble inputs
saveRDS(list(Constants = NimbleConsts, Data = NimbleData, NamesX = names(X), NamesY = names(Y_temp), ScalingX = ScaleParamsX, ScalingY = ScaleParamsY, FirstDate = FirstDate, LastDate = LastDate, FirstDateID_Orig = FirstDateID_Orig, LastDateID_Orig = LastDateID_Orig, SoilPCA = Soil_PCA), paste0("output/nimble_data/data_order", Order, "_lag", Lag, "_spcar", SpCAR, "_vartrend", VarTrend, "_firstdate", FirstDate, ".rds"))
