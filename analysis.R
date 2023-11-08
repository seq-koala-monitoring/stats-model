# load libraries
library(tidyverse)

# load input data
Surveys <- readRDS("input/survey_data/master.rds")
GridFrac <- readRDS("input/survey_data/grid_fractions.rds")
CovCons <- readRDS("input/survey_data/cov_constant_array.rds")
CovConsSurv <- readRDS("input/survey_data/cov_constant_array_surveylocations.rds")
CovTemp <- readRDS("input/survey_data/cov_temporal_array.rds")
CovTempSurv <- readRDS("input/survey_data/cov_temporal_array_surveylocations.rds")


