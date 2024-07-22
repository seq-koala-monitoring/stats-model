# load libraries
library(tidyverse)
library(abind)
library(nimble)
library(coda)
library(extraDistr)
library(parallel)
library(MCMCvis)
library(terra)
library(tidyterra)
library(foreach)
library(doParallel)

# get functions
source("functions.r")

# set up data frame to store model selection results

