# INTRODUCTION

This code fits spatio-temporal models to the South East Queensland Koala Monitoring Data.

What follows are Step-by-step instructions on how to set-up and run the analysis.

The analysis approach consist of the following ordered steps:

1. Software and folder set-up
2. Data set-up and updates
3. Model set-up and initialisation
4. Predictor variable processing
5. Data processing for analysis
6. Co-linearity checks and predictor variable selection
7. Model fitting
8. Model checks
9. Model outputs and inference

## 1. Software and folder set-up

### Download and install software

Ensure that you have R version 4.4.1 or higher, RStudio, and Rtools installed on your computer. The necessary software for this analysis can be found at:

* R: https://cran.r-project.org/
* RStudio: https://posit.co/downloads/
* Rtools: https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html

### Code and folder set-up

To create the correct folder structure and download the code follow these steps (note that if are just updating from a previous analysis and already have the folder structures in place start from step 2):

1. Create a new folder where you intend to run your analysis (this will be your working directory for the analysis)
2. Download `setup_script.R` from https://github.com/seq-koala-monitoring/stats-model/blob/main/setup_script.R and copy the file to your working directory
3. Open RStudio or R and set the working directory to the folder you copied the script to in step 2
4. Open and run the `setup_script.R` script

## 2. Data set-up and updates

### Koala survey data

STILL TO UPDATE

### Data initialisation

If this is the first time you have set up the analysis then follow the steps below to obtain ansd set up the initialisation data (if you have already previously set up the analysis these steps can be ignored and you can go to the `Predictor variables updates` section):

1. Contact Jonathan Rhodes at <jon.rhodes@qut.edu.au> to request the initialisation data
2. Copy files in the `dcdb` folder into the folder `/input/covariates/raw_data/dcdb`
3. Copy files in the `dem` folder into the folder `/input/covariates/raw_data/dem`
4. Copy files in the `gde` folder into the folder `/input/covariates/raw_data/gde` 
5. Copy files in the `land_use` folder into the folder `/input/covariates/raw_data/land_use`
6. Copy files in the `preclear_koala_habitat` into the folder `/input/covariates/raw_data/preclear_koala_habitat`
7. Copy files in the `remnant_cover` folder into the folder `/input/covariates/raw_data/remnant_cover`
8. Copy files in the `study_area` folder into the folder `/input/covariates/raw_data/study_area`
9. Copy `basedata.gdb` into the folder `/input`
10. Copy files in the `databases` folder into the folder `/input/databses`
11. Copy files in the `group_observers` folder into the folder `/input/group_observers`
12. Copy `KoalaSurveyData.gdb` into the folder `/input`
13. Copy files in the `line_transects_2018` folder into the folder `/input/line_transects_2018`
14. Copy files in the `mask` folder into the folder `/input/mask`


### Predictor variables updates

Most spatio-temporal predictor variables are either static and do not need updating or are updated automatically from online data source locations. However, some will need checking or updating before running the model. The predictor variables that need updating are listed below with update instructions and you shoulf follow these instructions before proceeding.

#### Land-use

Given there is unlikely to be any new land-use data available for South East Queensland in the very near future no updates are required for land-use currently. An update process for land-use when new data becomes available will be developed.

#### Lot size

The cadastre data that is used to calculate the lot soze variable needs to be updaqted once per year and prior to running the analysis via the following process:  

1. Locate the folder `/input/covariates/raw_data/dcdb`
2. Find the `cadastre_superimposed.gdb` geodatabase and open it in a GIS
3. Identify the year of the latest cadastre layer (the layer names all have the format `hlpsz_yyyymm`) and which years are misssing
3. Obtain cadastre layers for the study area for October of each missing year and make sure all layers are in `GDA 1994 MGA Zone 56` projected coordinate system
4. Ensure that each new cadastre layer contains a field called `LOTPLAN` which contains the LOT number and PLAN number concanated together with a space delimeter
5. Name each new cadastre layer using the format `hlpsz_yyyy10` (note all cadaste layers are for October in each year)
6. Add the new cadastre layers to the `cadastre_superimposed.gdb` geodatabase located in the folder `/input/covariates/raw_data/dcdb`

Note: currently the October 2022 cadastre is missing and needs to be be added when the next update is conducted.

#### Koala habitat

Updating koala habitat requires updating a pre-clearing koala habitat layer and updating historical remnant cover layers via the following process:

1. Obtain a shapefile of the new koala habitat preclearing layer (see `/input/covariates/raw_data/preclear_koala_habitat/KoalaPreClear.shp` for format required) and make sure this layer is in `GDA 1994 MGA Zone 56` projected coordinate system
2. The fields contained in this shapefile are not critical apart from the field `HSM_RULE` which must contain the habitat code taken from Figure 6 in "Spatial modelling for koalas in South East Queensland - Report Version 4.0" with the letter "R" in front of the habitat code number)
3. Ensure the preclearing koala habitat layer has complete coverage across the study region
4. Delete the existing `KoalaPreClear.shp` files in the folder `/input/covariates/raw_data/preclear_koala_habitat` and replace with the new `KoalaPreClear.shp` files
5. Obtain shapefiles of the new historical remnant cover layers since 1997 (see files in the folder `/input/covariates/raw_data/remnant_cover` for format required) and make sure all layers are in `GDA 1994 MGA Zone 56` projected coordinate system
6. Ensure filenames of these shapefiles have the format  `QLD_REMVEGCOV_yyyy`
7. There are no specific requirements for fields in these shapefiles
8. Delete the existing remnant cover shapefile files in the folder `/input/covariates/raw_data/remnant_cover` and replace with the new remnant cover shapefile files

#### Forest and sparse woody cover

Since an update of the forest and sparse woody cover data requires software modification email Jonathan Rhodes at <jon.rhodes@qut.edu.au> to request an update.

#### Groundwater dependent ecosystems

The groundwater dependent ecosystems data is a bespoke product provided by the Queensland Department of Environment, Tourism, Science and Innovation. if this ever requires updating email Jonathan Rhodes at <jon.rhodes@qut.edu.au> to discuss.

#### Canopy height

The existing model uses static canopy height data from 2009 as that was the only data available at the time of development. However, there are now time series of canopy height available (e.g., https://portal.tern.org.au/metadata/36c98155-39c8-4eec-9070-a978933f3fa3) and future model developments should consider integrating that data for model improvement. This would require some new software development. 

#### Understory fraction

The existing model uses a static understory fraction data from 2009 as that was the only data available at the time of development. However, if time series become available in the future, model developments should consider integrating that data for model improvement. This would require some new software development.

### TERN API

You will need to provide an API code that allows R to download some datasets from TERN. To do this follow the steps below (these steps can be ignored if you have already set this up for a previous analysis):

1. Navigate to https://portal.tern.org.au/browse/theme
2. Click "Sign in" in the right-hand side
3. Sign in with an account
4. Once logged in, click you name and TERN account
5. Click API Keys > Create API Keys > Type a name > Request API key
6. Open a text editor, create a new text document, and add to the first line of that docuement the text `TERN_API <- "YOUR API CODE"` (where YOUR API CODE is the API key you just generated from the TERN website - it will be a series of letters and numbers)
7. Save the text document to the folder `/keys` and name the file `apis.R`

## 3. Model set-up

We designed this project to give users flexibility to adjust parameters while still achieving reliable results. If you need to modify any spatial parameter from the default, adjust the arguments in the _parameters_init.txt_ file.
The parameters are:
* primary_grid_size: A numeric value representing the spatial dimensions (in metres) of each cell in the analysis and mapping grid.
    `primary_grid_size <- 500`
* secondary_grid_multiple: A numeric value that multiplies the primary grid size to create a coarser secondary grid.
    `secondary_grid_multiple <- 10`
* line_transect_buffer: A numeric value (in metres) for the estimated width of line transects (i.e., one side from the line transect).
    `line_transect_buffer <- 28.7`
* cov_impute_buffer: A numeric value (in metres) for the buffer around areas for covariate imputation.
    `cov_impute_buffer <- 0`
* area_buffer: A numeric value (in metres) defining the buffer around the study area.
    `area_buffer <- 0`
* genetic_units: The filename (with extension and within quotation marks) of the spatial layer containing boundaries of genetic populations.
    `genetic_units <- "Population_boundaries_v2.shp"`
* gen_pop_column_id: A string (within quotation marks) specifying the column name that identifies each genetic unit in the spatial layer.
    `gen_pop_column_id <- "GENPOP_ID"`
* whether to aggregate the genetic population into three regions or not
    `gen_pop_agg <- TRUE`
* update_database: Should the koala survey database be updated with new data?
    * This integrates new survey data into the modelling database by performing data integrity checks, matching the Site_ID with the spatial locations of survey sites, assigning a genetic unit (i.e., "genetic" population), estimating missing perpendicular distances, and formatting the data to meet modelling requirements. Optionally, the function can also incorporate the spatial representation of the new surveys into the existing spatial file containing all koala surveys (default is TRUE).
    `update_database <- TRUE`
* specifies whether to run data processing in parallel (requires RStudio API if true)
    `use_parallel <- TRUE`
* specify whether to run the covariate extraction algorithm; if not, read from disc in the output folder (computationally intensive, 30-minute run if TRUE)
    `run_cov_extraction <- TRUE`
* specify wheher to use spatial data imputation (default = FALSE - unlikely to be needed)
    `use_imputation <- FALSE`
* set a start date for modelling using the format dd/mm/yyyy. If you want to model the entire survey period (i.e., from 1996 to current), leave it as NULL and the model will then start from the date of the first survey.
    `FirstDate <- NULL`
* list of static predictor variables to use,
    * Variables removed by default due to high correlations (> 0.6 or < -0.6) were:
    * elevation,
    * terrain ruggedness index,
    * persistent green,
    * intensive land-use in a 2 km buffer,
    * maximum temperature,
    `static_variables <- c("htslo", "hspc1", "hspc2", "hcltp", "hcltt", "hhgde")`
* list of time dynamic predictor variables to use,
    * Variables removed by default due to high correlations (> 0.6 or < -0.6) were:
    * elevation,
    * terrain ruggedness index,
    * persistent green,
    * intensive land-use in a 2 km buffer,
    * maximum temperature,
    `dynamic_variables <- c("hhfwc", "hhpgr2km", "htpls2km", "hcpre", "hctmn", "hseas", "hhkha", "htlus")`
* list of observer error predictor variables,
    * Variables removed by default due to too much missing data were:
    * weather,
    * cloud cover,
    * wind,
    * canopy cover,
    * subcanopy cover,
    `obs_variables <- c("hhcht", "hhunf", "hhchtunf", "hdtma", "hdpre")`
* When extracting ObserverID as a predictor for observer error, should observers be grouped (from 2021 onwards)?
    * This parameter lets the user choose between two options: (1) Using each observer ID as a variable, which may result in models that are easier to converge but require long computing times; (2) Grouping observers based on their affiliation, as specified in the group_observers.csv file located in input/group_observers/. This option reduces computing time but may lead to convergence issues in the models.
    `obs_groups <- FALSE`
* set whether to mask rainforest,
    `RainMask <- TRUE`

## 4. Predictor variable processing


## 5. Data processing for analysis


## 6. Co-linearity checks and predictor variable selection


## 7. Model fitting


## 8. Model checks


## 9. Model outputs and inference
