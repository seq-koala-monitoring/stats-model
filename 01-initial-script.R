# --------------------------------------------------------------------------------
#                          GETTING STARTED  
# --------------------------------------------------------------------------------

# The initial function will download the latest code from GitHub (https://github.com/seq-koala-monitoring/stats-model/tree/main) 
# and configure the working directory to ensure that the analyses run properly. 
#
# You also need to provide APIs - code that allows different software programs to communicate - to download
# a few soil variables from TERN () and GitHub. Please, follow the steps below:
#   1) Navigate to https://geonetwork.tern.org.au/geonetwork/srv/eng/new.account
#   2) Create an account
#
# We designed this project to give users flexibility to adjust parameters while still achieving reliable results.
# If you need to modify any spatial parameter from the default, adjust the arguments in the fcn_set_parameters function.
# The parameters are:
#   - primary_grid_size: A numeric value representing the spatial dimensions (in metres) of each cell in the analysis and mapping grid. Default: 500
#   - secondary_grid_size: A numeric value that multiplies the primary grid size to create a coarser grid. Default: 10
#   - line_transect_buffer: A numeric value (in metres) for the estimated width of line transects. Default: 57.38
#   - cov_impute_buffer: A numeric value (in metres) for the buffer around areas for covariate imputation. Default: 0
#   - area_buffer: A numeric value (in metres) defining the buffer around the study area. Default: 0
#   - monitoring_units: The filename (with extension and within quotation marks) of the spatial layer containing boundaries of monitoring units, such as previous genetic populations. Default: Population_boundaries_v2.shp
#   - gen_pop_column_id: A string (within quotation marks) specifying the column name that identifies each monitoring unit in the spatial layer. Default: "GENPOP_ID"
# For example:
#   Example 1) fcn_set_parameters(primary_grid_size = 250) modifies the primary grid size from 500m to 250m
#   Example 2) fcn_set_parameters(primary_grid_size = 250, monitoring_units = "new_monit_units.shp") modifies the primary grid size from 500m to 250m AND the shapefile containing the boundaries of monitoring units
#
# Lastly, you will be redirected to a new tab with a file named covariate_processing.R. 
# This file will gather and prepare all the information needed for the Bayesian state-space model 
# that estimates densities of koalas across South East Queensland. 
#
# To get started, run the line below by pressing Ctrl + Enter on a Windows PC or Command + Return on a Mac.
{source("code/download_code_github.R")
  
  source("code/fcn_set_parameters.R")
  fcn_set_parameters() # If you want to change any default settings, add the parameter names and their new values inside the brackets. If you're changing more than one, separate them with commas.
  
  file.edit('code/covariate_processing.R')}





