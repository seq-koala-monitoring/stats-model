# We designed this project to give users flexibility to adjust parameters while still achieving reliable results.
# If you need to modify any spatial parameter from the default, adjust the parameters below.

# primary_grid_size: A numeric value representing the spatial dimensions (in metres) of each cell in the analysis and mapping grid. 
primary_grid_size <- 500


# secondary_grid_size: A numeric value that multiplies the primary grid size to create a coarser grid.
secondary_grid_size <- primary_grid_size*10


# line_transect_buffer: A numeric value (in metres) for the estimated width of line transects (i.e., one side from the line transect).
line_transect_buffer <- 28.7


# cov_impute_buffer: A numeric value (in metres) for the buffer around areas for covariate imputation. 
cov_impute_buffer <- 0


# area_buffer: A numeric value (in metres) defining the buffer around the study area.
area_buffer <- 0


# monitoring_units: The filename (with extension and within quotation marks) of the spatial layer containing boundaries of monitoring units, such as previous genetic populations.
monitoring_units <- "Population_boundaries_v2.shp"
gen_pop_file_path <- list.files(
  pattern = monitoring_units,
  recursive = T,
  ignore.case = T,
  full.names = F,
  include.dirs = T
)
gen_pop_file_path <- gen_pop_file_path[
  (substring(gen_pop_file_path, nchar(gen_pop_file_path) - 3 + 1) == "shp")]

if (length(gen_pop_file_path) == 0) {
  stop(sprintf("No % s file found. Please, make sure you have a file named % s in this working directory",
               monitoring_units, monitoring_units))} # close {} of error message
if(length(gen_pop_file_path) > 1){
  message("Multiple files found:")
  for (i in seq_along(gen_pop_file_path)) {
    cat(i, ":", gen_pop_file_path[i], "\n")}
  choice <- as.integer(readline("Enter the number corresponding to the file you want to use: "))
  if (is.na(choice) || choice < 1 || choice > length(gen_pop_file_path)) {
    stop("Invalid selection. Process terminated. Please, rerun this function and choose one of the provided options.")}
  gen_pop_file_path <- gen_pop_file_path[choice]}

# gen_pop_column_id: A string (within quotation marks) specifying the column name that identifies each monitoring unit in the spatial layer.
gen_pop_column_id <- "GENPOP_ID"

# update_database: Should the koala survey database be updated with new data?
# This function integrates new survey data into the modelling database by performing data integrity checks, matching the Site_ID with the spatial locations of survey sites, assigning a monitoring unit (i.e., "genetic" population), estimating missing perpendicular distances, and formatting the data to meet modelling requirements. Optionally, the function can also incorporate the spatial representation of the new surveys into the existing spatial file containing all koala surveys (default is TRUE).
update_database <- TRUE