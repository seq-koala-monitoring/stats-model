fcn_set_parameters <- function(primary_grid_size = 500,
                               secondary_grid_size = 10,
                               line_transect_buffer = 28.7,
                               cov_impute_buffer = 0,
                               area_buffer = 0,
                               monitoring_units = "Population_boundaries_v2.shp",
                               gen_pop_column_id = "GENPOP_ID",
                               update_database = TRUE) {
  
  gen_pop_file_path <- list.files(
    pattern = monitoring_units,
    recursive = T,
    ignore.case = T,
    full.names = F,
    include.dirs = T
  )
  gen_pop_file_path <- gen_pop_file_path[
    (substring(gen_pop_file_path, nchar(gen_pop_file_path) - 3 + 1) == "shp")]

  if (length(gen_pop_file_path) > 1) {
    stop(
      sprintf(
        "There are multiple files named %s. Please, make sure to keep only the most up-to-date information about survey site habitat in this working directory",
        monitoring.units
      )
    )
  } # close {} of error message
  
  
  param <- 
    list(
      primary_grid_size = primary_grid_size,
      secondary_grid_size = primary_grid_size*secondary_grid_size,
      line_transect_buffer = line_transect_buffer,
      cov_impute_buffer = cov_impute_buffer,
      area_buffer = area_buffer,
      gen_pop_file_path = gen_pop_file_path,
      gen_pop_column_id = gen_pop_column_id,
      update_database = update_database
    )
  
  saveRDS(param, file = "code/parameters_data_processing.rds")
}
fcn_set_parameters()
