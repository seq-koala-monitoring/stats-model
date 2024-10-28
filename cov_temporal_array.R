# Stand-alone to write the list of RDS objects into a three-dimensional array

use_large_ram = TRUE

fcn_temporal_covariate_rds_array <- function(out_dir = NULL) {
  out_dir_files <- list.files(out_dir, pattern = "cov_temporal_\\d{6}_\\d{6}.rds")
  pattern <- "(\\d{6}_\\d{6})"
  dates <- regmatches(out_dir_files, regexpr(pattern, out_dir_files))

  path_list <- c()
  for (i in 1:length(dates)) {
    date_id <- dates[i]
    print(paste("Processing covariate date: ", date_id))
    file <- paste0("cov_temporal_", date_id, ".rds")
    if (!(file %in% out_dir_files)) {
      stop(paste("Date", date_id, "does not exist in out_dir."))
    }
    path <- file.path(out_dir, file)
    path_list <- c(path_list, path)
    if (use_large_ram) {
      next
    } else {
      cov <- readr::read_rds(path)
      if (i==1) {
        out <- cov
      } else {
        out <- abind::abind(out, cov, along = 3)
      }
      rm(cov)
      gc()
    }
  }

  if (use_large_ram) {
    obj_list <- lapply(path_list, readr::read_rds)
    return(do.call(abind::abind, c(obj_list, list(along=3))))
  }

  return(out)
}

cov_temporal_array <- fcn_temporal_covariate_rds_array(paste0(out_dir,'/cov_raster'))
readr::write_rds(cov_temporal_array, file.path(out_dir, "cov_temporal_array.rds"))
