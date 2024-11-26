# Run data pipeline for inputs to the statistical model

dates <- fcn_get_date_intervals()

cov_constant_array <- fcn_cov_array('constant', write_path = out_dir)

if(exists("use_parallel") & use_parallel){
  library(parallel)
  n.cores <- detectCores() - 1
  cl <- makeCluster(n.cores)
  clusterEvalQ(cl, {
    library("SEQKoalaDataPipeline")
    library("dplyr")
  })
  clusterExport(cl, c("fcn_cov_array_detect",
                      "dates",
                      "fcn_covariate_layer_df_detect",
                      "fcn_covariate_interval_mean_detect",
                      "fcn_extract_cov_date_detect",
                      "fcn_cov_array_detect"))
  
  parLapply(cl, dates, function(d){
    fcn_cov_array_detect('temporal', d, write_path = paste0(out_dir, "/cov_raster"))
    })
  stopCluster(cl)
} else {
  source('code/run_cov_date.R')
}
