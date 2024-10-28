# Run data pipeline for inputs to the statistical model

dates <- fcn_get_date_intervals()

cov_constant_array <- fcn_cov_array('constant', write_path = out_dir)

for (i in 1:length(dates)) {
  d <- dates[[i]]
  id <- d$id
  d <- list(d)
  if (exists("use_parallel") & use_parallel){
    rstudioapi::jobRunScript(
    'run_cov_date.R',
      name = id,
      importEnv = TRUE,
      workingDir = getwd()
    )
    if (i %% 8 == 0) Sys.sleep(180)
    Sys.sleep(10)
  } else {
    source('run_cov_date.R')
  }

}
if (exists("use_parallel") & use_parallel) Sys.sleep(200)
