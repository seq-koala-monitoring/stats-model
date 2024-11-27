#' Prepare R to model koala popualtion densities across South East Queensland
#' This initial function will download the latest code from GitHub (https://github.com/seq-koala-monitoring/stats-model/tree/main) and configure the working directory to ensure that the analyses run smoothly.

  # install required packages if not already installed
required.packages <- c(
  "ghapps",
  "gh",
  "httr2"
)
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, quiet = T)

  # load packages
  library(ghapps)
  library(gh)
  library(httr2)
  
  # Generate access token
  jwt <- gh_app_jwt(app_id = "980056",
                    app_key = "seq-token.2024-08-26.private-key.pem")
  
  token <- gh_app_token("seq-koala-monitoring", 
                        jwt = jwt)
  
  
  # list stats-model contents
  repo_contents_stats <- gh(
    "GET /repos/:owner/:repo/contents",
    owner = "seq-koala-monitoring",
    repo = "stats-model",
    .token = token
  )
  
  # delete all files from the "code" directory except for initial_setup.R
  file.path <- list.files(path = "code", include.dirs = T, full.names = T)
  file.remove(file.path[!file.path %in% c("code/initial_setup.R",
                                          "code/fcn_update_db.R", 
                                          "code/fcn_set_parameters.R", 
                                          "code/line_transect.R")])
  
  
  # download r scripts
  for (file in repo_contents_stats) {
    if (file$type == "file") {
      file_extension <- tolower(tools::file_ext(file$name))
      if (file_extension == "r") {
        cat("Downloading:", file$name, "\n")
        download.file(
          url = file$download_url,
          destfile = file.path("code", file$name),
          headers = c(Authorization = paste("token", token)),
          mode = "wb"  # Binary mode to ensure proper file handling
        )
      }
    }
  }
  
  rm(list=ls())
  gc()

    # end
  