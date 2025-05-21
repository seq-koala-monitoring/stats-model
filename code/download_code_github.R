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
  # jwt <- gh_app_jwt(app_id = "980056",
  #                   app_key = "seq-token.2024-08-26.private-key.pem")
  # 
  # token <- gh_app_token("seq-koala-monitoring", 
  #                       jwt = jwt)
  
  
# list stats-model contents
owner <- "seq-koala-monitoring"
repo <- "stats-model"
branch <- "main"

repo_contents_stats <- gh(
  "GET /repos/:owner/:repo/git/trees/:branch?recursive=1",
  owner = owner,
  repo = repo,
  branch = branch
)
  
# delete all files from the "code" directory except for initial_setup.R
file.path <- list.files(path = "code",
                        include.dirs = T,
                        full.names = T)
file.remove(file.path[!file.path %in% c("code/download_code_github.R")])
  
# extract R script paths
r.files <- repo_contents_stats$tree %>%
  purrr::keep( ~ .x$type == "blob" && grepl("\\.R$", .x$path)) %>%
  purrr::map_chr( ~ .x$path)
  
  
# download r scripts
for (file in r.files) {
  cat("Downloading:", file, "\n")
  raw_url <- sprintf("https://raw.githubusercontent.com/%s/%s/%s/%s",
                     owner,
                     repo,
                     branch,
                     file)
  download.file(
    url = raw_url,
    destfile = file.path(file),
    mode = "wb"  # Binary mode to ensure proper file handling
  )
}

# extract parameter's path
par.file <- repo_contents_stats$tree %>%
  purrr::keep( ~ .x$type == "blob" && grepl("parameters_init.txt$", .x$path)) %>%
  purrr::map_chr( ~ .x$path)

raw_url <- sprintf("https://raw.githubusercontent.com/%s/%s/%s/%s",
                   owner,
                   repo,
                   branch, par.file)

download.file(
  url = raw_url,
  destfile = par.file,
  mode = "wb"  # Binary mode to ensure proper file handling
)
  
rm(list = ls())
gc()

# end
  
