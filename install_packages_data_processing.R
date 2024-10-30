# install packages if required
required.packages <- c(
  "terra",
  "sf",
  "rvest",
  "pbapply",
  "httr",
  "stringi",
  "foreach",
  "doParallel",
  "tidyverse",
  "exactextractr",
  "tidyterra",
  "abind",
  "mice",
  "factoextra",
  "nimble",
  "devtools",
  "readr",
  "rstudioapi"
)

new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, quiet = T)
rm(list=ls())
gc()
