# install packages if required
required.packages <- c("terra",
                       "sf",
                       "rvest",
                       "pbapply",
                       "httr",
                       "stringi",
                       "foreach",
                       "doParallel",
                       "tidyverse")

new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, quiet = T)
rm(list=ls())
gc()
