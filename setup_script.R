# --------------------------------------------------------------------------------
#                          GETTING STARTED  
# --------------------------------------------------------------------------------

# The initial function will download the latest code from GitHub (https://github.com/seq-koala-monitoring/stats-model/tree/main) 
# and configure the working directory to ensure that the analyses run properly. 
#
# You also need to provide APIs - code that allows different software programs to communicate - to download
# a few soil variables from TERN. Please, follow the steps below:
#   1) Navigate to https://geonetwork.tern.org.au/geonetwork/srv/eng/new.account
#   2) Create an account
#
# After downloading the code, open the next script to continue with the analysis. Go to File > Open file... > locate the file named covariate_processing.R. > Open.
# This next script will gather and prepare all the information needed for the Bayesian state-space model 
# that estimates densities of koalas across South East Queensland. 
#
# To get started, run the line below by pressing Ctrl + Enter on a Windows PC or Command + Return on a Mac.
source("code/download_code_github.R")

# You can now close this script by clicking the X next to its name in the script tab.



