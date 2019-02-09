# Reproduce results from Campbell et. al. (2019) "Latent influence networks in global environmental politics."
#
# Users must:
#   1. Specify parent directory
#   2. Make sure all reuqired libraries in code/libraries.R are installed
#
# Warning: some calls may require large amounts of memory and CPU time to complete
#


rm(list=ls())   # clear current environment


#### 
# Specify parent directory
maindir <- "~/User/desktop" # "Specify directory in this string, ~/User/desktop, e.g."
####


#### Source functions
source(file.path(maindir, "code/libraries.R"))
source(file.path(maindir, "code/MLE_functions.R"))
source(file.path(maindir, "code/fit_blin.R"))
source(file.path(maindir, "code/post_functions.R"))
source(file.path(maindir, "code/plot_functions.R"))
####


### 
# Fit blin model using 5 cores. Will create "results" subdirectory and write results there
fit_blin(maindir, resultsfilename="results", ncores=10, write=TRUE, verbose=TRUE)
####


####
# Run counterfactual by looping through 10 simulation sets of size 100; (total in paper is 1e4 simulations)
# recommend running in parallel
# writes to "results" subdirectory
run_counterfactual(maindir, resultsdir="results", writedir="results", nsims=1e3, ncores=10)
####


####
# Make figures from the paper;
# writes plots in "figures" subdirectory
make_figs(maindir, resultsdir="results/final_run_6", plotdir="figures")
####




