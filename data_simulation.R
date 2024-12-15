#### simul data ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/master_thesis")
source("script_simu_supercomputer.R")
source("source_file_common_all.R")
source("script_CVM.R")
source("script_CVMSE.R")
source("script_lepskichen.R")
source("script_lepskiboot.R")

library(splines)
library(MASS)
library(caret)
library(expm)
library(foreach)
library(doParallel)
library(ggplot2)

rhouv <- 0.5
rhozw <- 0.9
case <- 2
n_values <- 1000
data_param = c(n_values, rhouv, rhozw)

n_MC = 2000

set.seed(47820)

simul_all <- list()

for (n in 1:n_MC){
  simul_inter <- simulate_data_3(data_param, g_sim_3, case)
  simul_all[[n]] <- simul_inter
}

filename <- paste("data_", n_MC, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, ".R" ,sep = "")
save(simul_all, file = filename)

