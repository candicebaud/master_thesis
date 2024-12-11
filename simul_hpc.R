#setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/master_thesis")
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


#### Only one parameter framework ####

source("/softs/R/createCluster.R")
cl <- createCluster()
registerDoParallel(cl)
rhouv <- 0.5
rhozw <- 0.9
case <- 2
n_values <- 1000
data_param = c(n_values, rhouv, rhozw)

degree = 3
x_evaluation = seq(-2, 2, length.out = 100)
n_MC = 2000
J_bs <- c(5, 7, 11, 19, 35)
J_ns <- c(3, 5, 9, 17, 33)
p_train = 0.5
n_boot = 100


set.seed(47820)

filename <- paste("data_", n_MC, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, ".R" ,sep = "")
load(filename)

# Parallelize the MC simulations
results <-foreach(n = 1:n_MC, .packages = c("splines", "MASS", "caret", "expm")) %dopar% {
  tryCatch({
  # Simulate data for each iteration
  #simul <- simulate_data_3(data_param, g_sim_3, case)
  simul <- simul_all[[n]]
  
  # Perform the J-optimality computation
  res <- J_opt_data_fixed(simul, J_bs, J_ns, p_train, x_evaluation,
                          degree, n_boot)
  filen <- paste("opt_", n_MC, "_degree", degree, "_ptrain", p_train, "_nboot", n_boot, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, "_simu", n, ".R" ,sep = "")
  save(res, file = filen)
  
  # Store the results from the J-optimality computations
  list(
    J_opt_CV_M_bs = res$J_opt_CV_M_bs,
    J_opt_CV_M_ns = res$J_opt_CV_M_ns,
    J_opt_CVMSE_bs = res$J_opt_CVMSE_bs,
    J_opt_CVMSE_ns = res$J_opt_CVMSE_ns,
    J_opt_lepski_bs = res$J_opt_lepski_bs,
    J_opt_lepski_ns = res$J_opt_lepski_ns,
    J_opt_lepskiboot_bs = res$J_opt_lepskiboot_bs,
    J_opt_lepskiboot_ns = res$J_opt_lepskiboot_ns,
    gamma_bs = res$list_gamma_bs,
    gamma_ns = res$list_gamma_ns,
    g_hat_J_bs = res$list_g_hat_J_bs,
    g_hat_J_ns = res$list_g_hat_J_ns
  ) }, error = function(e) {
    # In case of error, return NULL or an empty list, so the iteration is skipped
    return(NULL)
  })
}
filename = paste("opt_", n_MC , "_degree", degree, "_ptrain", p_train, "_nboot", n_boot, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, ".R" ,sep = "")
save(results, file = filename)


stopCluster(cl)








#source("/softs/R/createCluster.R")
#cl <- createCluster()
#registerDoParallel(cl)
#N_values <- c(200) #pour commencer
#Cases <- c(2, 3)
#Rhouv_Rhozw <- list(c(0.5, 0.9), c(0.8, 0.9), c(0.8, 0.7))

#set.seed(47820)

#parameter_combinations <- expand.grid(
#  N = N_values,
#  Case = Cases,
#  Rhouv_Rhozw = seq_along(Rhouv_Rhozw) # Use indices for combinations
#)

#degree = 3
#x_evaluation = seq(-2, 2, length.out = 100)
#n_MC = 50
#J_bs <- c(5, 7, 11, 19, 35)
#J_ns <- c(3, 5, 9, 17, 33)
#p_train = 0.8
#n_boot = 100

#foreach (j=1:nrow(parameter_combinations))%dopar%{
#  params <- parameter_combinations[j,]
#  rhouv <- as.numeric(Rhouv_Rhozw[[params$Rhouv_Rhozw]][1])
#  rhozw <- as.numeric(Rhouv_Rhozw[[params$Rhouv_Rhozw]][2])
#  case <- as.numeric(params$Case)
#  n_values <- as.numeric(params$N)
#  data_param = c(n_values, rhouv, rhozw)
  
  #simulate data 
#  simul <- simulate_data_3(data_param, g_sim_3, case)
  
  #compute the J_opt and return all the results
#  res_opt <- MC_j_opt_parallelized(n_MC, data_param, case, J_bs, J_ns, p_train, x_evaluation, degree, n_boot)
#  filename = paste("opt_2000", "_degree", degree, "_ptrain", p_train, "_nboot", n_boot, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
#  save(res_opt,file=filename)
  
#}

#stopCluster(cl)
