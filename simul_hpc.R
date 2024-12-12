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
case <- 3
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

indices_list <- seq(1, n_MC, by = 1)
already_done <- integer(0)

#check which things are already computed
for (n in 1:n_MC){
  filen <- paste("opt_", n_MC, "_degree", degree, "_ptrain", p_train, "_nboot", n_boot, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, "_simu", n, ".R" ,sep = "")
  if (file.exists(filen)) {
    # If the file exists, it means the simulation is already done
    already_done <- c(already_done, n)
  }
  }

new_indices_list <-  setdiff(indices_list, already_done)
n_MC_new <- length(new_indices_list)

# Parallelize the MC simulations
results <-foreach(i = 1:n_MC_new, .packages = c("splines", "MASS", "caret", "expm")) %dopar% {
  n = new_indices_list[i]
  
  # Simulate data for each iteration
  #simul <- simulate_data_3(data_param, g_sim_3, case)
  simul <- simul_all[[n]]
  tryCatch({
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
    filen <- paste("opt_", n_MC, "_degree", degree, "_ptrain", p_train, "_nboot", n_boot, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, "_simu", n, ".R" ,sep = "")
    res <- NULL
    save(res, file = filen)
    return(NULL)
  })
}
filename = paste("opt_", n_MC , "_degree", degree, "_ptrain", p_train, "_nboot", n_boot, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, ".R" ,sep = "")
save(results, file = filename)


stopCluster(cl)




