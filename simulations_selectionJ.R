#### Import code ####
#setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/master_thesis")
source("code_b_splines_monte_carlo.R")
library(splines)
library(MASS)
library(caret)
library(expm)
library(foreach)
library(doParallel)
library(ggplot2)

#### Simulations ####
# on fait la sélection de paramètre, quand on plotera on regardera qui on enlève pcq le gamma = 0

source("/softs/R/createCluster.R")
cl <- createCluster()
registerDoParallel(cl)
N_values <- c(200, 400, 1000, 2500)
Cases <- c(2, 3)
Rhouv_Rhozw <- list(c(0.5, 0.9), c(0.8, 0.9), c(0.8, 0.7))

parameter_combinations <- expand.grid(
  N = N_values,
  Case = Cases,
  Rhouv_Rhozw = seq_along(Rhouv_Rhozw) # Use indices for combinations
)

degree = 3
x_evaluation = seq(-2, 2, length.out = 100)
n_MC = 2000

J_val_CV <- c(4, 6, 10, 18, 34) 

p_train = 0.8

foreach (j=1:nrow(parameter_combinations))%dopar%{
  params <- parameter_combinations[j,]
  rhouv <- as.numeric(Rhouv_Rhozw[[params$Rhouv_Rhozw]][1])
  rhozw <- as.numeric(Rhouv_Rhozw[[params$Rhouv_Rhozw]][2])
  case <- as.numeric(params$Case)
  n_val <- as.numeric(params$N)
  data_param = c(n_val, rhouv, rhozw)
  
  #selection par CV_M
  opt_CV_M <- MC_CV('CV_M', n_MC, J_val_CV, p_train, degree, x_evaluation, g_sim_3, case, data_param)
  new_opt_CV_M <- compute_new_MC_selection(opt_CV_M)
  filename = paste("opt_CV_M_2000", "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
  save(new_opt_CV_M,file=filename)
  perf_CV_M <- rep(0, 5)
  perf_CV_M[1] = compute_perf(new_opt_CV_M, 'M')
  perf_CV_M[2] = compute_perf(new_opt_CV_M, 'supnorm')
  perf_CV_M[3] = compute_perf(new_opt_CV_M, 'Var')
  perf_CV_M[4] = compute_perf(new_opt_CV_M, 'MSE')
  perf_CV_M[5] = compute_perf(new_opt_CV_M, 'bias')
  filename_perf <- paste ("perf2000_CV_M", "_J", J, "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
  save(perf_CV_M,file=filename_perf)
  
  #selection par CV_MSE
  opt_CV_MSE <- MC_CV('CV_MSE', n_MC, J_val_CV, p_train, degree, x_evaluation, g_sim_3, case, data_param)
  new_opt_CV_MSE <- compute_new_MC_selection(opt_CV_MSE)
  filename = paste("opt_CV_MSE_2000", "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
  save(new_opt_CV_MSE,file=filename)
  perf_CV_MSE <- rep(0, 5)
  perf_CV_MSE[1] = compute_perf(new_opt_CV_MSE, 'M')
  perf_CV_MSE[2] = compute_perf(new_opt_CV_MSE, 'supnorm')
  perf_CV_MSE[3] = compute_perf(new_opt_CV_MSE, 'Var')
  perf_CV_MSE[4] = compute_perf(new_opt_CV_MSE, 'MSE')
  perf_CV_MSE[5] = compute_perf(new_opt_CV_MSE, 'bias')
  filename_perf <- paste ("perf2000_CV_MSE", "_J", J, "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
  save(perf_CV_MSE,file=filename_perf)
  
  #selection by lepski boot
  opt_lepski_boot <- MC_lepski_boot(n_MC, 1000, x_evaluation, valid_dim, degree, g_sim_3, case, data_param)
  new_opt_lepski_boot <- compute_new_MC_selection(opt_lepski_boot)
  filename = paste ("opt_lepski_boot_2000", "_grid", 1000, "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
  save(new_opt_lepski_boot,file=filename)
  perf_lepski_boot <- rep(0, 5)
  perf_lepski_boot[1] = compute_perf(new_opt_lepski_boot, 'M')
  perf_lepski_boot[2] = compute_perf(new_opt_lepski_boot, 'supnorm')
  perf_lepski_boot[3] = compute_perf(new_opt_lepski_boot, 'Var')
  perf_lepski_boot[4] = compute_perf(new_opt_lepski_boot, 'MSE')
  perf_lepski_boot[5] = compute_perf(new_opt_lepski_boot, 'bias')
  filename_perf <- paste ("perf2000_lepski_boot", "_J", J, "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
  save(perf_lepski_boot,file=filename_perf)
  
  
  #selection by lepski
  opt_lepski <- MC_lepski(n_MC, x_evaluation, degree, valid_dim, g_sim_3, case, data_param)
  new_opt_lepski <- compute_new_MC_selection(opt_lepski)
  filename = paste ("opt_lepski_2000", "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
  save(new_opt_lepski,file=filename)
  perf_lepski <- rep(0, 5)
  perf_lepski[1] = compute_perf(new_opt_lepski, 'M')
  perf_lepski[2] = compute_perf(new_opt_lepski, 'supnorm')
  perf_lepski[3] = compute_perf(new_opt_lepski, 'Var')
  perf_lepski[4] = compute_perf(new_opt_lepski, 'MSE')
  perf_lepski[5] = compute_perf(new_opt_lepski, 'bias')
  filename_perf <- paste ("perf2000_lepski", "_J", J, "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
  save(perf_lepski,file=filename_perf)

}

stopCluster(cl)


#perf_MC <- rep(0, 5)
#perf_MC[1] = compute_perf(opt_CV_M, 'M')
#perf_MC[2] = compute_perf(opt_CV_M, 'supnorm')
#perf_MC[3] = compute_perf(opt_CV_M, 'Var')
#perf_MC[4] = compute_perf(opt_CV_M, 'MSE')
#perf_MC[5] = compute_perf(opt_CV_M, 'bias')
#filename = paste ("perf_CV_M_2000", "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
#save(perf_MC, file=filename)


