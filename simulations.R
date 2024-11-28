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

#### Simulations avec le cluster (J fixé) ####
N_values <- c(200, 400, 1000, 2500)
Cases <- c(2, 3)
Rhouv_Rhozw <- list(c(0.5, 0.9), c(0.8, 0.9), c(0.8, 0.7))

parameter_combinations <- expand.grid(
  N = N_values,
  Case = Cases,
  Rhouv_Rhozw = seq_along(Rhouv_Rhozw) # Use indices for combinations
)

J_val <- c(4, 6, 10, 18, 34) 

source("/softs/R/createCluster .R")
cl <- createCluster()
registerDoParallel(cl)
degree = 3
x_evaluation = seq(-2, 2, length.out = 100)
n_MC = 2000


foreach (j=1:nrow(parameter_combinations))%dopar%{
  params <- parameter_combinations[j,]
  rhouv <- as.numeric(Rhouv_Rhozw[[params$Rhouv_Rhozw]][1])
  rhozw <- as.numeric(Rhouv_Rhozw[[params$Rhouv_Rhozw]][2])
  case <- as.numeric(params$Case)
  n_val <- as.numeric(params$N)
  data_param = c(n_val, rhouv, rhozw)
  for (i in 1:length(J_val)){
    J <- J_val[i]
    filename_MC <- paste ("MC2000_fixedJ", "_J", J, "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
    filename_perf <- paste ("perf2000_fixedJ", "_J", J, "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
    MC <- MC_fixed_J(J, n_MC, degree, x_evaluation, g_sim_3, case, data_param)
    perf_MC <- rep(0, 5)
    perf_MC[1] = compute_perf(MC, 'M')
    perf_MC[2] = compute_perf(MC, 'supnorm')
    perf_MC[3] = compute_perf(MC, 'Var')
    perf_MC[4] = compute_perf(MC, 'MSE')
    perf_MC[5] = compute_perf(MC, 'bias')
    save(MC,file=filename_MC)
    save(perf_MC,file=filename_perf)
  }
  
}

stopCluster(cl)








#### Simulation : degree = 3, J = 4, n = 200,rho_1 = 0.5, rho_2 = 0.9, case2 ####
# First, evaluate performances for different values of J 
#data_param = c(200, 0.5, 0.9)
#n_MC = 2000
#J = 4
#degree = 3
#x_evaluation = seq(-2, 2, length.out = 100)
#MC_J4_n200_deg3_rhocase1_gcase2 <- MC_fixed_J(J, n_MC, degree, x_evaluation, g_sim_3, 2, data_param)

#perf_MC_J4_n200_deg3_rhocase1_gcase2 <- rep(0, 5)
#perf_MC_J4_n200_deg3_rhocase1_gcase2[1] = compute_perf(MC_J4_n200_deg3_rhocase1_gcase2, 'M')
#perf_MC_J4_n200_deg3_rhocase1_gcase2[2] = compute_perf(MC_J4_n200_deg3_rhocase1_gcase2, 'supnorm')
#perf_MC_J4_n200_deg3_rhocase1_gcase2[3] = compute_perf(MC_J4_n200_deg3_rhocase1_gcase2, 'Var')
#perf_MC_J4_n200_deg3_rhocase1_gcase2[4] = compute_perf(MC_J4_n200_deg3_rhocase1_gcase2, 'MSE')
#perf_MC_J4_n200_deg3_rhocase1_gcase2[5] = compute_perf(MC_J4_n200_deg3_rhocase1_gcase2, 'bias')

#save.image(file = "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/perf_MC_J4_n200_deg3_rhocase1_gcase2.RData")

#load("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/perf_MC_J4_n200_deg3_rhocase1_gcase2.RData")

#### Simulation : degree = 3, J = 4, n = 400, rho_1 = 0.5, rho_2 = 0.9, case2 ####
#data_param = c(400, 0.5, 0.9)
#n_MC = 2000
#J = 4
#degree = 3
#x_evaluation = seq(-2, 2, length.out = 100)
#MC_J4_n400_deg3_rhocase1_gcase2 <- MC_fixed_J(J, n_MC, degree, x_evaluation, g_sim_3, 2, data_param)

#perf_MC_J4_n400_deg3_rhocase1_gcase2 <- rep(0, 5)
#perf_MC_J4_n400_deg3_rhocase1_gcase2[1] = compute_perf(MC_J4_n400_deg3_rhocase1_gcase2, 'M')
#perf_MC_J4_n400_deg3_rhocase1_gcase2[2] = compute_perf(MC_J4_n400_deg3_rhocase1_gcase2, 'supnorm')
#perf_MC_J4_n400_deg3_rhocase1_gcase2[3] = compute_perf(MC_J4_n400_deg3_rhocase1_gcase2, 'Var')
#perf_MC_J4_n400_deg3_rhocase1_gcase2[4] = compute_perf(MC_J4_n400_deg3_rhocase1_gcase2, 'MSE')
#perf_MC_J4_n400_deg3_rhocase1_gcase2[5] = compute_perf(MC_J4_n400_deg3_rhocase1_gcase2, 'bias')

#save.image(file = "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/perf_MC_J4_n400_deg3_rhocase1_gcase2.RData")



#### Simulation : degree = 3, J = 4, n = 1000, rho_1 = 0.5, rho_2 = 0.9, case2 ####
#data_param = c(1000, 0.5, 0.9)
#n_MC = 2000
#J = 4
#degree = 3
#x_evaluation = seq(-2, 2, length.out = 100)
#MC_J4_n1000_deg3_rhocase1_gcase2 <- MC_fixed_J(J, n_MC, degree, x_evaluation, g_sim_3, 2, data_param)

#perf_MC_J4_n1000_deg3_rhocase1_gcase2 <- rep(0, 5)
#perf_MC_J4_n1000_deg3_rhocase1_gcase2[1] = compute_perf(MC_J4_n1000_deg3_rhocase1_gcase2, 'M')
#perf_MC_J4_n1000_deg3_rhocase1_gcase2[2] = compute_perf(MC_J4_n1000_deg3_rhocase1_gcase2, 'supnorm')
#perf_MC_J4_n1000_deg3_rhocase1_gcase2[3] = compute_perf(MC_J4_n1000_deg3_rhocase1_gcase2, 'Var')
#perf_MC_J4_n1000_deg3_rhocase1_gcase2[4] = compute_perf(MC_J4_n1000_deg3_rhocase1_gcase2, 'MSE')
#perf_MC_J4_n1000_deg3_rhocase1_gcase2[5] = compute_perf(MC_J4_n1000_deg3_rhocase1_gcase2, 'bias')

#save.image(file = "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/perf_MC_J4_n1000_deg3_rhocase1_gcase2.RData")

#### Simulation : degree = 3, J = 4, n = 2500, rho_1 = 0.5, rho_2 = 0.9, case2 ####
#data_param = c(2500, 0.5, 0.9)
#n_MC = 2000
#J = 4
#degree = 3
#x_evaluation = seq(-2, 2, length.out = 100)
#MC_J4_n2500_deg3_rhocase1_gcase2 <- MC_fixed_J(J, n_MC, degree, x_evaluation, g_sim_3, 2, data_param)

#perf_MC_J4_n2500_deg3_rhocase1_gcase2 <- rep(0, 5)
#perf_MC_J4_n2500_deg3_rhocase1_gcase2[1] = compute_perf(MC_J4_n2500_deg3_rhocase1_gcase2, 'M')
#perf_MC_J4_n2500_deg3_rhocase1_gcase2[2] = compute_perf(MC_J4_n2500_deg3_rhocase1_gcase2, 'supnorm')
#perf_MC_J4_n2500_deg3_rhocase1_gcase2[3] = compute_perf(MC_J4_n2500_deg3_rhocase1_gcase2, 'Var')
#perf_MC_J4_n2500_deg3_rhocase1_gcase2[4] = compute_perf(MC_J4_n2500_deg3_rhocase1_gcase2, 'MSE')
#perf_MC_J4_n2500_deg3_rhocase1_gcase2[5] = compute_perf(MC_J4_n2500_deg3_rhocase1_gcase2, 'bias')

#save.image(file = "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/perf_MC_J4_n2500_deg3_rhocase1_gcase2.RData")

#### Simulation : degree = 3, J = 6, n = 200, rho_1 = 0.5, rho_2 = 0.9, case2 ####
#data_param = c(200, 0.5, 0.9)
#n_MC = 2000
#J = 6
#degree = 3
#x_evaluation = seq(-2, 2, length.out = 100)
#MC_J6_n200_deg3_rhocase1_gcase2 <- MC_fixed_J(J, n_MC, degree, x_evaluation, g_sim_3, 2, data_param)

#perf_MC_J6_n200_deg3_rhocase1_gcase2 <- rep(0, 5)
#perf_MC_J6_n200_deg3_rhocase1_gcase2[1] = compute_perf(MC_J6_n200_deg3_rhocase1_gcase2, 'M')
#perf_MC_J6_n200_deg3_rhocase1_gcase2[2] = compute_perf(MC_J6_n200_deg3_rhocase1_gcase2, 'supnorm')
#perf_MC_J6_n200_deg3_rhocase1_gcase2[3] = compute_perf(MC_J6_n200_deg3_rhocase1_gcase2, 'Var')
#perf_MC_J6_n200_deg3_rhocase1_gcase2[4] = compute_perf(MC_J6_n200_deg3_rhocase1_gcase2, 'MSE')
#perf_MC_J6_n200_deg3_rhocase1_gcase2[5] = compute_perf(MC_J6_n200_deg3_rhocase1_gcase2, 'bias')

#save.image(file = "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/perf_MC_J6_n200_deg3_rhocase1_gcase2.RData")


#### Simulation : degree = 3, J = 6, n = 400, rho_1 = 0.5, rho_2 = 0.9, case2 ####
#data_param = c(400, 0.5, 0.9)
#n_MC = 2000
#J = 6
#degree = 3
#x_evaluation = seq(-2, 2, length.out = 100)
#MC_J6_n400_deg3_rhocase1_gcase2 <- MC_fixed_J(J, n_MC, degree, x_evaluation, g_sim_3, 2, data_param)

#perf_MC_J6_n400_deg3_rhocase1_gcase2 <- rep(0, 5)
#perf_MC_J6_n400_deg3_rhocase1_gcase2[1] = compute_perf(MC_J6_n400_deg3_rhocase1_gcase2, 'M')
#perf_MC_J6_n400_deg3_rhocase1_gcase2[2] = compute_perf(MC_J6_n400_deg3_rhocase1_gcase2, 'supnorm')
#perf_MC_J6_n400_deg3_rhocase1_gcase2[3] = compute_perf(MC_J6_n400_deg3_rhocase1_gcase2, 'Var')
#perf_MC_J6_n400_deg3_rhocase1_gcase2[4] = compute_perf(MC_J6_n400_deg3_rhocase1_gcase2, 'MSE')
#perf_MC_J6_n400_deg3_rhocase1_gcase2[5] = compute_perf(MC_J6_n400_deg3_rhocase1_gcase2, 'bias')

#save.image(file = "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/perf_MC_J6_n400_deg3_rhocase1_gcase2.RData")


#### Simulation : degree = 3, J = 6, n = 1000, rho_1 = 0.5, rho_2 = 0.9, case2 ####
#data_param = c(1000, 0.5, 0.9)
#n_MC = 2000
#J = 6
#degree = 3
#x_evaluation = seq(-2, 2, length.out = 100)
#MC_J6_n1000_deg3_rhocase1_gcase2 <- MC_fixed_J(J, n_MC, degree, x_evaluation, g_sim_3, 2, data_param)

#perf_MC_J6_n1000_deg3_rhocase1_gcase2 <- rep(0, 5)
#perf_MC_J6_n1000_deg3_rhocase1_gcase2[1] = compute_perf(MC_J6_n1000_deg3_rhocase1_gcase2, 'M')
#perf_MC_J6_n1000_deg3_rhocase1_gcase2[2] = compute_perf(MC_J6_n1000_deg3_rhocase1_gcase2, 'supnorm')
#perf_MC_J6_n1000_deg3_rhocase1_gcase2[3] = compute_perf(MC_J6_n1000_deg3_rhocase1_gcase2, 'Var')
#perf_MC_J6_n1000_deg3_rhocase1_gcase2[4] = compute_perf(MC_J6_n1000_deg3_rhocase1_gcase2, 'MSE')
#perf_MC_J6_n1000_deg3_rhocase1_gcase2[5] = compute_perf(MC_J6_n1000_deg3_rhocase1_gcase2, 'bias')

#save.image(file = "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/perf_MC_J6_n1000_deg3_rhocase1_gcase2.RData")


#### Simulation : degree = 3, J = 6, n = 2500, rho_1 = 0.5, rho_2 = 0.9, case2 ####
#data_param = c(2500, 0.5, 0.9)
#n_MC = 2000
#J = 6
#degree = 3
#x_evaluation = seq(-2, 2, length.out = 100)
#MC_J6_n2500_deg3_rhocase1_gcase2 <- MC_fixed_J(J, n_MC, degree, x_evaluation, g_sim_3, 2, data_param)

#perf_MC_J6_n2500_deg3_rhocase1_gcase2 <- rep(0, 5)
#perf_MC_J6_n2500_deg3_rhocase1_gcase2[1] = compute_perf(MC_J6_n2500_deg3_rhocase1_gcase2, 'M')
#perf_MC_J6_n2500_deg3_rhocase1_gcase2[2] = compute_perf(MC_J6_n2500_deg3_rhocase1_gcase2, 'supnorm')
#perf_MC_J6_n2500_deg3_rhocase1_gcase2[3] = compute_perf(MC_J6_n2500_deg3_rhocase1_gcase2, 'Var')
#perf_MC_J6_n2500_deg3_rhocase1_gcase2[4] = compute_perf(MC_J6_n2500_deg3_rhocase1_gcase2, 'MSE')
#perf_MC_J6_n2500_deg3_rhocase1_gcase2[5] = compute_perf(MC_J6_n2500_deg3_rhocase1_gcase2, 'bias')

#save.image(file = "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/perf_MC_J6_n2500_deg3_rhocase1_gcase2.RData")






