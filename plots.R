#### Source code ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/master_thesis")
source("code_b_splines_monte_carlo.R")
library(splines)
library(MASS)
library(caret)
library(expm)
library(foreach)
library(ggplot2)


#### Import the data and plot, for fixed J ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_441463_complete_allJ_benchmark")
load("MC2000_fixedJ_J4_degree3_rhozw0.9_rhouv0.8_case3_n200.R")

# filter the things that we can not reuse 

x_evaluation = seq(-2, 2, length.out = 100)
J = 4
degree = 3
case = 3
rhozw = 0.9
rhouv = 0.8
plot_mean_true(new_MC, x_evaluation,J, degree, rhozw, rhouv, case) #: plots the average estimations and the true function
plot_allcurves_true(new_MC, x_evaluation,J, degree, rhozw, rhouv, case) #: plots all the estimated curves, the true one and the average one


#### Matrix of performances, fixed J ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_441463_complete_allJ_benchmark")

#load the data
degree = 3
x_evaluation = seq(-2, 2, length.out = 100)
n_MC = 2000

N_values <- c(200, 400, 1000, 2500)
Cases <- c(2, 3)
Rhouv_Rhozw <- list(c(0.5, 0.9), c(0.8, 0.9), c(0.8, 0.7))

parameter_combinations <- expand.grid(
  N = N_values,
  Case = Cases,
  Rhouv_Rhozw = seq_along(Rhouv_Rhozw) # Use indices for combinations
)

J_val <- c(4, 6, 10, 18, 34) 

for (j in 1:nrow(parameter_combinations)){
  params <- parameter_combinations[j,]
  rhouv <- as.numeric(Rhouv_Rhozw[[params$Rhouv_Rhozw]][1])
  rhozw <- as.numeric(Rhouv_Rhozw[[params$Rhouv_Rhozw]][2])
  case <- as.numeric(params$Case)
  n_val <- as.numeric(params$N)
  data_param = c(n_val, rhouv, rhozw)
  for (i in 1:length(J_val)){
    J <- J_val[i]
    filename_perf <- paste("perf2000_fixedJ", "_J", J, "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
    load(filename_perf)
    perf_MC <- c(filename_perf, perf_MC)
    assign(filename_perf, perf_MC)
    }
}

#create a df with aggregated performances
vector_names <- ls(pattern = "^perf2000")
vector_list <- mget(vector_names)
df <- as.data.frame(do.call(rbind, vector_list))


library(dplyr)
library(stringr)

df_ <- df %>%
  mutate(
    J = str_extract(V1, "(?<=_J)\\d+"),
    degree = str_extract(V1, "(?<=_degree)\\d+"),
    rhozw = str_extract(V1, "(?<=_rhozw)\\d*\\.?\\d+"),
    rhouv = str_extract(V1, "(?<=_rhouv)\\d*\\.?\\d+"),
    case = str_extract(V1, "(?<=_case)\\d+"),
    n_val = str_extract(V1, "(?<=_n)\\d+")
  )



  
#### Plot the distribution of the selected values of J_opt ####
load("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_450196_running_selectionJ/selection_timeout/opt_lepski_boot_2000_grid1000_degree3_rhozw0.7_rhouv0.8_case3_n200.R")
plot_mean_true(new_opt_lepski_boot, x_evaluation, mean(new_opt_lepski_boot$list_J_opt), 3, 0.7, 0.8, 2)
plot_allcurves_true(new_opt_lepski_boot, x_evaluation,mean(new_opt_lepski_boot$list_J_opt), 3, 0.7, 0.8, 2)

hist(new_opt_lepski_boot$list_J_opt)
load("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_450196_running_selectionJ/selection_timeout/perf2000_lepski_boot_degree3_rhozw0.7_rhouv0.8_case3_n200.R")
