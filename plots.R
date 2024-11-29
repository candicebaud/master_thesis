#### Source code ####
source("code_b_splines_monte_carlo.R")
library(splines)
library(MASS)
library(caret)
library(expm)
library(foreach)
library(doParallel)
library(ggplot2)

#### Import the data and plot, for fixed J ####
load("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_431821_incomplete/MC2000_fixedJ_J10_degree3_rhozw0.9_rhouv0.8_case3_n200.R")

# filter the things that we can not reuse 
zero_indices <- which(sapply(MC$list_gamma, function(x) all(x == 0)))
filtered_gamma <- MC$list_gamma[-zero_indices]
filtered_g_hat_on_x <- MC$list_g_hat_on_x[-zero_indices]
filtered_W <- MC$list_W[-zero_indices]
filtered_Y <- MC$list_Y[-zero_indices]
filtered_Z <- MC$list_Z[-zero_indices]

new_MC <- list(list_gamma = filtered_gamma, list_g_hat_on_x = filtered_g_hat_on_x, 
               list_W = filtered_W, list_Y = filtered_Y, list_Z = filtered_Z, g_0_on_x = MC$g_0_on_x)


x_evaluation = seq(-2, 2, length.out = 100)
J = 10
degree = 3
case = 3
rhozw = 0.9
rhouv = 0.8
plot_mean_true(new_MC, x_evaluation,J, degree, rhozw, rhouv, case) #: plots the average estimations and the true function
plot_allcurves_true(new_MC, x_evaluation,J, degree, rhozw, rhouv, case) #: plots all the estimated curves, the true one and the average one

load("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_431821_incomplete/perf2000_fixedJ_J10_degree3_rhozw0.9_rhouv0.8_case3_n2500.R")




#### Plot the distribution of the selected values of J_opt ####