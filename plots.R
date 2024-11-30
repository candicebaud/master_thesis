#### Source code ####
source("code_b_splines_monte_carlo.R")
library(splines)
library(MASS)
library(caret)
library(expm)
library(foreach)
library(ggplot2)

#### Import the data and plot, for fixed J ####
load("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_441463_complete/MC2000_fixedJ_J4_degree3_rhozw0.9_rhouv0.8_case3_n2500.R")

# filter the things that we can not reuse 

x_evaluation = seq(-2, 2, length.out = 100)
J = 4
degree = 3
case = 3
rhozw = 0.9
rhouv = 0.8
plot_mean_true(new_MC, x_evaluation,J, degree, rhozw, rhouv, case) #: plots the average estimations and the true function
plot_allcurves_true(new_MC, x_evaluation,J, degree, rhozw, rhouv, case) #: plots all the estimated curves, the true one and the average one

load("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_431821_incomplete/perf2000_fixedJ_J10_degree3_rhozw0.9_rhouv0.8_case3_n2500.R")



#### Plot the distribution of the selected values of J_opt ####