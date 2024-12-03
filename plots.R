#### Source code ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/master_thesis")
source("code_b_splines_monte_carlo.R")
library(splines)
library(MASS)
library(caret)
library(expm)
library(foreach)
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
library(xtable)
library(kableExtra)


#### Courbes CVM et CVM NS ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451438_selectionJ_CVM_NS_complete")
load("NS_opt_CV_M_2000_degree3_rhozw0.7_rhouv0.8_case3_n2500.R")
plot_mean_true(new_opt_CV_M, seq(-2, 2, length.out = 100), 'none', 3, 0.7, 0.8, 3)
plot_allcurves_true(new_opt_CV_M, seq(-2, 2, length.out = 100), 'none', 3, 0.7, 0.8, 3)

load("NS_opt_CV_M_2000_degree3_rhozw0.7_rhouv0.8_case2_n200.R")
plot_mean_true(new_opt_CV_M, seq(-2, 2, length.out = 100), 'none', 3, 0.7, 0.8, 2)
plot_allcurves_true(new_opt_CV_M, seq(-2, 2, length.out = 100), 'none', 3, 0.7, 0.8, 2)


setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451416_selectionJ_CVM_complete")
load("opt_CV_M_2000_degree3_rhozw0.7_rhouv0.8_case2_n200.R")
plot_mean_true(new_opt_CV_M, seq(-2, 2, length.out = 100), 'none', 3, 0.7, 0.8, 2)
plot_allcurves_true(new_opt_CV_M, seq(-2, 2, length.out = 100), 'none', 3, 0.7, 0.8, 2)

#p1 <- plot_mean_true(new_opt_CV_M, seq(-2, 2, length.out = 100), 'none', 3, 0.7, 0.8, 2)
#p2 <- plot_allcurves_true(new_opt_CV_M, seq(-2, 2, length.out = 100), 'none', 3, 0.7, 0.8, 2)
#combined_plot <- p1 + p2
#print(combined_plot)


#### Import the data and plot, for fixed J ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451423_allJ_NS")
load("NS_MC2000_fixedJ_J4_degree3_rhozw0.9_rhouv0.5_case2_n200.R")

x_evaluation = seq(-2, 2, length.out = 100)
J = 4
degree = 3
case = 2
rhozw = 0.9
rhouv = 0.5
plot_mean_true(new_MC, x_evaluation,J, degree, rhozw, rhouv, case) #: plots the average estimations and the true function
plot_allcurves_true(new_MC, x_evaluation,J, degree, rhozw, rhouv, case) #: plots all the estimated curves, the true one and the average one







#### Compare selection methods ####

#### Plot the distribution of the selected values of J_opt ####
load("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451416_selectionJ_CVM_complete/opt_CV_M_2000_degree3_rhozw0.7_rhouv0.8_case2_n2500.R")
plot_mean_true(new_opt_CV_M, x_evaluation, mean(new_opt_CV_M$list_J_opt), 3, 0.7, 0.8, 2)
plot_allcurves_true(new_opt_CV_M, x_evaluation,mean(new_opt_CV_M$list_J_opt), 3, 0.7, 0.8, 2)

hist(new_opt_CV_M$list_J_opt)

length(new_opt_CV_M$list_gamma)

load("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_450196_running_selectionJ/selection_timeout/perf2000_lepski_boot_degree3_rhozw0.7_rhouv0.8_case3_n200.R")




#### Performance Lepski ####
load("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451418_selectionJ_lepski/opt_lepski_2000_degree3_rhozw0.9_rhouv0.8_case3_n200.R")

unique(new_opt_lepski$list_J_opt)
length(new_opt_lepski$list_gamma)
plot_mean_true(new_opt_lepski, seq(-2, 2, length.out = 100), 'non', 3, 0.7, 0.8, 2)
plot_allcurves_true(new_opt_lepski, seq(-2, 2, length.out = 100), 'none', 3, 0.7, 0.8, 2)


load("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451418_selectionJ_lepski/opt_lepski_2000_degree3_rhozw0.7_rhouv0.8_case3_n200.R")
length(new_opt_lepski$list_gamma)
plot_mean_true(new_opt_lepski, seq(-2, 2, length.out = 100), 'non', 3, 0.7, 0.8, 2)
plot_allcurves_true(new_opt_lepski, seq(-2, 2, length.out = 100), 'none', 3, 0.7, 0.8, 2)

hist(new_opt_lepski$list_J_opt)

unique(new_opt_lepski$list_J_opt)