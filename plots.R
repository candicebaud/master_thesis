#### Source code ####
#setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/master_thesis")
#source("code_b_splines_monte_carlo.R")
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

#### Functions to plot ####
library(ggplot2)
plot_mean_true <- function(res_MC, x_evaluation,J, degree, rhozw, rhouv, case){
  n_MC = length(res_MC$list_gamma)
  n_val = length(res_MC$list_g_hat_on_x[[1]])
  n_data = length(res_MC$list_Z[[1]])
  
  #compute the avg
  avg <- rep(0, n_val)
  for (x in 1:n_val){
    for (n in 1:n_MC){
      avg[x] <- avg[x] + res_MC$list_g_hat_on_x[[n]][x]/n_MC
    }}
  
  data <- data.frame(x_evaluation, avg, g_0_on_x = res_MC$g_0_on_x)
  graph <- ggplot(data, aes(x = x_evaluation)) +
    geom_line(aes(y = g_0_on_x, colour = 'True Function')) + # Changed the label
    geom_line(aes(y = avg, colour = 'Estimate by MC')) + # Changed the label
    ylab("Function value") +
    xlab("Evaluation points") + 
    scale_colour_manual(
      name = "Functions", # Change this to your desired legend title
      values = c("True Function" = "blue", "Estimate by MC" = "red") # Customize colors
    ) +
    ggtitle(paste("True vs avg, (n_MC =", n_MC, ", n_data =", n_data, ", rhowz =", rhozw, ", rhouv =", rhouv, ", degree =", degree, ", case = ", case, ", J =", J, ")"))
  return(graph)
}

#plot_mean_true(test, seq(-2, 2, by = 0.1), 10)

plot_allcurves_true <- function(res_MC, x_evaluation,J, degree, rhozw, rhouv, case){
  n_MC = length(res_MC$list_gamma)
  data_allcurves <- do.call(rbind, lapply(1:n_MC, function(i) data.frame(
    x = x_evaluation,
    y = res_MC$list_g_hat_on_x[[i]],
    id = i
  )))
  
  true_values <- data.frame(x = x_evaluation, y = res_MC$g_0_on_x)
  
  n_val = length(res_MC$list_g_hat_on_x[[1]])
  #compute the avg
  avg <- rep(0, n_val)
  for (x in 1:n_val){
    for (n in 1:n_MC){
      avg[x] <- avg[x] + res_MC$list_g_hat_on_x[[n]][x]/n_MC
    }}
  
  n_data = length(res_MC$list_Z[[1]])
  
  data_avg <- data.frame(x = x_evaluation, y = avg)
  
  ggplot() +
    geom_line(data = data_allcurves, aes(x = x, y = y, group = id,  color = "Estimated functions"), alpha = 0.5) +  # Estimations
    geom_line(data = true_values, aes(x = x, y = y, color = "True Function"), size = 1.2) +    # True function
    geom_line(data = data_avg, aes(x = x, y = y, color = "Average Estimate"), size = 1) +
    scale_color_manual(name = "Legend",
                       values = c("Estimated functions" = "grey", "True Function" = "green", "Average Estimate" = "black")) +
    theme_minimal() +
    labs(x = "Grid for estimation",
         y = "Estimated values of g")+
    ggtitle(paste("MC results, (n_MC =", n_MC, ", n_data =", n_data, ", rhowz =", rhozw, ", rhouv =", rhouv, ", degree =", degree, ", case = ", case, ", J =", J, ")"))
  
}


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








#### Plot the distribution of the selected values of J_opt ####
load("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451416_selectionJ_CVM_complete/opt_CV_M_2000_degree3_rhozw0.7_rhouv0.8_case2_n2500.R")
plot_mean_true(new_opt_CV_M, x_evaluation, mean(new_opt_CV_M$list_J_opt), 3, 0.7, 0.8, 2)
plot_allcurves_true(new_opt_CV_M, x_evaluation,mean(new_opt_CV_M$list_J_opt), 3, 0.7, 0.8, 2)

hist(new_opt_CV_M$list_J_opt)

length(new_opt_CV_M$list_gamma)

load("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_450196_running_selectionJ/selection_timeout/perf2000_lepski_boot_degree3_rhozw0.7_rhouv0.8_case3_n200.R")



#### Lepski boot ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/2")

load("simu_414886_lepskiboot/opt_lepski_boot_2000_grid50_degree3_rhozw0.9_rhouv0.8_case3_n2500.R")
plot_mean_true(new_opt_lepski_boot, seq(-2, 2, length.out = 100), mean(new_opt_lepski_boot$list_J_opt), 3, 0.7, 0.8, 2)
plot_allcurves_true(new_opt_lepski_boot, seq(-2, 2, length.out = 100),mean(new_opt_lepski_boot$list_J_opt), 3, 0.7, 0.8, 2)

#### Lepski ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/2")

load("simu_414887_lepski/opt_lepski_2000_degree3_rhozw0.9_rhouv0.8_case3_n200.R")
plot_mean_true(new_opt_lepski, seq(-2, 2, length.out = 100), mean(new_opt_lepski$list_J_opt), 3, 0.7, 0.8, 2)
plot_allcurves_true(new_opt_lepski, seq(-2, 2, length.out = 100),mean(new_opt_lepski), 3, 0.7, 0.8, 2)



