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

source("source_file_common_all.R")

#### Plot with the output of the algo ####
plot_mean_true <- function(for_curves, x_evaluation, case, method){
  if (method == "CVM_bs"){
    matrix_est <- for_curves$res_CVM_bs$list_est_values_CV_M_bs
  }
  if (method == "CVMSE_bs"){
    matrix_est <- for_curves$res_CVMSE_bs$list_est_values_CVMSE_bs
  }
  if (method == "lepski_bs"){
    matrix_est <- for_curves$res_lepski_bs$list_est_values_lepski_bs
  }
  if (method == "lepskiboot_bs"){
    matrix_est <- for_curves$res_lepskiboot_bs$list_est_values_lepskiboot_bs
  }
  
  if (method == "CVM_ns"){
    matrix_est <- for_curves$res_CVM_ns$list_est_values_CV_M_ns
  }
  if (method == "CVMSE_ns"){
    matrix_est <- for_curves$res_CVMSE_ns$list_est_values_CVMSE_ns
  }
  if (method == "lepski_ns"){
    matrix_est <- for_curves$res_lepski_ns$list_est_values_lepski_ns
  }
  if (method == "lepskiboot_ns"){
    matrix_est <- for_curves$res_lepskiboot_ns$list_est_values_lepskiboot_ns
  }
  
  avg <- colMeans(matrix_est)
  g_0_on_x = g_sim_3(x_evaluation, case)
  
  #n_MC = length(res_MC$list_gamma)
  #n_val = length(res_MC$list_g_hat_on_x[[1]])
  #n_data = length(res_MC$list_Z[[1]])
  
  data <- data.frame(x_evaluation, avg, g_0_on_x)
  graph <- ggplot(data, aes(x = x_evaluation)) +
    geom_line(aes(y = g_0_on_x, colour = 'True Function')) + # Changed the label
    geom_line(aes(y = avg, colour = 'Estimate by MC')) + # Changed the label
    ylab("Function value") +
    xlab("Evaluation points") + 
    scale_colour_manual(
      name = "Functions", # Change this to your desired legend title
      values = c("True Function" = "blue", "Estimate by MC" = "red") # Customize colors
    ) +
    ggtitle(paste("True vs avg"))
  return(graph)
}

plot_mean_true(for_curves, seq(-2, 2, length.out = 100), 2, "lepskiboot_bs")


plot_J_fixed <- function(j_index, J, bs_bool, res, x_evaluation, case){
  if (bs_bool == 1){
    matrix_to_plot <- res$values_bs_allJ[[j_index]]
  }
  else{
    matrix_to_plot <- res$values_ns_allJ[[j_index]]
  }
  
  if (nrow(matrix_to_plot)>0){
    avg <- colMeans(matrix_to_plot)
  }
  else{
    avg <- rep(0, 100)
  }
  g_0_on_x = g_sim_3(x_evaluation, case)
  
  data <- data.frame(x_evaluation, avg, g_0_on_x)
  graph <- ggplot(data, aes(x = x_evaluation)) +
    geom_line(aes(y = g_0_on_x, colour = 'True Function')) + # Changed the label
    geom_line(aes(y = avg, colour = 'Estimate by MC')) + # Changed the label
    ylab("Function value") +
    xlab("Evaluation points") + 
    scale_colour_manual(
      name = "Functions", # Change this to your desired legend title
      values = c("True Function" = "blue", "Estimate by MC" = "red") # Customize colors
    ) +
    ggtitle(paste("True vs avg, (J =", J, ", b =", bs_bool, ")"))
  return(graph)
}

plot_J_fixed(1, 5, 0, res, seq(-2, 2, length.out = 100), 2)

