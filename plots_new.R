setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/master_thesis")
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
library(patchwork)

source("source_file_common_all.R")

setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/5_final")
load("res_1")
load("res_2")
load("res_3")
load("res_4")
load("res_5")

#### Plot the true function and the different optimization methods ####
plot_methods_true <- function(for_curves, x_evaluation, case){
  matrix_CVM_ns <- for_curves$res_CVM_ns$list_est_values_CV_M_ns
  matrix_CVMSE_ns <- for_curves$res_CVMSE_ns$list_est_values_CVMSE_ns
  matrix_lepski_ns <- for_curves$res_lepski_ns$list_est_values_lepski_ns
  matrix_lepskiboot_ns <- for_curves$res_lepskiboot_ns$list_est_values_lepskiboot_ns
  
  matrix_CVM_bs <- for_curves$res_CVM_bs$list_est_values_CV_M_bs
  matrix_CVMSE_bs <- for_curves$res_CVMSE_bs$list_est_values_CVMSE_bs
  matrix_lepski_bs <- for_curves$res_lepski_bs$list_est_values_lepski_bs
  matrix_lepskiboot_bs <- for_curves$res_lepskiboot_bs$list_est_values_lepskiboot_bs
  
  avg_CVM_ns <- colMeans(matrix_CVM_ns)
  avg_CVMSE_ns <- colMeans(matrix_CVMSE_ns)
  avg_lepski_ns <- colMeans(matrix_lepski_ns)
  avg_lepskiboot_ns <- colMeans(matrix_lepskiboot_ns)
  
  avg_CVM_bs <- colMeans(matrix_CVM_bs)
  avg_CVMSE_bs <- colMeans(matrix_CVMSE_bs)
  avg_lepski_bs <- colMeans(matrix_lepski_bs)
  avg_lepskiboot_bs <- colMeans(matrix_lepskiboot_bs)

  g_0_on_x <- g_sim_3(x_evaluation, case)
  
  data_1 <- data.frame(x_evaluation, avg_CVM_bs, avg_CVMSE_bs, avg_lepski_bs, avg_lepskiboot_bs, g_0_on_x)
  graph_1 <- ggplot(data_1, aes(x = x_evaluation)) +
    geom_line(aes(y = g_0_on_x, colour = 'True Function'), size=1.1, show.legend = FALSE) + # Changed the label
    geom_line(aes(y = avg_CVM_bs, colour = 'Estimate by CVM'), show.legend = FALSE) + # Changed the label
    geom_line(aes(y = avg_CVMSE_bs, colour = 'Estimate by CVMSE'),show.legend = FALSE) + 
    geom_line(aes(y = avg_lepski_bs, colour = 'Estimate by Lepski'), show.legend = FALSE) + 
    geom_line(aes(y = avg_lepskiboot_bs, colour = 'Estimate by Lepski Boot'), show.legend = FALSE) + 
    ylab("Function value") +
    xlab("Evaluation points") + 
    scale_colour_manual(
      name = "Functions", # Change this to your desired legend title
      values = c("True Function" = "#009E73", "Estimate by CVM" = "#E69F00",
                 "Estimate by CVMSE" = "#56B4E9", "Estimate by Lepski" = "#661100",
                 "Estimate by Lepski Boot" = "#CC79A7") # Customize colors
    ) +
    ggtitle("B splines")
  
  data_2 <- data.frame(x_evaluation, avg_CVM_ns, avg_CVMSE_ns, avg_lepski_ns, avg_lepskiboot_ns, g_0_on_x)
  graph_2 <- ggplot(data_2, aes(x = x_evaluation)) +
    geom_line(aes(y = g_0_on_x, colour = 'True Function'), size=1.1) + # Changed the label
    geom_line(aes(y = avg_CVM_ns, colour = 'Estimate by CVM')) + # Changed the label
    geom_line(aes(y = avg_CVMSE_ns, colour = 'Estimate by CVMSE')) + 
    geom_line(aes(y = avg_lepski_ns, colour = 'Estimate by Lepski')) + 
    geom_line(aes(y = avg_lepskiboot_ns, colour = 'Estimate by Lepski Boot')) + 
    ylab("Function value") +
    xlab("Evaluation points") + 
    scale_colour_manual(
      name = "Functions", # Change this to your desired legend title
      values = c("True Function" = "#009E73", "Estimate by CVM" = "#E69F00",
                 "Estimate by CVMSE" = "#56B4E9", "Estimate by Lepski" = "#661100",
                 "Estimate by Lepski Boot" = "#CC79A7") # Customize colors
    ) +
    ggtitle("Natural splines")
  
  combined_graph <- graph_1 + graph_2 + plot_layout(ncol = 2)

  return(combined_graph)
}

x_eval = seq(-2, 2, length.out = 100)
plot_methods_true(res_1$for_curves, x_eval, 2)

x_eval = seq(-2, 2, length.out = 100)
plot_methods_true(res_2$for_curves, x_eval, 2)

x_eval = seq(-2, 2, length.out = 100)
plot_methods_true(res_3$for_curves, x_eval, 3)

x_eval = seq(-2, 2, length.out = 100)
plot_methods_true(res_4$for_curves, x_eval, 3)

x_eval = seq(-2, 2, length.out = 100)
plot_methods_true(res_5$for_curves, x_eval, 3)

#res_1
par(mfrow=c(2,2))
breaks = seq(1, 25, by = 1)
hist(res_1$all_lists$list_J_opt_CV_M_bs, breaks = breaks, main = "CV M", xlab = "J values", xlim = c(0, 20))
hist(res_1$all_lists$list_J_opt_CVMSE_bs, breaks = breaks, main = "CV MSE", xlab = "J values",  xlim = c(0, 20))
hist(res_1$all_lists$list_J_opt_lepski_bs, breaks = breaks, main = "Lepski", xlab = "J values", xlim = c(0, 20))
hist(res_1$all_lists$list_J_opt_lepskiboot_bs, breaks = breaks, main = "Lepskiboot", xlab = "J values", xlim = c(0, 20))


par(mfrow=c(2,2))
hist(res_1$all_lists$list_J_opt_CV_M_ns, breaks = breaks, main = "CV M", xlab = "J values", xlim = c(0, 20))
hist(res_1$all_lists$list_J_opt_CVMSE_ns, breaks = breaks, main = "CV MSE", xlab = "J values", xlim = c(0, 20))
hist(res_1$all_lists$list_J_opt_lepski_ns,breaks = breaks, main = "Lepski", xlab = "J values", xlim = c(0, 20))
hist(res_1$all_lists$list_J_opt_lepskiboot_ns,breaks = breaks, main = "Lepskiboot", xlab = "J values", xlim = c(0, 20))

#res_2
par(mfrow=c(2,2))
breaks = seq(1, 25, by = 1)
hist(res_2$all_lists$list_J_opt_CV_M_bs, breaks = breaks, main = "CV M", xlab = "J values", xlim = c(0, 20))
hist(res_2$all_lists$list_J_opt_CVMSE_bs, breaks = breaks, main = "CV MSE", xlab = "J values",  xlim = c(0, 20))
hist(res_2$all_lists$list_J_opt_lepski_bs, breaks = breaks, main = "Lepski", xlab = "J values", xlim = c(0, 20))
hist(res_2$all_lists$list_J_opt_lepskiboot_bs, breaks = breaks, main = "Lepskiboot", xlab = "J values", xlim = c(0, 20))


par(mfrow=c(2,2))
hist(res_2$all_lists$list_J_opt_CV_M_ns, breaks = breaks, main = "CV M", xlab = "J values", xlim = c(0, 20))
hist(res_2$all_lists$list_J_opt_CVMSE_ns, breaks = breaks, main = "CV MSE", xlab = "J values", xlim = c(0, 20))
hist(res_2$all_lists$list_J_opt_lepski_ns,breaks = breaks, main = "Lepski", xlab = "J values", xlim = c(0, 20))
hist(res_2$all_lists$list_J_opt_lepskiboot_ns,breaks = breaks, main = "Lepskiboot", xlab = "J values", xlim = c(0, 20))


#res_3
par(mfrow=c(2,2))
breaks = seq(1, 25, by = 1)
hist(res_3$all_lists$list_J_opt_CV_M_bs, breaks = breaks, main = "CV M", xlab = "J values", xlim = c(0, 20))
hist(res_3$all_lists$list_J_opt_CVMSE_bs, breaks = breaks, main = "CV MSE", xlab = "J values",  xlim = c(0, 20))
hist(res_3$all_lists$list_J_opt_lepski_bs, breaks = breaks, main = "Lepski", xlab = "J values", xlim = c(0, 20))
hist(res_3$all_lists$list_J_opt_lepskiboot_bs, breaks = breaks, main = "Lepskiboot", xlab = "J values", xlim = c(0, 20))


par(mfrow=c(2,2))
hist(res_3$all_lists$list_J_opt_CV_M_ns, breaks = breaks, main = "CV M", xlab = "J values", xlim = c(0, 20))
hist(res_3$all_lists$list_J_opt_CVMSE_ns, breaks = breaks, main = "CV MSE", xlab = "J values", xlim = c(0, 20))
hist(res_3$all_lists$list_J_opt_lepski_ns,breaks = breaks, main = "Lepski", xlab = "J values", xlim = c(0, 20))
hist(res_3$all_lists$list_J_opt_lepskiboot_ns,breaks = breaks, main = "Lepskiboot", xlab = "J values", xlim = c(0, 20))


#res_4
par(mfrow=c(2,2))
breaks = seq(1, 25, by = 1)
hist(res_4$all_lists$list_J_opt_CV_M_bs, breaks = breaks, main = "CV M", xlab = "J values", xlim = c(0, 20))
hist(res_4$all_lists$list_J_opt_CVMSE_bs, breaks = breaks, main = "CV MSE", xlab = "J values",  xlim = c(0, 20))
hist(res_4$all_lists$list_J_opt_lepski_bs, breaks = breaks, main = "Lepski", xlab = "J values", xlim = c(0, 20))
hist(res_4$all_lists$list_J_opt_lepskiboot_bs, breaks = breaks, main = "Lepskiboot", xlab = "J values", xlim = c(0, 20))


par(mfrow=c(2,2))
hist(res_4$all_lists$list_J_opt_CV_M_ns, breaks = breaks, main = "CV M", xlab = "J values", xlim = c(0, 20))
hist(res_4$all_lists$list_J_opt_CVMSE_ns, breaks = breaks, main = "CV MSE", xlab = "J values", xlim = c(0, 20))
hist(res_4$all_lists$list_J_opt_lepski_ns,breaks = breaks, main = "Lepski", xlab = "J values", xlim = c(0, 20))
hist(res_4$all_lists$list_J_opt_lepskiboot_ns,breaks = breaks, main = "Lepskiboot", xlab = "J values", xlim = c(0, 20))


#res_5
par(mfrow=c(2,2))
breaks = seq(1, 25, by = 1)
hist(res_5$all_lists$list_J_opt_CV_M_bs, breaks = breaks, main = "CV M", xlab = "J values", xlim = c(0, 20))
hist(res_5$all_lists$list_J_opt_CVMSE_bs, breaks = breaks, main = "CV MSE", xlab = "J values",  xlim = c(0, 20))
hist(res_5$all_lists$list_J_opt_lepski_bs, breaks = breaks, main = "Lepski", xlab = "J values", xlim = c(0, 20))
hist(res_5$all_lists$list_J_opt_lepskiboot_bs, breaks = breaks, main = "Lepskiboot", xlab = "J values", xlim = c(0, 20))


par(mfrow=c(2,2))
hist(res_5$all_lists$list_J_opt_CV_M_ns, breaks = breaks, main = "CV M", xlab = "J values", xlim = c(0, 20))
hist(res_5$all_lists$list_J_opt_CVMSE_ns, breaks = breaks, main = "CV MSE", xlab = "J values", xlim = c(0, 20))
hist(res_5$all_lists$list_J_opt_lepski_ns,breaks = breaks, main = "Lepski", xlab = "J values", xlim = c(0, 20))
hist(res_5$all_lists$list_J_opt_lepskiboot_ns,breaks = breaks, main = "Lepskiboot", xlab = "J values", xlim = c(0, 20))




