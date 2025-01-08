setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/master_thesis")
source("source_file_common_all.R")

##### Create subintervals #####
x_eval = seq(-2, 2, length.out = 100)

# on divise en 10 sous intervalles [-a, a]
x_eval_neg = x_eval[1:50]
x_eval_pos = x_eval[51:100]

a_values = x_eval_pos[10 * (1:5)]

subint <- matrix(0, nrow = 5, ncol = 2)
subint_indices <- matrix(0, nrow = 5, ncol = 2)
for (i in 1:5){
  subint[i,] <- c(-a_values[i], a_values[i])
  #minus_a = -a_values[i]
  plus_a = a_values[i]
  index_plus_a = which(x_eval == plus_a)
  index_minus_a = 100-index_plus_a+1
  subint_indices[i,]<- c(index_minus_a, index_plus_a)
}

#### Compute bias #### 
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/5_final")

compute_MSE <- function(vec_1, vec_2){
  return (mean((vec_1 - vec_2)**2))
}

compute_vec <- function(matrix, subint_index){
  avg_est_values = colMeans(matrix)[subint_index[1]:subint_index[2]]
  return(avg_est_values)
}

compute_all_bias_subint <- function(all_lists, subint_indices, g_0_on_x){
  n_subint = nrow(subint_indices)
  bias_CVM_bs <- c()
  bias_CVMSE_bs <- c()
  bias_lepski_bs <- c()
  bias_lepskiboot_bs <- c()
  
  bias_CVM_ns <- c()
  bias_CVMSE_ns <- c()
  bias_lepski_ns <- c()
  bias_lepskiboot_ns <- c()
  
  for (i in 1:n_subint){
    subint_index = subint_indices[i,]
    
    bias_CVM_bs[i] <- compute_MSE(compute_vec(all_lists$list_est_values_CV_M_bs, subint_index), g_0_on_x[subint_index[1]:subint_index[2]])
    bias_CVM_ns[i] <- compute_MSE(compute_vec(all_lists$list_est_values_CV_M_ns, subint_index), g_0_on_x[subint_index[1]:subint_index[2]])
    
    bias_CVMSE_bs[i] <- compute_MSE(compute_vec(all_lists$list_est_values_CVMSE_bs, subint_index), g_0_on_x[subint_index[1]:subint_index[2]])
    bias_CVMSE_ns[i] <- compute_MSE(compute_vec(all_lists$list_est_values_CVMSE_ns, subint_index), g_0_on_x[subint_index[1]:subint_index[2]])
    
    bias_lepski_bs[i] <- compute_MSE(compute_vec(all_lists$list_est_values_lepski_bs, subint_index), g_0_on_x[subint_index[1]:subint_index[2]])
    bias_lepski_ns[i] <- compute_MSE(compute_vec(all_lists$list_est_values_lepski_ns, subint_index), g_0_on_x[subint_index[1]:subint_index[2]])
    
    bias_lepskiboot_bs[i] <- compute_MSE(compute_vec(all_lists$list_est_values_lepskiboot_bs, subint_index), g_0_on_x[subint_index[1]:subint_index[2]])
    bias_lepskiboot_ns[i] <- compute_MSE(compute_vec(all_lists$list_est_values_lepskiboot_ns, subint_index), g_0_on_x[subint_index[1]:subint_index[2]])
  }
  
  return(list(bias_CVM_bs = bias_CVM_bs, bias_CVM_ns = bias_CVM_ns,
              bias_CVMSE_bs = bias_CVMSE_bs, bias_CVMSE_ns = bias_CVMSE_ns, 
              bias_lepski_bs = bias_lepski_bs, bias_lepski_ns = bias_lepski_ns,
              bias_lepskiboot_bs = bias_lepskiboot_bs, bias_lepskiboot_ns = bias_lepskiboot_ns))
}


load('res_1')
bias_res_1 = compute_all_bias_subint(res_1$all_lists, subint_indices, g_sim_3(x_eval, 2))
save(bias_res_1, file = "bias_res_1_intervals")

load('res_2')
bias_res_2 = compute_all_bias_subint(res_2$all_lists, subint_indices, g_sim_3(x_eval, 2))
save(bias_res_2, file = "bias_res_2_intervals")

load('res_3')
bias_res_3 = compute_all_bias_subint(res_3$all_lists, subint_indices, g_sim_3(x_eval, 3))
save(bias_res_3, file = "bias_res_3_intervals")

load('res_4')
bias_res_4 = compute_all_bias_subint(res_4$all_lists, subint_indices, g_sim_3(x_eval, 3))
save(bias_res_4, file = "bias_res_4_intervals")

load('res_5')
bias_res_5 = compute_all_bias_subint(res_5$all_lists, subint_indices, g_sim_3(x_eval, 3))
save(bias_res_5, file = "bias_res_5_intervals")


##### Compute MSE ####
compute_all_MSE_subint <- function(all_lists, subint_indices, g_0_on_x){
  n_subint = nrow(subint_indices)
  n_MC <- 2000
  MSE_CVM_bs <- matrix(0, nrow = n_MC, ncol = n_subint)
  MSE_CVMSE_bs <- matrix(0, nrow = n_MC, ncol = n_subint)
  MSE_lepski_bs <- matrix(0, nrow = n_MC, ncol = n_subint)
  MSE_lepskiboot_bs <- matrix(0, nrow = n_MC, ncol = n_subint)
  
  MSE_CVM_ns <- matrix(0, nrow = n_MC, ncol = n_subint)
  MSE_CVMSE_ns <- matrix(0, nrow = n_MC, ncol = n_subint)
  MSE_lepski_ns <- matrix(0, nrow = n_MC, ncol = n_subint)
  MSE_lepskiboot_ns <- matrix(0, nrow = n_MC, ncol = n_subint)
  
  for (i in 1:n_subint){
    subint_index = subint_indices[i,]
    for (n in 1:n_MC){
      MSE_CVM_bs[n,i] = compute_MSE(all_lists$list_est_values_CV_M_bs[n,][subint_index[1]:subint_index[2]], g_0_on_x[subint_index[1]:subint_index[2]])
      MSE_CVM_ns[n,i] = compute_MSE(all_lists$list_est_values_CV_M_ns[n,][subint_index[1]:subint_index[2]], g_0_on_x[subint_index[1]:subint_index[2]])
      
      MSE_CVMSE_bs[n,i] = compute_MSE(all_lists$list_est_values_CVMSE_bs[n,][subint_index[1]:subint_index[2]], g_0_on_x[subint_index[1]:subint_index[2]])
      MSE_CVMSE_ns[n,i] = compute_MSE(all_lists$list_est_values_CVMSE_ns[n,][subint_index[1]:subint_index[2]], g_0_on_x[subint_index[1]:subint_index[2]])
      
      MSE_lepski_bs[n,i] = compute_MSE(all_lists$list_est_values_lepski_bs[n,][subint_index[1]:subint_index[2]], g_0_on_x[subint_index[1]:subint_index[2]])
      MSE_lepski_ns[n,i] = compute_MSE(all_lists$list_est_values_lepski_ns[n,][subint_index[1]:subint_index[2]], g_0_on_x[subint_index[1]:subint_index[2]])
      
      MSE_lepskiboot_bs[n,i] = compute_MSE(all_lists$list_est_values_lepskiboot_bs[n,][subint_index[1]:subint_index[2]], g_0_on_x[subint_index[1]:subint_index[2]])
      MSE_lepskiboot_ns[n,i] = compute_MSE(all_lists$list_est_values_lepskiboot_ns[n,][subint_index[1]:subint_index[2]], g_0_on_x[subint_index[1]:subint_index[2]])
      
    }}
  
  
  return(list(MSE_CVM_bs = colMeans(MSE_CVM_bs), MSE_CVM_ns = colMeans(MSE_CVM_ns),
              MSE_CVMSE_bs = colMeans(MSE_CVMSE_bs), MSE_CVMSE_ns = colMeans(MSE_CVMSE_ns), 
              MSE_lepski_bs = colMeans(MSE_lepski_bs), MSE_lepski_ns = colMeans(MSE_lepski_ns),
              MSE_lepskiboot_bs = colMeans(MSE_lepskiboot_bs), MSE_lepskiboot_ns = colMeans(MSE_lepskiboot_ns)))
}


MSE_res_1 <- compute_all_MSE_subint(res_1$all_lists, subint_indices, g_sim_3(x_eval, 2))
save(MSE_res_1, file = "MSE_res_1_interval")
MSE_res_2 <- compute_all_MSE_subint(res_2$all_lists, subint_indices, g_sim_3(x_eval, 2))
save(MSE_res_2, file = "MSE_res_2_interval")
MSE_res_3 <- compute_all_MSE_subint(res_3$all_lists, subint_indices, g_sim_3(x_eval, 3))
save(MSE_res_3, file = "MSE_res_3_interval")
MSE_res_4 <- compute_all_MSE_subint(res_4$all_lists, subint_indices, g_sim_3(x_eval, 3))
save(MSE_res_4, file = "MSE_res_4_interval")
MSE_res_5 <- compute_all_MSE_subint(res_5$all_lists, subint_indices, g_sim_3(x_eval, 3))
save(MSE_res_5, file = "MSE_res_5_interval")


##### Variance ####
diff_function <- function(list_1, list_2){
  n_elem = length(list_1)
  diff_list <- list()
  for (i in 1:n_elem){
    diff_list[[i]] <- list_2[[i]] - list_1[[i]]
  }
  return(diff_list)
}

Var_res_1 <- diff_function(bias_res_1, MSE_res_1)
Var_res_2 <- diff_function(bias_res_2, MSE_res_2)
Var_res_3 <- diff_function(bias_res_3, MSE_res_3)
Var_res_4 <- diff_function(bias_res_4, MSE_res_4)
Var_res_5 <- diff_function(bias_res_5, MSE_res_5)

save(Var_res_1, file = "Var_res_1_intervals")
save(Var_res_2, file = "Var_res_2_intervals")
save(Var_res_3, file = "Var_res_3_intervals")
save(Var_res_4, file = "Var_res_4_intervals")
save(Var_res_5, file = "Var_res_5_intervals")

#### Plots ####
library(ggplot2)

df_res1_2 <- data.frame(MSE_res_1$MSE_CVM_bs, MSE_res_2$MSE_CVM_bs, 
           MSE_res_1$MSE_CVMSE_bs, MSE_res_2$MSE_CVMSE_bs,
           MSE_res_1$MSE_lepski_bs, MSE_res_2$MSE_lepski_bs,
           MSE_res_1$MSE_lepskiboot_bs, MSE_res_2$MSE_lepskiboot_bs,
           MSE_res_1$MSE_CVM_ns, MSE_res_2$MSE_CVM_ns, 
           MSE_res_1$MSE_CVMSE_ns, MSE_res_2$MSE_CVMSE_ns,
           MSE_res_1$MSE_lepski_ns, MSE_res_2$MSE_lepski_ns,
           MSE_res_1$MSE_lepskiboot_ns, MSE_res_2$MSE_lepskiboot_ns)

colnames(df_res1_2) <- c("MSE_CVM_bs_1", "MSE_CVM_bs_2",
                  "MSE_CVMSE_bs_1", "MSE_CVMSE_bs_2",
                  "MSE_lepski_bs_1", "MSE_lepski_bs_2",
                  "MSE_lepskiboot_bs_1", "MSE_lepskiboot_bs_2",
                  "MSE_CVM_ns_1", "MSE_CVM_ns_2",
                  "MSE_CVMSE_ns_1", "MSE_CVMSE_ns_2",
                  "MSE_lepski_ns_1", "MSE_lepski_ns_2",
                  "MSE_lepskiboot_ns_1", "MSE_lepskiboot_ns_2")


library(ggplot2)
library(patchwork)

graph_1 <- ggplot(df_res1_2, aes(x = c(0.77, 1.58, 2.39, 3.20, 4))) +
  geom_line(aes(y = MSE_CVM_bs_1, colour = 'CVM  - sim 1', linetype = 'CVM  - sim 1'), show.legend = FALSE) +
  geom_line(aes(y = MSE_CVM_bs_2, colour = 'CVM  - sim 2', linetype = 'CVM  - sim 2'), show.legend = FALSE) +
  
  geom_line(aes(y = MSE_CVMSE_bs_1, colour = 'CVMSE  - sim 1', linetype = 'CVMSE  - sim 1'), show.legend = FALSE) +
  geom_line(aes(y = MSE_CVMSE_bs_2, colour = 'CVMSE  - sim 2', linetype = 'CVMSE  - sim 2'), show.legend = FALSE) +
  
  geom_line(aes(y = MSE_lepski_bs_1, colour = 'Lepski  - sim 1', linetype = 'Lepski  - sim 1'), show.legend = FALSE) +
  geom_line(aes(y = MSE_lepski_bs_2, colour = 'Lepski  - sim 2', linetype = 'Lepski  - sim 2'), show.legend = FALSE) +
  
  geom_line(aes(y = MSE_lepskiboot_bs_1, colour = 'Lepski boot  - sim 1', linetype = 'Lepski boot  - sim 1'), show.legend = FALSE) +
  geom_line(aes(y = MSE_lepskiboot_bs_2, colour = 'Lepski boot  - sim 2', linetype = 'Lepski boot  - sim 2'), show.legend = FALSE) +
  
  ylab("MSE") +
  xlab("Interval length") + 
  scale_colour_manual(
    name = "Functions",
    values = c(
      "CVM  - sim 1" = "#E69F00", "CVM  - sim 2" = "#E69F00",
      "CVMSE  - sim 1" = "#56B4E9", "CVMSE  - sim 2" = "#56B4E9",
      "Lepski  - sim 1" = "#661100", "Lepski  - sim 2" = "#661100",
      "Lepski boot  - sim 1" = "#CC79A7", "Lepski boot  - sim 2" = "#CC79A7"
    )
  ) +
  scale_linetype_manual(
    name = "Functions",
    values = c(
      "CVM  - sim 1" = "solid", "CVM  - sim 2" = "dashed",
      "CVMSE  - sim 1" = "solid", "CVMSE  - sim 2" = "dashed",
      "Lepski  - sim 1" = "solid", "Lepski  - sim 2" = "dashed",
      "Lepski boot  - sim 1" = "solid", "Lepski boot  - sim 2" = "dashed"
    )
  ) +
  ggtitle("MSE comparison - Simulation 1 vs 2 - B splines")


graph_2 <- ggplot(df_res1_2, aes(x = c(0.77, 1.58, 2.39, 3.20, 4))) +
  geom_line(aes(y = MSE_CVM_ns_1, colour = 'CVM  - sim 1', linetype = 'CVM  - sim 1'), show.legend = TRUE) +
  geom_line(aes(y = MSE_CVM_ns_2, colour = 'CVM  - sim 2', linetype = 'CVM  - sim 2'), show.legend = TRUE) +
  
  geom_line(aes(y = MSE_CVMSE_ns_1, colour = 'CVMSE  - sim 1', linetype = 'CVMSE  - sim 1'), show.legend = TRUE) +
  geom_line(aes(y = MSE_CVMSE_ns_2, colour = 'CVMSE  - sim 2', linetype = 'CVMSE  - sim 2'), show.legend = TRUE) +
  
  geom_line(aes(y = MSE_lepski_ns_1, colour = 'Lepski  - sim 1', linetype = 'Lepski  - sim 1'), show.legend = TRUE) +
  geom_line(aes(y = MSE_lepski_ns_2, colour = 'Lepski  - sim 2', linetype = 'Lepski  - sim 2'), show.legend = TRUE) +
  
  geom_line(aes(y = MSE_lepskiboot_ns_1, colour = 'Lepski boot  - sim 1', linetype = 'Lepski boot  - sim 1'), show.legend = TRUE) +
  geom_line(aes(y = MSE_lepskiboot_ns_2, colour = 'Lepski boot  - sim 2', linetype = 'Lepski boot  - sim 2'), show.legend = TRUE) +
  
  ylab("MSE") +
  xlab("Interval length") + 
  scale_colour_manual(
    name = "Functions",
    values = c(
      "CVM  - sim 1" = "#E69F00", "CVM  - sim 2" = "#E69F00",
      "CVMSE  - sim 1" = "#56B4E9", "CVMSE  - sim 2" = "#56B4E9",
      "Lepski  - sim 1" = "#661100", "Lepski  - sim 2" = "#661100",
      "Lepski boot  - sim 1" = "#CC79A7", "Lepski boot  - sim 2" = "#CC79A7"
    )
  ) +
  scale_linetype_manual(
    name = "Functions",
    values = c(
      "CVM  - sim 1" = "solid", "CVM  - sim 2" = "dashed",
      "CVMSE  - sim 1" = "solid", "CVMSE  - sim 2" = "dashed",
      "Lepski  - sim 1" = "solid", "Lepski  - sim 2" = "dashed",
      "Lepski boot  - sim 1" = "solid", "Lepski boot  - sim 2" = "dashed"
    )
  ) +
  ggtitle("MSE comparison - Simulation 1 vs 2 - Natural splines")


combined_graph <- graph_1 + graph_2 + plot_layout(ncol = 2)
plot(combined_graph)




df_res1_2 <- data.frame(bias_res_1$bias_CVM_bs, bias_res_2$bias_CVM_bs, 
                        bias_res_1$bias_CVMSE_bs, bias_res_2$bias_CVMSE_bs,
                        bias_res_1$bias_lepski_bs, bias_res_2$bias_lepski_bs,
                        bias_res_1$bias_lepskiboot_bs, bias_res_2$bias_lepskiboot_bs,
                        bias_res_1$bias_CVM_ns, bias_res_2$bias_CVM_ns, 
                        bias_res_1$bias_CVMSE_ns, bias_res_2$bias_CVMSE_ns,
                        bias_res_1$bias_lepski_ns, bias_res_2$bias_lepski_ns,
                        bias_res_1$bias_lepskiboot_ns, bias_res_2$bias_lepskiboot_ns)

colnames(df_res1_2) <- c("bias_CVM_bs_1", "bias_CVM_bs_2",
                         "bias_CVMSE_bs_1", "bias_CVMSE_bs_2",
                         "bias_lepski_bs_1", "bias_lepski_bs_2",
                         "bias_lepskiboot_bs_1", "bias_lepskiboot_bs_2",
                         "bias_CVM_ns_1", "bias_CVM_ns_2",
                         "bias_CVMSE_ns_1", "bias_CVMSE_ns_2",
                         "bias_lepski_ns_1", "bias_lepski_ns_2",
                         "bias_lepskiboot_ns_1", "bias_lepskiboot_ns_2")


library(ggplot2)
library(patchwork)

graph_1 <- ggplot(df_res1_2, aes(x = c(0.77, 1.58, 2.39, 3.20, 4))) +
  geom_line(aes(y = bias_CVM_bs_1, colour = 'CVM  - sim 1', linetype = 'CVM  - sim 1'), show.legend = FALSE) +
  geom_line(aes(y = bias_CVM_bs_2, colour = 'CVM  - sim 2', linetype = 'CVM  - sim 2'), show.legend = FALSE) +
  
  geom_line(aes(y = bias_CVMSE_bs_1, colour = 'CVMSE  - sim 1', linetype = 'CVMSE  - sim 1'), show.legend = FALSE) +
  geom_line(aes(y = bias_CVMSE_bs_2, colour = 'CVMSE  - sim 2', linetype = 'CVMSE  - sim 2'), show.legend = FALSE) +
  
  geom_line(aes(y = bias_lepski_bs_1, colour = 'Lepski  - sim 1', linetype = 'Lepski  - sim 1'), show.legend = FALSE) +
  geom_line(aes(y = bias_lepski_bs_2, colour = 'Lepski  - sim 2', linetype = 'Lepski  - sim 2'), show.legend = FALSE) +
  
  geom_line(aes(y = bias_lepskiboot_bs_1, colour = 'Lepski boot  - sim 1', linetype = 'Lepski boot  - sim 1'), show.legend = FALSE) +
  geom_line(aes(y = bias_lepskiboot_bs_2, colour = 'Lepski boot  - sim 2', linetype = 'Lepski boot  - sim 2'), show.legend = FALSE) +
  
  ylab("bias") +
  xlab("Interval length") + 
  scale_colour_manual(
    name = "Functions",
    values = c(
      "CVM  - sim 1" = "#E69F00", "CVM  - sim 2" = "#E69F00",
      "CVMSE  - sim 1" = "#56B4E9", "CVMSE  - sim 2" = "#56B4E9",
      "Lepski  - sim 1" = "#661100", "Lepski  - sim 2" = "#661100",
      "Lepski boot  - sim 1" = "#CC79A7", "Lepski boot  - sim 2" = "#CC79A7"
    )
  ) +
  scale_linetype_manual(
    name = "Functions",
    values = c(
      "CVM  - sim 1" = "solid", "CVM  - sim 2" = "dashed",
      "CVMSE  - sim 1" = "solid", "CVMSE  - sim 2" = "dashed",
      "Lepski  - sim 1" = "solid", "Lepski  - sim 2" = "dashed",
      "Lepski boot  - sim 1" = "solid", "Lepski boot  - sim 2" = "dashed"
    )
  ) +
  ggtitle("bias comparison - Simulation 1 vs 2 - B splines")


graph_2 <- ggplot(df_res1_2, aes(x = c(0.77, 1.58, 2.39, 3.20, 4))) +
  geom_line(aes(y = bias_CVM_ns_1, colour = 'CVM  - sim 1', linetype = 'CVM  - sim 1'), show.legend = TRUE) +
  geom_line(aes(y = bias_CVM_ns_2, colour = 'CVM  - sim 2', linetype = 'CVM  - sim 2'), show.legend = TRUE) +
  
  geom_line(aes(y = bias_CVMSE_ns_1, colour = 'CVMSE  - sim 1', linetype = 'CVMSE  - sim 1'), show.legend = TRUE) +
  geom_line(aes(y = bias_CVMSE_ns_2, colour = 'CVMSE  - sim 2', linetype = 'CVMSE  - sim 2'), show.legend = TRUE) +
  
  geom_line(aes(y = bias_lepski_ns_1, colour = 'Lepski  - sim 1', linetype = 'Lepski  - sim 1'), show.legend = TRUE) +
  geom_line(aes(y = bias_lepski_ns_2, colour = 'Lepski  - sim 2', linetype = 'Lepski  - sim 2'), show.legend = TRUE) +
  
  geom_line(aes(y = bias_lepskiboot_ns_1, colour = 'Lepski boot  - sim 1', linetype = 'Lepski boot  - sim 1'), show.legend = TRUE) +
  geom_line(aes(y = bias_lepskiboot_ns_2, colour = 'Lepski boot  - sim 2', linetype = 'Lepski boot  - sim 2'), show.legend = TRUE) +
  
  ylab("bias") +
  xlab("Interval length") + 
  scale_colour_manual(
    name = "Functions",
    values = c(
      "CVM  - sim 1" = "#E69F00", "CVM  - sim 2" = "#E69F00",
      "CVMSE  - sim 1" = "#56B4E9", "CVMSE  - sim 2" = "#56B4E9",
      "Lepski  - sim 1" = "#661100", "Lepski  - sim 2" = "#661100",
      "Lepski boot  - sim 1" = "#CC79A7", "Lepski boot  - sim 2" = "#CC79A7"
    )
  ) +
  scale_linetype_manual(
    name = "Functions",
    values = c(
      "CVM  - sim 1" = "solid", "CVM  - sim 2" = "dashed",
      "CVMSE  - sim 1" = "solid", "CVMSE  - sim 2" = "dashed",
      "Lepski  - sim 1" = "solid", "Lepski  - sim 2" = "dashed",
      "Lepski boot  - sim 1" = "solid", "Lepski boot  - sim 2" = "dashed"
    )
  ) +
  ggtitle("bias comparison - Simulation 1 vs 2 - Natural splines")


combined_graph <- graph_1 + graph_2 + plot_layout(ncol = 2)
plot(combined_graph)

