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


#### Number of lost simulations, fixed J  ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_441463_complete_allJ_benchmark")
degree = 3
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

name <- c()
number <- c()

for (j in 1:nrow(parameter_combinations)){
  params <- parameter_combinations[j,]
  rhouv <- as.numeric(Rhouv_Rhozw[[params$Rhouv_Rhozw]][1])
  rhozw <- as.numeric(Rhouv_Rhozw[[params$Rhouv_Rhozw]][2])
  case <- as.numeric(params$Case)
  n_val <- as.numeric(params$N)
  data_param = c(n_val, rhouv, rhozw)
  for (i in 1:length(J_val)){
    J <- J_val[i]
    filename <- paste ("MC2000_fixedJ", "_J", J, "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
    load(filename)
    new_MC <- c(filename, new_MC)
    assign(filename, new_MC)
    name <- c(name, filename)
    number <- c(number, length(new_MC$list_gamma))
    rm(new_MC)
  }
}

df_lost <- data.frame(name, number)
df_lost <- df_lost %>% mutate(
  J = str_extract(name, "(?<=_J)\\d+"),
  degree = str_extract(name, "(?<=_degree)\\d+"),
  rhozw = str_extract(name, "(?<=_rhozw)\\d*\\.?\\d+"),
  rhouv = str_extract(name, "(?<=_rhouv)\\d*\\.?\\d+"),
  case = str_extract(name, "(?<=_case)\\d+"),
  n_val = str_extract(name, "(?<=_n)\\d+")
)

df_lost <- df_lost %>% select(-c(name))

install.packages("writexl")
library(writexl)
write_xlsx(df_lost, "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/all_J_benchmark_lost.xlsx")

#### Number of lost simulations, fixed J, NS  ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451423_allJ_NS")
degree = 3
n_MC = 2000

N_values <- c(200, 400, 1000, 2500)
Cases <- c(2, 3)
Rhouv_Rhozw <- list(c(0.5, 0.9), c(0.8, 0.9), c(0.8, 0.7))

parameter_combinations <- expand.grid(
  N = N_values,
  Case = Cases,
  Rhouv_Rhozw = seq_along(Rhouv_Rhozw) # Use indices for combinations
)

J_val <- c(4, 8, 16, 32) 

name <- c()
number <- c()

for (j in 1:nrow(parameter_combinations)){
  params <- parameter_combinations[j,]
  rhouv <- as.numeric(Rhouv_Rhozw[[params$Rhouv_Rhozw]][1])
  rhozw <- as.numeric(Rhouv_Rhozw[[params$Rhouv_Rhozw]][2])
  case <- as.numeric(params$Case)
  n_val <- as.numeric(params$N)
  data_param = c(n_val, rhouv, rhozw)
  for (i in 1:length(J_val)){
    J <- J_val[i]
    filename <- paste ("NS_MC2000_fixedJ", "_J", J, "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
    load(filename)
    new_MC <- c(filename, new_MC)
    assign(filename, new_MC)
    name <- c(name, filename)
    number <- c(number, length(new_MC$list_gamma))
    rm(new_MC)
  }
}

df_lost <- data.frame(name, number)
df_lost <- df_lost %>% mutate(
  J = str_extract(name, "(?<=_J)\\d+"),
  degree = str_extract(name, "(?<=_degree)\\d+"),
  rhozw = str_extract(name, "(?<=_rhozw)\\d*\\.?\\d+"),
  rhouv = str_extract(name, "(?<=_rhouv)\\d*\\.?\\d+"),
  case = str_extract(name, "(?<=_case)\\d+"),
  n_val = str_extract(name, "(?<=_n)\\d+")
)

df_lost <- df_lost %>% select(-c(name))

#install.packages("writexl")
#library(writexl)
write_xlsx(df_lost, "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/all_J_benchmark_lost_NS.xlsx")




#### Matrix of performances, fixed J ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_441463_complete_allJ_benchmark")

#load the data
degree = 3
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

colnames(df) <- c("Name", "M", "supnorm", "Var", "MSE", "bias")
df <- df %>%
  mutate(
    J = str_extract(Name, "(?<=_J)\\d+"),
    degree = str_extract(Name, "(?<=_degree)\\d+"),
    rhozw = str_extract(Name, "(?<=_rhozw)\\d*\\.?\\d+"),
    rhouv = str_extract(Name, "(?<=_rhouv)\\d*\\.?\\d+"),
    case = str_extract(Name, "(?<=_case)\\d+"),
    n_val = str_extract(Name, "(?<=_n)\\d+")
  )


df <- df %>%
  rownames_to_column(var = "Temp") %>%   # Create a temporary column for row names
  select(-Temp)

df <- df%>%select(-1)
df <- df%>%select(n_val, J, case, rhozw, rhouv, degree, M, supnorm, Var, MSE, bias)
df <- data.frame(lapply(df, as.numeric))
df <- df[order(df$n_val, df$J), ]

df <- data.frame(lapply(df, as.character))
df <- subset(df, select=-c(degree))

df_2 <- df %>% filter(case == 2) 
df_3 <- df %>% filter(case == 3)

xtable(df_2)
xtable(df_3)

#install.packages("writexl")
#library(writexl)
write_xlsx(df, "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/perf_all_J.xlsx")

#### Matrix of performances, fixed J, NS ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451423_allJ_NS")

#load the data
degree = 3
n_MC = 2000

N_values <- c(200, 400, 1000, 2500)
Cases <- c(2, 3)
Rhouv_Rhozw <- list(c(0.5, 0.9), c(0.8, 0.9), c(0.8, 0.7))

parameter_combinations <- expand.grid(
  N = N_values,
  Case = Cases,
  Rhouv_Rhozw = seq_along(Rhouv_Rhozw) # Use indices for combinations
)

J_val <- c(4, 8, 16, 32)

for (j in 1:nrow(parameter_combinations)){
  params <- parameter_combinations[j,]
  rhouv <- as.numeric(Rhouv_Rhozw[[params$Rhouv_Rhozw]][1])
  rhozw <- as.numeric(Rhouv_Rhozw[[params$Rhouv_Rhozw]][2])
  case <- as.numeric(params$Case)
  n_val <- as.numeric(params$N)
  data_param = c(n_val, rhouv, rhozw)
  for (i in 1:length(J_val)){
    J <- J_val[i]
    filename_perf <- paste ("NS_perf2000_fixedJ", "_J", J, "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
    load(filename_perf)
    perf_MC <- c(filename_perf, perf_MC)
    assign(filename_perf, perf_MC)
  }
}

#create a df with aggregated performances
vector_names <- ls(pattern = "^NS_perf2000")
vector_list <- mget(vector_names)
df <- as.data.frame(do.call(rbind, vector_list))

colnames(df) <- c("Name", "M", "supnorm", "Var", "MSE", "bias")
df <- df %>%
  mutate(
    J = str_extract(Name, "(?<=_J)\\d+"),
    degree = str_extract(Name, "(?<=_degree)\\d+"),
    rhozw = str_extract(Name, "(?<=_rhozw)\\d*\\.?\\d+"),
    rhouv = str_extract(Name, "(?<=_rhouv)\\d*\\.?\\d+"),
    case = str_extract(Name, "(?<=_case)\\d+"),
    n_val = str_extract(Name, "(?<=_n)\\d+")
  )


df <- df %>%
  rownames_to_column(var = "Temp") %>%   # Create a temporary column for row names
  select(-Temp)

df <- df%>%select(-1) #remove the name
df <- df%>%select(n_val, J, case, rhozw, rhouv, degree, M, supnorm, Var, MSE, bias) #order the columns
df <- data.frame(lapply(df, as.numeric))
df <- df[order(df$n_val, df$J), ]

df <- data.frame(lapply(df, as.character))
df <- subset(df, select=-c(degree))

df_2 <- df %>% filter(case == 2) 
df_3 <- df %>% filter(case == 3)

xtable(df_2)
xtable(df_3)


#install.packages("writexl")
#library(writexl)
write_xlsx(df, "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/perf_all_J_NS.xlsx")





#### Performance CVM ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451416_selectionJ_CVM_complete")

#load("perf2000_CV_M_degree3_rhozw0.9_rhouv0.5_case2_n1000.R")
#load the data
degree = 3
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
    filename_perf <- paste ("perf2000_CV_M", "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
    load(filename_perf)
    perf_CV_M <- c(filename_perf, perf_CV_M)
    assign(filename_perf, perf_CV_M)
  }
}

#create a df with aggregated performances
vector_names <- ls(pattern = "^perf2000")
vector_list <- mget(vector_names)
df <- as.data.frame(do.call(rbind, vector_list))

colnames(df) <- c("Name", "M", "supnorm", "Var", "MSE", "bias")
df <- df %>%
  mutate(
    degree = str_extract(Name, "(?<=_degree)\\d+"),
    rhozw = str_extract(Name, "(?<=_rhozw)\\d*\\.?\\d+"),
    rhouv = str_extract(Name, "(?<=_rhouv)\\d*\\.?\\d+"),
    case = str_extract(Name, "(?<=_case)\\d+"),
    n_val = str_extract(Name, "(?<=_n)\\d+")
  )


df <- df %>%
  rownames_to_column(var = "Temp") %>%   # Create a temporary column for row names
  select(-Temp)

df <- df%>%select(-1)
df <- df%>%select(n_val, case, rhozw, rhouv, degree, M, supnorm, Var, MSE, bias)
df <- data.frame(lapply(df, as.numeric))
df <- df[order(df$n_val, df$case), ]

df <- data.frame(lapply(df, as.character))
df <- subset(df, select=-c(degree))

df_2 <- df %>% filter(case == 2) 
df_3 <- df %>% filter(case == 3)

xtable(df_2)
xtable(df_3)

#install.packages("writexl")
#library(writexl)
write_xlsx(df, "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/perf_CV_M.xlsx")



#### Performance CVMSE ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451417_selectionJ_CVMSE_complete")

#load("perf2000_CV_M_degree3_rhozw0.9_rhouv0.5_case2_n1000.R")
#load the data
degree = 3
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
    filename_perf <- paste ("perf2000_CV_MSE", "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
    load(filename_perf)
    perf_CV_MSE <- c(filename_perf, perf_CV_MSE)
    assign(filename_perf, perf_CV_MSE)
  }
}

#create a df with aggregated performances
vector_names <- ls(pattern = "^perf2000")
vector_list <- mget(vector_names)
df <- as.data.frame(do.call(rbind, vector_list))

colnames(df) <- c("Name", "M", "supnorm", "Var", "MSE", "bias")
df <- df %>%
  mutate(
    degree = str_extract(Name, "(?<=_degree)\\d+"),
    rhozw = str_extract(Name, "(?<=_rhozw)\\d*\\.?\\d+"),
    rhouv = str_extract(Name, "(?<=_rhouv)\\d*\\.?\\d+"),
    case = str_extract(Name, "(?<=_case)\\d+"),
    n_val = str_extract(Name, "(?<=_n)\\d+")
  )


df <- df %>%
  rownames_to_column(var = "Temp") %>%   # Create a temporary column for row names
  select(-Temp)

df <- df%>%select(-1)
df <- df%>%select(n_val, case, rhozw, rhouv, degree, M, supnorm, Var, MSE, bias)
df <- data.frame(lapply(df, as.numeric))
df <- df[order(df$n_val, df$case), ]

df <- data.frame(lapply(df, as.character))
df <- subset(df, select=-c(degree))

df_2 <- df %>% filter(case == 2) 
df_3 <- df %>% filter(case == 3)

xtable(df_2)
xtable(df_3)

#install.packages("writexl")
#library(writexl)
write_xlsx(df, "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/perf_CV_MSE.xlsx")




#### Test ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451416_selectionJ_CVM_complete")
load("opt_CV_M_2000_degree3_rhozw0.7_rhouv0.8_case2_n2500.R")
length(new_opt_CV_M$list_gamma)

plot_mean_true(new_opt_CV_M, seq(-2, 2, length.out = 100), 'none', 3, 0.7, 0.8,2)





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
