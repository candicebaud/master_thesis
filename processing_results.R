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
    rm(list = filename)
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
    rm(list = filename)
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





#### Lost simulations, CVM ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451416_selectionJ_CVM_complete")
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
  filename <- paste("opt_CV_M_2000", "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
  load(filename)
  new_opt_CV_M <- c(filename, new_opt_CV_M)
  assign(filename, new_opt_CV_M)
  name <- c(name, filename)
  number <- c(number, length(new_opt_CV_M$list_gamma))
  rm(list = filename)
}

df_lost <- data.frame(name, number)
df_lost <- df_lost %>% mutate(
  degree = str_extract(name, "(?<=_degree)\\d+"),
  rhozw = str_extract(name, "(?<=_rhozw)\\d*\\.?\\d+"),
  rhouv = str_extract(name, "(?<=_rhouv)\\d*\\.?\\d+"),
  case = str_extract(name, "(?<=_case)\\d+"),
  n_val = str_extract(name, "(?<=_n)\\d+"))

df_lost <- df_lost %>% select(-c(name))

#install.packages("writexl")
library(writexl)
write_xlsx(df_lost, "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/CVM_lost.xlsx")

#### Lost simulations, CVM NS ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451438_selectionJ_CVM_NS_complete")
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
  filename <- paste("NS_opt_CV_M_2000", "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
  load(filename)
  new_opt_CV_M <- c(filename, new_opt_CV_M)
  assign(filename, new_opt_CV_M)
  name <- c(name, filename)
  number <- c(number, length(new_opt_CV_M$list_gamma))
  rm(list = filename)
}

df_lost <- data.frame(name, number)
df_lost <- df_lost %>% mutate(
  degree = str_extract(name, "(?<=_degree)\\d+"),
  rhozw = str_extract(name, "(?<=_rhozw)\\d*\\.?\\d+"),
  rhouv = str_extract(name, "(?<=_rhouv)\\d*\\.?\\d+"),
  case = str_extract(name, "(?<=_case)\\d+"),
  n_val = str_extract(name, "(?<=_n)\\d+"))

df_lost <- df_lost %>% select(-c(name))

#install.packages("writexl")
library(writexl)
write_xlsx(df_lost, "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/CVM_NS_lost.xlsx")



#### Lost simulations, CVMSE ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451417_selectionJ_CVMSE_complete")
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
  filename <- paste("opt_CV_MSE_2000", "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
  load(filename)
  new_opt_CV_MSE <- c(filename, new_opt_CV_MSE)
  assign(filename, new_opt_CV_MSE)
  name <- c(name, filename)
  number <- c(number, length(new_opt_CV_MSE$list_gamma))
  rm(list = filename)
}

df_lost <- data.frame(name, number)
df_lost <- df_lost %>% mutate(
  degree = str_extract(name, "(?<=_degree)\\d+"),
  rhozw = str_extract(name, "(?<=_rhozw)\\d*\\.?\\d+"),
  rhouv = str_extract(name, "(?<=_rhouv)\\d*\\.?\\d+"),
  case = str_extract(name, "(?<=_case)\\d+"),
  n_val = str_extract(name, "(?<=_n)\\d+"))

df_lost <- df_lost %>% select(-c(name))

#install.packages("writexl")
library(writexl)
write_xlsx(df_lost, "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/CVMSE_lost.xlsx")


#### Lost simulations, CVMSE NS ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451439_selectionJ_CVMSE_NS_complete")
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
  filename <- paste("NS_opt_CV_MSE_2000", "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
  load(filename)
  new_opt_CV_MSE <- c(filename, new_opt_CV_MSE)
  assign(filename, new_opt_CV_MSE)
  name <- c(name, filename)
  number <- c(number, length(new_opt_CV_MSE$list_gamma))
  rm(list = filename)
}

df_lost <- data.frame(name, number)
df_lost <- df_lost %>% mutate(
  degree = str_extract(name, "(?<=_degree)\\d+"),
  rhozw = str_extract(name, "(?<=_rhozw)\\d*\\.?\\d+"),
  rhouv = str_extract(name, "(?<=_rhouv)\\d*\\.?\\d+"),
  case = str_extract(name, "(?<=_case)\\d+"),
  n_val = str_extract(name, "(?<=_n)\\d+"))

df_lost <- df_lost %>% select(-c(name))

#install.packages("writexl")
library(writexl)
write_xlsx(df_lost, "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/CVMSE_NS_lost.xlsx")



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





#### Performance CVM NS ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451438_selectionJ_CVM_NS_complete")

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
    filename_perf <- paste ("NS_perf2000_CV_M", "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
    load(filename_perf)
    perf_CV_M <- c(filename_perf, perf_CV_M)
    assign(filename_perf, perf_CV_M)
  }
}

#create a df with aggregated performances
vector_names <- ls(pattern = "^NS_perf2000")
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
write_xlsx(df, "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/perf_CV_M_NS.xlsx")





#### Performance CVMSE NS ####
setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/simu_451439_selectionJ_CVMSE_NS_complete")

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
    filename_perf <- paste("NS_perf2000_CV_MSE", "_degree", degree, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_val, ".R" ,sep = "")
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
write_xlsx(df, "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/perf_CV_MSE_NS.xlsx")





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