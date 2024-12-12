# data simulation with the cluster
source("source_file_common_all.R")

#simulation 1
source("/softs/R/createCluster.R")
cl <- createCluster()
registerDoParallel(cl)
set.seed(47820)

rhouv <- 0.5
rhozw <- 0.9
case <- 2
n_values <- 1000
data_param = c(n_values, rhouv, rhozw)
n_MC = 2000

simul_all <- vector("list", n_MC)

simul_all <-foreach(n = 1:n_MC, .packages = c("splines", "MASS", "caret", "expm")) %dopar% {
    # Simulate data for each iteration
    simul_inter <- simulate_data_3(data_param, g_sim_3, case)
    simul_all[[n]] <- simul_inter
}

filename <- paste("data_", n_MC, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, ".R" ,sep = "")
save(simul_all, file = filename)


stopCluster(cl)

#simulation 2
source("/softs/R/createCluster.R")
cl <- createCluster()
registerDoParallel(cl)
set.seed(47820)

rhouv <- 0.5
rhozw <- 0.9
case <- 2
n_values <- 2500
data_param = c(n_values, rhouv, rhozw)
n_MC = 2000

simul_all <- vector("list", n_MC)

simul_all <-foreach(n = 1:n_MC, .packages = c("splines", "MASS", "caret", "expm")) %dopar% {
  # Simulate data for each iteration
  simul_inter <- simulate_data_3(data_param, g_sim_3, case)
  simul_all[[n]] <- simul_inter
}

filename <- paste("data_", n_MC, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, ".R" ,sep = "")
save(simul_all, file = filename)

stopCluster(cl)


#simulation 3
source("/softs/R/createCluster.R")
cl <- createCluster()
registerDoParallel(cl)
set.seed(47820)

rhouv <- 0.5
rhozw <- 0.9
case <- 3
n_values <- 1000
data_param = c(n_values, rhouv, rhozw)
n_MC = 2000

simul_all <- vector("list", n_MC)

simul_all <-foreach(n = 1:n_MC, .packages = c("splines", "MASS", "caret", "expm")) %dopar% {
  # Simulate data for each iteration
  simul_inter <- simulate_data_3(data_param, g_sim_3, case)
  simul_all[[n]] <- simul_inter
}

filename <- paste("data_", n_MC, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, ".R" ,sep = "")
save(simul_all, file = filename)

stopCluster(cl)


#simulation 4
source("/softs/R/createCluster.R")
cl <- createCluster()
registerDoParallel(cl)
set.seed(47820)

rhouv <- 0.8
rhozw <- 0.9
case <- 3
n_values <- 1000
data_param = c(n_values, rhouv, rhozw)
n_MC = 2000

simul_all <- vector("list", n_MC)

simul_all <-foreach(n = 1:n_MC, .packages = c("splines", "MASS", "caret", "expm")) %dopar% {
  # Simulate data for each iteration
  simul_inter <- simulate_data_3(data_param, g_sim_3, case)
  simul_all[[n]] <- simul_inter
}

filename <- paste("data_", n_MC, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, ".R" ,sep = "")
save(simul_all, file = filename)

stopCluster(cl)



#simulation 5
source("/softs/R/createCluster.R")
cl <- createCluster()
registerDoParallel(cl)
set.seed(47820)

rhouv <- 0.8
rhozw <- 0.7
case <- 3
n_values <- 1000
data_param = c(n_values, rhouv, rhozw)
n_MC = 2000

simul_all <- vector("list", n_MC)

simul_all <-foreach(n = 1:n_MC, .packages = c("splines", "MASS", "caret", "expm")) %dopar% {
  # Simulate data for each iteration
  simul_inter <- simulate_data_3(data_param, g_sim_3, case)
  simul_all[[n]] <- simul_inter
}

filename <- paste("data_", n_MC, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, ".R" ,sep = "")
save(simul_all, file = filename)

stopCluster(cl)


