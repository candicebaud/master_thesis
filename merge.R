#splitted simulations

setwd("C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/simulation_results/3/nMC500_n1000_case2_rho1")

results <- list()

n_MC = 500
degree = 3
p_train = 0.5
n_boot = 100
rhozw = 0.9
rhouv = 0.5
case = 2
n_values = 1000

J_bs <- c(5, 7, 11, 19, 35)
J_ns <- c(3, 5, 9, 17, 33)

n_eval = 100


for (i in 1:n_MC) {
  # Construct the file name
  filename <- paste("opt_", n_MC, "_degree", degree, "_ptrain", p_train,
                    "_nboot", n_boot, "_rhozw", rhozw, "_rhouv", rhouv,
                    "_case", case, "_n", n_values, "_simu", i, ".R", sep = "")
  if (file.exists(filename)) {
    # Load the file
    loaded_objects <- load(filename)
    
    # Assuming the loaded object contains a single variable, retrieve it
    results[[i]] <- get(loaded_objects)
  } else {
    message(paste("File not found:", filename))
    results[[i]] <- NULL
  }
}

filename = paste("opt_", n_MC , "_degree", degree, "_ptrain", p_train, "_nboot", n_boot, "_rhozw" , rhozw,"_rhouv", rhouv , "_case", case, "_n", n_values, ".R" ,sep = "")
save(results, file = filename)