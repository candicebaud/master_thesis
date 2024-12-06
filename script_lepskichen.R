#### Lepski Chen ####
V_hat_J <- function(J, n, s_J_hat){
  return (sqrt(log(n)/n)/s_J_hat) #the one in our case 
}

#difference_norm <- function(gamma_J, gamma_J_prime, J, J_prime, Z, degree,create_P, x_grid){
#  P_x_grid_J <- create_P(x_grid, Z, J,degree)
#  P_x_grid_J_prime <- create_P(x_grid, Z, J_prime,degree)
#  g_hat_J_x_grid <- P_x_grid_J%*%gamma_J
#  g_hat_J_prime_x_grid <- P_x_grid_J_prime%*%gamma_J_prime
#  m_m = max(Z) - min(Z)
#  n = length(Z)
#  
#  rm(P_x_grid_J, P_x_grid_J_prime)
#  gc()
  
#  return(m_m*sum((g_hat_J_prime_x_grid - g_hat_J_x_grid)^{2})/n)}

lepski_chen <- function(I_hat, c_0, m_m, list_g_hat_J, W, Z, Y, degree, valid_dim, create_P){
  n = length(Z)
  
  length_I_hat = length(I_hat)
  
  J_index = 1
  J = I_hat[J_index]
  condition = FALSE
  while ((J_index<length_I_hat && condition == FALSE)){
    s_J <- calcul_s_J(J, W, Z, Y, degree, create_P)
    V_J = V_hat_J(J, n, s_J)
    J_prime_index = J_index + 1 
    J_prime = I_hat[J_prime_index]
    while ((J_prime_index < length_I_hat && condition == FALSE)){
      s_J_prime <- calcul_s_J(J_prime, W, Z, Y, degree, create_P)
      V_J_prime = V_hat_J(J, n, s_J_prime)
      a = m_m*sum((list_g_hat_J[[J_prime_index]] - list_g_hat_J[[J_index]])^{2})/n
      b = c_0*(V_J + V_J_prime)
      if(a<b){
        condition = TRUE}
      else{
        condition = FALSE}
      J_prime_index = J_prime_index + 1
      J_prime = I_hat[J_prime_index]}
    J_index = J_index + 1
    J = I_hat[J_index]
  }
  rm(condition, J_index)
  gc()
  J_opt = J
  return(J_opt)
}


#MC_lepski <-  function(n_MC, x_eval, degree, valid_dim, g_0, case, data_param){
#  c_0 = 10
#  list_W <- list()
#  list_Y <- list()
#  list_Z <- list()
  
#  list_J_opt <- rep(0, n_MC)
#  list_gamma_opt <- list()
#  list_g_hat_on_x <- list()
  
#  for (n in 1:n_MC){
    #generate data 
#    simul <- simulate_data_3(data_param, g_0, case)
#    W <- simul$W
#    Y <- simul$Y
#    Z <- simul$Z
    
#    list_W[[n]] <- W
#    list_Y[[n]] <- Y
#    list_Z[[n]] <- Z
    
    #compute the optimal J 
#    J_opt <- lepski_chen(c_0, W, Z, Y, degree, valid_dim)
#    list_J_opt[n] <- J_opt
    
    #compute the gamma_J_opt
#    gamma_hat_J_opt <- estimation_gamma(J_opt, W, Z, Y, degree)
#    list_gamma_opt[[n]] <- gamma_hat_J_opt
    
    #compute the function on the grid
#    basis <- create_P(x_eval, Z, J_opt, degree)
#    g_hat_on_x <- basis%*%gamma_hat_J_opt
#    list_g_hat_on_x[[n]] <- g_hat_on_x}
  
#  g_0_on_x <- g_0(x_eval, case)
  
#  return(list(list_J_opt = list_J_opt, list_gamma = list_gamma_opt, 
#              list_g_hat_on_x = list_g_hat_on_x, g_0_on_x = g_0_on_x,
#              list_W = list_W, list_Y = list_Y, list_Z = list_Z))}

