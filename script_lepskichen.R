#### Lepski Chen ####
compute_J_max <- function(W, Z, Y, degree, valid_dim){
  J_init = degree - 1 + 2
  J_next = degree - 1 + 2^{2}
  n = length(Z)
  s_1 = 1/calcul_s_J(J_init, W, Z, Y, degree)
  s_2 = 1/calcul_s_J(J_next, W, Z, Y, degree)
  prev_ratio = s_1/sqrt(n) 
  new_ratio = s_2/sqrt(n)
  
  J = J_next
  s_J_hat = s_2
  if ((prev_ratio <= 10 && new_ratio > 10 && s_J_hat !=999)) {
    rm(J_init, J_next, n, s_1, s_2, prev_ratio, new_ratio)
    gc()
    return(J)
  } else {
    while (!(prev_ratio <= 10 && new_ratio > 10 && s_J_hat != 999)){
      J = J + 1
      
      # Look for a valid J
      while (!valid_dim(J, degree)) {  
        J = J + 1}
      
      # Calculate s_hat_J and update ratios
      s_hat_J = calcul_s_J(J, W, Z, Y, degree)
      prev_ratio = new_ratio
      new_ratio = s_hat_J / sqrt(n)
    }
    rm(J_init, J_next, n, s_1, s_2, prev_ratio, new_ratio, s_hat_J)
    return(J)}}

calcul_s_J <- function(J, W, Z, Y, degree){
  n = length(W)
  
  Omega <- create_W(W)
  Phi <- create_dyadic_P_splines(Z, Z, J, degree)
  tPhi_Phi <- t(Phi)%*%Phi
  
  bool = verify_solve_is_working(tPhi_Phi)
  
  if (bool == TRUE){
    inverse_tPhi_Phi <- solve(tPhi_Phi)
    sqrt_inverse_tPhi_Phi <- sqrtm(inverse_tPhi_Phi)
    res_mat <- n*sqrt_inverse_tPhi_Phi%*%t(Phi)%*%Omega%*%Phi%*%sqrt_inverse_tPhi_Phi
    
    sing_values <- svd(res_mat)$d
    
    rm(Omega, Phi, tPhi_Phi, inverse_tPhi_Phi, sqrt_inverse_tPhi_Phi, res_mat)
    gc()
    return(min(sing_values))
  }else{
    rm(Omega, tPhi_Phi, Phi)
    gc()
    return(999)
  }}


valid_dim_b_splines <- function(J, degree){ #attention modifier pour NS
  if (J-degree+1 <= 1){
    return(FALSE)
  }
  else{
    l = log(J - degree +1)/log(2)
    if (as.integer(l)!=l){
      return(FALSE)
    }else{return(TRUE)}}}

V_hat_J <- function(J, n, s_J_hat){
  return (sqrt(log(n)/n)/s_J_hat) #the one in our case 
}

difference_norm <- function(gamma_J, gamma_J_prime, J, J_prime, Z, degree, x_grid){
  P_x_grid_J <- create_dyadic_P_splines(x_grid, Z, J,degree)
  P_x_grid_J_prime <- create_dyadic_P_splines(x_grid, Z, J_prime,degree)
  g_hat_J_x_grid <- P_x_grid_J%*%gamma_J
  g_hat_J_prime_x_grid <- P_x_grid_J_prime%*%gamma_J_prime
  m_m = max(Z) - min(Z)
  n = length(Z)
  
  rm(P_x_grid_J, P_x_grid_J_prime)
  gc()
  
  return(m_m*sum((g_hat_J_prime_x_grid - g_hat_J_x_grid)^{2})/n)
}

lepski_chen <- function(c_0,W,Z,Y,degree, valid_dim){
  n = length(Z)
  x_grid = seq(min(Z), max(Z), length.out = 100)
  
  J_max = compute_J_max(W, Z, Y, degree, valid_dim)
  I_hat = seq(as.integer(0.1 * (log(J_max)^2)), as.integer(J_max), by = 1) #on commence au moins à 2 parce que 1 c'est pas une bonne solution
  I_hat = sort(I_hat[sapply(I_hat,valid_dim, degree)]) #select only the valid dimensions
  length_I_hat = length(I_hat)
  
  J_index = 1
  J = I_hat[J_index]
  condition = FALSE
  while ((J_index<length(I_hat) && condition == FALSE)){
    gamma_J <- estimation_gamma(J,W,Z,Y, degree)
    s_J <- calcul_s_J(J, W, Z, Y, degree)
    V_J = V_hat_J(J, n, s_J)
    J_prime_index = J_index + 1 
    J_prime = I_hat[J_prime_index]
    while ((J_prime_index < length(I_hat) && condition == FALSE)){
      gamma_J_prime <- estimation_gamma(J_prime, W, Z, Y, degree)
      s_J_prime <- calcul_s_J(J_prime, W, Z, Y, degree)
      V_J_prime = V_hat_J(J, n, s_J_prime)
      a = difference_norm(gamma_J, gamma_J_prime, J, J_prime, Z, degree, x_grid)
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
  rm(J_max, condition, s_J, V_J, J_prime_index, J_prime,
     J_index)
  gc()
  J_opt = J
  return(J_opt)
}


MC_lepski <-  function(n_MC, x_eval, degree, valid_dim, g_0, case, data_param){
  c_0 = 10
  list_W <- list()
  list_Y <- list()
  list_Z <- list()
  
  list_J_opt <- rep(0, n_MC)
  list_gamma_opt <- list()
  list_g_hat_on_x <- list()
  
  for (n in 1:n_MC){
    #generate data 
    simul <- simulate_data_3(data_param, g_0, case)
    W <- simul$W
    Y <- simul$Y
    Z <- simul$Z
    
    list_W[[n]] <- W
    list_Y[[n]] <- Y
    list_Z[[n]] <- Z
    
    #compute the optimal J 
    J_opt <- lepski_chen(c_0, W, Z, Y, degree, valid_dim)
    list_J_opt[n] <- J_opt
    
    #compute the gamma_J_opt
    gamma_hat_J_opt <- estimation_gamma(J_opt, W, Z, Y, degree)
    list_gamma_opt[[n]] <- gamma_hat_J_opt
    
    #compute the function on the grid
    basis <- create_dyadic_P_splines(x_eval, Z, J_opt, degree)
    g_hat_on_x <- basis%*%gamma_hat_J_opt
    list_g_hat_on_x[[n]] <- g_hat_on_x
  }
  
  g_0_on_x <- g_0(x_eval, case)
  
  return(list(list_J_opt = list_J_opt, list_gamma = list_gamma_opt, 
              list_g_hat_on_x = list_g_hat_on_x, g_0_on_x = g_0_on_x,
              list_W = list_W, list_Y = list_Y, list_Z = list_Z))
}

