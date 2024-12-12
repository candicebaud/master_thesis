#### Lepski bootstrap ####
calcul_hat_g_on_Z <- function(j, Z, gamma_hat, a, b, Z_to_evaluate, degree, create_P) {
  # Initialize the result vector
  g_hat_on_x <- numeric(length(Z_to_evaluate))
  # Compute the basis matrix
  basis <- create_P(Z_to_evaluate, Z, j, degree)
  # Compute g_hat_on_Z
  if (length(gamma_hat) > 1) {
    g_hat_on_Z <- basis %*% gamma_hat
  } else {
    g_hat_on_Z <- gamma_hat * sum(basis)
  }
  # Clean up large objects
  rm(basis)
  gc()
  
  return(g_hat_on_Z)
}

create_matrix_U <- function(u_j, u_j_prime){
  n = length(u_j)
  res_mat <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n){
    res_mat[i,i] = u_j[i]*u_j_prime[i]
  }
  return(res_mat)
}


T_stat_D <- function(x_grid, Y, Z, W, degree, create_P, I_hat, w_vector, list_M_boot, matrix_sigma_all_j_all_x, matrix_all_u, matrix_sigma_all_j_j2_all_x){
  I_hat = sort(I_hat)
  length_I_hat = length(I_hat)
  length_x_grid = length(x_grid)
  
  diff_matrix_allx <- replicate(length_x_grid, matrix(0, length_I_hat, length_I_hat), simplify = FALSE)#matrice des diff de D_j
  matrix_ratio_allx <-  replicate(length_x_grid, matrix(0, length_I_hat, length_I_hat), simplify = FALSE)
  
  for (j_index in 1:length(I_hat)){
    j = I_hat[j_index]
    u_hat_j = matrix_all_u[j_index,] #chaque ligne correspond au vecteur des erreurs de prédiction pour J = valeur de la ligne
    u_star_j = u_hat_j*w_vector #on calcule le vecteur bootstrap
    
    P_j_on_x <- create_P(x_grid, Z, j, degree)
    D_j_star = P_j_on_x%*%list_M_boot[[j_index]]%*%u_star_j #vecteur de taille nx 
    
    for (j_prime_index in j_index:length(I_hat)){
      j_prime = I_hat[j_prime_index]
      
      #compute u_star_j_prime
      u_hat_j_prime = matrix_all_u[j_prime_index,]
      u_star_j_prime = u_hat_j_prime*w_vector
      
      #compute D_j_prime sur les x
      P_j_prime_on_x <- create_P(x_grid, Z, j_prime, degree)
      D_j_prime_star = P_j_prime_on_x%*%list_M_boot[[j_prime_index]]%*%u_star_j_prime #vecteur de taille x
      
      #Stockage des résultats
      for (x_index in 1:length(x_grid)){
        diff_matrix_allx[[x_index]][j_index, j_prime_index] = abs(D_j_prime_star[x_index] - D_j_star[x_index]) 
        
        sigma_j_x = matrix_sigma_all_j_all_x[j_index, x_index]
        sigma_j_prime_x = matrix_sigma_all_j_all_x[j_prime_index, x_index]
        sigma_j_j_prime_x = matrix_sigma_all_j_j2_all_x[[x_index]][j_index, j_prime_index]
        sigma_val_x = sigma_j_x + sigma_j_prime_x - 2*sigma_j_j_prime_x
        
        if (sigma_val_x !=0){ #attention on a des sigma = 0 (bizarre) et donc pb
          matrix_ratio_allx[[x_index]][j_index, j_prime_index] = diff_matrix_allx[[x_index]][j_index,j_prime_index]/sqrt(sigma_val_x)
        }}}}
  
  vec_max <- rep(0, length(x_grid))
  for (x_index in 1:length(x_grid)){
    vec_max[x_index] = max(matrix_ratio_allx[[x_index]])
  }
  m = max(vec_max)
  
  rm(I_hat, diff_matrix_allx, matrix_ratio_allx, vec_max, u_hat_j, u_star_j, P_j_on_x, D_j_star, 
     u_hat_j_prime, u_star_j_prime, P_j_prime_on_x, D_j_prime_star)
  gc()
  return(m)}



compute_J_hat <- function(theta_star, I_hat, x_grid, degree,create_P, Z, list_g_hat_J, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x){
  #j_index = 1
  #J = I_hat[j_index]
  #ratio = calcul_ratio(J, j_index, I_hat, x_grid, degree, create_P, Z, list_g_hat_J, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x)
  
  vec_ratio <- numeric(length(I_hat))
  for (j_index in 1:length(I_hat)){
    vec_ratio[j_index] <- calcul_ratio(I_hat[j_index], j_index, I_hat, x_grid, degree, create_P, Z, 
                                       list_g_hat_J, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x)
  }
  
  threshold <- 10 * theta_star
  J <- vector[which(vec_ratio > threshold)[1]]

  # while ((ratio > 10*theta_star && j_index<= length(I_hat))){
  #   j_index <- j_index + 1
  #   if (j_index > length(I_hat)) break
  #   J <- I_hat[j_index]
  #   ratio <- calcul_ratio(J, j_index, I_hat, x_grid, degree, create_P, Z, 
  #                         list_g_hat_J, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x)
  # }
  
  rm(vec_ratio, treshold)
  gc()
  return(J)
}

calcul_ratio <- function(J, J_index, I_hat, x_grid, degree, create_P, Z, list_g_hat_J, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x){
  length_I_hat = length(I_hat)
  res_mat <- matrix(0, length_I_hat, length(x_grid))
  
  #P_j_on_x <- create_P(x_grid, Z, J, degree) #on crée la base pour J donné
  #if (J>1){
  #  g_hat_J_x <- P_j_on_x%*%list_gamma[[J_index]]
  #}else{
  #  g_hat_J_x <- list_gamma[[J_index]]*P_j_on_x}
  
  g_hat_J_x <- list_g_hat_J[[J_index]]
  sigma_J_vect <- matrix_sigma_all_j_all_x[J_index,]
  
  for (J_prime_index in J_index:length(I_hat)){
    J_prime <- I_hat[J_prime_index]
    #P_j_prime_on_x <- create_P(x_grid, Z, J_prime, degree)
    #if (J>1){
    #  g_hat_J_prime_x <- P_j_prime_on_x%*%list_gamma[[J_prime_index]]
    #}else{
    #  g_hat_J_prime_x <- list_gamma[[J_prime_index]]*P_j_prime_on_x}
    
    g_hat_J_prime_x <- list_g_hat_J[[J_prime_index]]
    sigma_J_prime_vect <- matrix_sigma_all_j_all_x[J_prime_index,]
    
    for (x_index in 1:length(x_grid)){
      sigma_J_J_prime <- abs(
        sigma_J_vect[x_index]^2 + 
          sigma_J_prime_vect[x_index]^2 - 
          2 * matrix_sigma_all_j_j2_all_x[[x_index]][J_index, J_prime_index]
      )
      denom = sqrt(sigma_J_J_prime)
      if (denom!=0){
        res_mat[J_prime_index, x_index] = abs(g_hat_J_prime_x[x_index] - g_hat_J_x[x_index])/(denom)
      }}
  }
  
  m = max(res_mat)
  
  rm(sigma_J_vect, res_mat, g_hat_J_x)
  gc()
  return(m)

}


lepski_bootstrap <- function(I_hat, n_boot, list_gamma, list_M_boot, list_g_hat_J, valid_dim,
                             x_grid, W, Z, Y, degree, create_P){#attention : ici x_grid juste pour calculer les sup norm
  
  length_I_hat = length(I_hat)
  n = length(Z)
  
  #draw the w_i
  W_boot_matrix <- matrix(rnorm(n_boot * n, mean = 0, sd = 1), nrow = n_boot, ncol = n)
  
  #compute the sigma_J(x) for all values of J and all values of x
  matrix_sigma_all_j_all_x <- matrix(0, nrow = length_I_hat, ncol = length(x_grid))
  matrix_all_u <- matrix(0, nrow = length_I_hat, ncol = n) #n number of data points (ie length(Z))
  
  for (j_index in 1:length_I_hat){
    j = I_hat[j_index]
    h_hat_j_on_Z = calcul_hat_g_on_Z(j, Z, list_gamma[[j_index]], a, b, Z, degree, create_P)
    u_hat_j = Y - h_hat_j_on_Z 
    matrix_all_u[j_index,] = u_hat_j
    
    P_j_on_x <- create_P(x_grid, Z, j, degree) #create for all x_grid
    U_j_j <- create_matrix_U(u_hat_j, u_hat_j) 
    middle_matrix <- list_M_boot[[j_index]]%*%U_j_j%*%t(list_M_boot[[j_index]]) #ok
    
    vec_sigma_j <- rep(0, length(x_grid))
    if (length(middle_matrix)>1){
      for (x_index in 1:length(x_grid)){
        vec_sigma_j[x_index] <- P_j_on_x[x_index,]%*%(middle_matrix)%*%P_j_on_x[x_index,]
      }
    }
    else{
      for (x_index in 1:length(x_grid)){
        vec_sigma_j[x_index] = middle_matrix*P_j_on_x[x_index]^{2}
      }
    }
    matrix_sigma_all_j_all_x[j_index,] = vec_sigma_j #ligne j prend les valeurs des sigma_j(x) pour toutes les valeurs de x_grid
  }
  
  #compute the sigma_J,J_2(x) for all values of x
  matrix_sigma_all_j_j2_all_x <- replicate(length(x_grid), matrix(0, length_I_hat, length_I_hat), simplify = FALSE)
  for (j_index in 1:length_I_hat){
    j = I_hat[j_index]
    P_j_on_x <- create_P(x_grid, Z, j, degree)
    
    for (j_prime_index in (j_index):length_I_hat){
      j_prime = I_hat[j_prime_index]
      
      #calcul de sigma_j,jprime(x) pour tous les x_grid et ensuite on met dans les matrices
      U_j_jp <- create_matrix_U(matrix_all_u[j_index,], matrix_all_u[j_prime_index,])
      P_j_prime_on_x <- create_P(x_grid, Z,j_prime,degree)
      
      middle_matrix <- list_M_boot[[j_index]]%*%U_j_jp%*%t(list_M_boot[[j_prime_index]])
      
      vec_sigma_j_jprime <- numeric(length(x_grid))
      
      if (j == 1 && j_prime == 1){
        vec_sigma_j_jprime <- rep(0, length(x_grid))
        for (x_index in 1:length(x_grid)){
          vec_sigma_j_jprime[x_index] = middle_matrix*P_j_on_x[x_index]*P_j_prime_on_x[x_index]
          matrix_sigma_all_j_j2_all_x[[x_index]][j_index, j_prime_index] = vec_sigma_j_jprime[x_index]
        }}
      if (j == 1 && j_prime >1){
        vec_sigma_j_jprime <- rep(0, length(x_grid))
        for (x_index in 1:length(x_grid)){
          vec_sigma_j_jprime[x_index] = P_j_on_x[x_index]*(middle_matrix%*%P_j_prime_on_x[x_index,])
          matrix_sigma_all_j_j2_all_x[[x_index]][j_index, j_prime_index] = vec_sigma_j_jprime[x_index]
        }}
      if (j>1 && j_prime == 1){
        vec_sigma_j_jprime <- rep(0, length(x_grid))
        for (x_index in 1:length(x_grid)){
          vec_sigma_j_jprime[x_index] = (P_j_on_x[x_index,]%*%middle_matrix)*P_j_prime_on_x[x_index]
          matrix_sigma_all_j_j2_all_x[[x_index]][j_index, j_prime_index] = vec_sigma_j_jprime[x_index]
        }}
      if (j>1 && j_prime >1){
        vec_sigma_j_jprime <- rep(0, length(x_grid))
        for (x_index in 1:length(x_grid)){
          vec_sigma_j_jprime[x_index] = P_j_on_x[x_index,]%*%(middle_matrix)%*%P_j_prime_on_x[x_index,]
          matrix_sigma_all_j_j2_all_x[[x_index]][j_index, j_prime_index] = vec_sigma_j_jprime[x_index]
        }}
    }}
  
  
  #compute the T stat for each draw of w
  T_i = rep(0, n_boot)
  for (i in 1:n_boot){
    T_i[i] = T_stat_D(x_grid, Y, Z, W, degree, create_P, I_hat, W_boot_matrix[i,], list_M_boot, matrix_sigma_all_j_all_x, matrix_all_u, matrix_sigma_all_j_j2_all_x) 
  }
  J_max = I_hat[length(I_hat)]
  alpha = min(0.5, 1/J_max)
  theta = quantile(T_i, probs = 1 - alpha)
  
  J_n_hat = I_hat[length(I_hat)-1] #on prend l'avant dernier pour être le plus petit strict
  
  rm(list_M_boot, W_boot_matrix, matrix_all_u, T_i)
  gc()
  J_hat = compute_J_hat(theta,I_hat, x_grid, degree, create_P, Z, list_g_hat_J, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x)
  
  rm(list_gamma,matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x, alpha, theta)
  gc()
  J_tilde = min(J_hat, J_n_hat)
  
  
  return(J_tilde)
}



