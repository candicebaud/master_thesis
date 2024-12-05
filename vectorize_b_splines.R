#### Packages ####
library(splines)
library(MASS)
library(caret)
library(expm)

#### Estimate the gamma in the series estimation ####
kernel_functions <- function(x,ker,knorm="sq")
{
  # Densities such as the integral of the square of the density is one if knorm is sq 
  # or such that the sd is one if knorm is sd.
  if (ker=='normal') 
  { 
    s <- switch(knorm,
                sd = 1,
                sq = 1/(2*sqrt(pi)))
    aux <- exp(- x^2 /(2*s^2)) / (sqrt(2*pi)*s)
  }
  if (ker=='triangle')
  {
    s <- switch(knorm,
                sd = sqrt(6),
                sq = (2/3))
    aux <- (s-abs(x))
    aux <- aux * (aux >0)/s^2
  }
  if (ker=='laplace')
  {
    s <- switch(knorm,
                sd = 1/sqrt(2),
                sq = 1/4)
    aux <- exp(-abs(x)/s)/(2*s)
  }
  return(aux)
}

create_W <- function(W, h=1 ,ker='laplace',knorm="sq",remove=FALSE)
{
  #   ICM and smoothing matrix 
  #   If no bandwidth is provided for the function, h=1 is used  (no smoothing)
  #   The principle diagonal is replaced with zeros if remove = TRUE.
  n_w <- length(W)
  mat <- sapply(W,function(e) W - e*pracma::ones(n_w,1))
  wmat <-  kernel_functions(mat / h,ker=ker,knorm=knorm)/h; # kernel smoothing 
  # principle diagonal of the matrix replaced with zeros
  if (remove==TRUE) wmat <-  wmat - diag(diag(wmat))
  return(wmat/(n_w**2))
}

#ancienne manière de générer la base avec des noeuds sur les quantiles
#create_P_splines <- function(x, Z, J, degree){
#  knots = quantile(Z, probs = seq(0, 1, length.out = J))
#  return(bs(x, knots = knots, degree = degree, intercept = FALSE)[,1:length(knots)+degree-1])}

#nouvelle manière : dyadic B splines
create_dyadic_interval <- function(nb_pts, a, b) {
  # Base case: if nb_pts is 0, return only the endpoints [a, b]
  if (nb_pts == 0) {
    return(c(a, b))
  }
  # Recursively calculate the inner points for nb_pts - 1
  inner_points <- create_dyadic_interval(nb_pts - 1, a, b)
  # Compute the midpoints between consecutive points
  mid_points <- (inner_points[-length(inner_points)] + inner_points[-1]) / 2
  # Combine the points: original inner points and midpoints, sorted
  combined_points <- sort(c(inner_points, mid_points))
  return(combined_points)
}

create_dyadic_P_splines <- function(x, Z, J, degree){
  n_dyadic = log(J - degree +1)/log(2)
  if (as.integer(n_dyadic)!=n_dyadic){
    print("Error, dimension is not valid")
  }
  else{
    a = min(Z)
    b = max(Z)
    n_pts = 2^{n_dyadic} +1 
    knots = create_dyadic_interval(n_dyadic, a, b)[2:2^{n_dyadic}]
    return(bs(x, knots = knots, degree = degree, intercept = FALSE))
  }}


# compute gamma_hat
verify_solve_is_working <- function(mat) {
  result <- tryCatch({
    solve(mat) # Attempt to invert the matrix
    TRUE       # If successful, return TRUE
  }, error = function(e) {
    FALSE      # If an error occurs, return FALSE
  })
  return(result)
}

estimation_gamma <- function(J, W, Z, Y, degree){ # J = (2^{l}-1) + degree
  n = length(Z)
  #compute Omega
  Omega <- create_W(W)
  
  #compute P
  P <- create_dyadic_P_splines(Z, Z, J, degree)
  #compute gamma
  gamma_step = t(P)%*%Omega%*%Y
  mat = t(P)%*%Omega%*%P
  
  bool = verify_solve_is_working(mat)
  
  if (bool == TRUE){
    gamma_step_invert = solve(t(P)%*%Omega%*%P)
    gamma_hat = gamma_step_invert%*%gamma_step}
  else{
    #gamma_step_invert = solve(t(P)%*%Omega%*%P + 0.1*diag(1, J)) #première solution
    gamma_hat = rep(0, J) #2e solution : pour identifier les cas où c'est tout pourri pcq ça marche pas
  }
  
  rm(Omega, gamma_step_invert, gamma_step, mat, P)
  gc()
  
  return(gamma_hat)
}

#### Data simulation #### 
simulate_data_3 <- function(data_param, g_sim, case){
  n = data_param[1]
  rhoev = data_param[2]
  rhowz = data_param[3]
  
  W = rnorm(n, 0, 1)
  V = rnorm(n, 0, 1)
  
  beta = sqrt((rhoev**2)/(1 - rhoev**2))
  Z = (beta*W + V)/sqrt(1+beta**2)
  
  eta = rnorm(n, 0, 1)
  a = sqrt((rhowz**2)/(1 - rhowz**2))
  eps = (a*V + eta)/sqrt(1+beta**2)
  
  Y = g_sim(Z,case) + eps
  
  return(list(W = W, V = V, Z = Z, Y = Y, eps = eps))
}

g_sim_3 <- function(x,case)  
{
  switch(case,  
         x,
         x^2/sqrt(2),
         sqrt(3*sqrt(3))*x*exp(-(x^2)/2),
         sqrt(3)*x*sin(pi*x/2),
         4*exp(-abs(x)) , 
         log(abs(x-1)+1)*ifelse(x>1,1,-1)*sqrt(10/3),
         sqrt(2)*x*cos(pi*x),
         (log(abs(x-1)+1)*ifelse(x>1,1,-1)*sqrt(10/3) - 0.6*x+ (x^3)*2)/8
  )
}


#### Selection of J by cross-validation on M ####
calcul_M_g_hat <- function(gamma_hat, J, Z, Y, W, degree){ #ok testé
  Omega <- create_W(W)
  a = min(Z)
  b = max(Z)
  n = length(Z)
  g_hat_on_Z <- rep(0, n)
  basis <- create_dyadic_P_splines(Z, Z, J, degree)
  g_hat_on_Z <- basis%*%gamma_hat
  
  d <- Y - g_hat_on_Z        # Vector of differences
  M <- sum((d %*% t(d)) * Omega)
  
  rm(Omega, g_hat_on_Z, d, basis)
  gc()
  return(M/(n**2))}

calcul_M_true_g <- function(g, J, Z, Y, W, case){
  M = 0
  n = length(Z)
  Omega <- create_W(W)
  g_hat <- sapply(Z, g, case = case)
  
  d <- Y - g_hat                    
  M <- sum((d %*% t(d)) * Omega) 
  rm(Omega, d, M)
  gc()
  
  return(M/(n**2))}

calcul_M_g_hat_test_sample <- function(g_hat_on_Z_test, Omega, n_test, Y_test){
  M = 0
  for (i in 1:n_test){
    for (j in 1:n_test){
      M = M + (Y_test[i] - g_hat_on_Z_test[i])*(Y_test[j] - g_hat_on_Z_test[j])*Omega[i,j]
  }}
  return(M/(n_test**2))
}

sample_train_test <- function(Z, Y, W, p_train){
  df <- data.frame(Z = Z, Y = Y, W = W)
  
  # Create bins based on the X variable
  bins <- cut(df$Z, breaks = 3, labels = FALSE)
  index_train <- createDataPartition(bins, p = p_train, list = FALSE)
  index_val <- setdiff(seq_len(nrow(df)), index_train)
  
  train_set <- df[index_train, ]
  val_set   <- df[index_val, ]
  
  return(list(
    Z_train = train_set$Z, Y_train = train_set$Y, W_train = train_set$W,
    Z_validation = val_set$Z, Y_validation = val_set$Y, W_validation = val_set$W
  ))
}

optimization_CV_M <- function(Z, W, Y, vect_J_to_test, p_train, degree, x_grid){
  gamma_hat_list <- vector("list", length(vect_J_to_test))
  a = min(Z)
  b = max(Z)
  M_values <- numeric(length(vect_J_to_test))
  n = length(Z)
  
  # split the data in training and test set
  sampled_data <- sample_train_test(Z, Y, W, p_train)
  Z_train <- sampled_data$Z_train
  Y_train <- sampled_data$Y_train
  W_train <- sampled_data$W_train
  Z_validation <- sampled_data$Z_validation
  Y_validation <- sampled_data$Y_validation
  W_validation <- sampled_data$W_validation
  n_train <- length(Z_train)
  
  # compute gamma_hat_train and g_hat_train for the training set
  for (j_index in 1:length(vect_J_to_test)){
    J = vect_J_to_test[j_index]
    #compute gamma_hat_train
    gamma_hat <- estimation_gamma(J,W_train,Z_train,Y_train,degree)
    gamma_hat_list[[j_index]] <- gamma_hat
    
    # compute the prediction of g_hat_train on the test set 
    n_val <- length(Z_validation)
    basis <- create_dyadic_P_splines(Z_validation, Z_train, J, degree)
    g_hat_on_Z_val <-  basis%*%gamma_hat
    
    # compute M and store it 
    Omega_val <- create_W(W_validation)
    M_values[j_index] <- calcul_M_g_hat_test_sample(g_hat_on_Z_val, Omega_val, n_val, Y_validation)
  }
  
  # return J corresponding to the smallest MSE
  J_opt = vect_J_to_test[which.min(M_values)]

  #evaluate hat_g_hat_J on a grid x
  gamma_hat <- estimation_gamma(J_opt,W,Z,Y,degree)
  
  g_hat_on_x_opt = rep(0, length(x_grid))
  basis <- create_dyadic_P_splines(x_grid, Z, J_opt, degree)
  g_hat_on_x_opt <- basis%*%gamma_hat
  
  rm(basis, Z_train, Z_validation, Y_train, Y_validation, W_train, W_validation)
  gc()
  return(list(J_opt = J_opt, gamma_hat_j_opt = gamma_hat, g_hat_on_x_opt = g_hat_on_x_opt, gamma_hat_list = gamma_hat_list, M_values = M_values))
  
}

#### Selection of J by cross-validation on the MSE ####
optimization_CV_MSE <- function(Z, W, Y, vect_J_to_test, p_train, degree, x_grid){
  gamma_hat_list <- vector("list", length(vect_J_to_test))
  a = min(Z)
  b = max(Z)
  MSE_values <- numeric(length(vect_J_to_test))
  n = length(Z)
  
  # split the data in training and test set
  sampled_data <- sample_train_test(Z, Y, W, p_train)
  Z_train <- sampled_data$Z_train
  Y_train <- sampled_data$Y_train
  W_train <- sampled_data$W_train
  Z_validation <- sampled_data$Z_validation
  Y_validation <- sampled_data$Y_validation
  W_validation <- sampled_data$W_validation
  n_train <- length(Z_train)
  
  # compute gamma_hat_train and g_hat_train for the training set
  for (j_index in 1:length(vect_J_to_test)){
    J = vect_J_to_test[j_index]
    #compute gamma_hat_train
    gamma_hat <- estimation_gamma(J,W_train,Z_train,Y_train, degree)
    gamma_hat_list[[j_index]] <- gamma_hat
    
    # compute the prediction of g_hat_train on the test set 
    n_val <- length(Z_validation)
    basis <- create_dyadic_P_splines(Z_validation, Z_train, J, degree)
    g_hat_on_Z_val <-  basis%*%gamma_hat
    
    # compute MSE and store it
    MSE = sum((g_hat_on_Z_val - Z_validation)**2)
    MSE_values[j_index] = MSE
  }
  
  # return J corresponding to the smallest MSE
  J_opt = vect_J_to_test[which.min(MSE_values)]
  gamma_hat <- estimation_gamma(J_opt,W,Z,Y,degree)
  
  basis <- create_dyadic_P_splines(x_grid, Z, J_opt, degree)
  g_hat_on_x_opt <- basis%*%gamma_hat
  
  rm(basis, Z_train, Z_validation, Y_train, Y_validation, W_train, W_validation)
  gc()
  
  return(list(J_opt = J_opt, gamma_hat_j_opt = gamma_hat, g_hat_on_x_opt = g_hat_on_x_opt, gamma_hat_list = gamma_hat_list, MSE_values = MSE_values))
  
}



#### Lepski bootstrap ####
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


compute_M_bootstrap <- function(J, W, Z, Y, degree){
  n = length(Z)
  Omega <- create_W(W)
  P <- create_dyadic_P_splines(Z, Z, J, degree)
  
  #compute gamma
  M_boot_step = t(P)%*%Omega
  mat = t(P)%*%Omega%*%P
  
  bool = verify_solve_is_working(mat)
  
  if (bool == TRUE){
    M_boot_step_invert = solve(t(P)%*%Omega%*%P)
    M_boot = M_boot_step_invert%*%M_boot_step}
  else{
    #gamma_step_invert = solve(t(P)%*%Omega%*%P + 0.1*diag(1, J)) #première solution
    M_boot_step_invert = matrix(0, nrow = J, ncol = J) #2e solution : pour identifier les cas où c'est tout pourri pcq ça marche pas
    M_boot = M_boot_step_invert%*%M_boot_step
  }
  
  rm(Omega, P, M_boot_step, mat, M_boot_step_invert)
  gc()
  return(M_boot)
}

calcul_hat_g_on_Z <- function(j, Z, gamma_hat, a, b, Z_to_evaluate, degree){
  g_hat_on_x <- rep(0, length(Z_to_evaluate))
  basis <- create_dyadic_P_splines(Z_to_evaluate, Z, j, degree)
  if (length(gamma_hat)>1){
    g_hat_on_Z <- basis%*%gamma_hat}
  else {
    g_hat_on_Z <- gamma_hat*sum(basis)
  }
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


T_stat_D <- function(x_grid, Y, Z, W, degree, I_hat, w_vector, list_M_boot, matrix_sigma_all_j_all_x, matrix_all_u, matrix_sigma_all_j_j2_all_x){
  I_hat = sort(I_hat)
  length_I_hat = length(I_hat)
  length_x_grid = length(x_grid)
  
  diff_matrix_allx <- replicate(length_x_grid, matrix(0, length_I_hat, length_I_hat), simplify = FALSE)#matrice des diff de D_j
  matrix_ratio_allx <-  replicate(length_x_grid, matrix(0, length_I_hat, length_I_hat), simplify = FALSE)
  
  for (j_index in 1:length(I_hat)){
    j = I_hat[j_index]
    u_hat_j = matrix_all_u[j_index,] #chaque ligne correspond au vecteur des erreurs de prédiction pour J = valeur de la ligne
    u_star_j = u_hat_j*w_vector #on calcule le vecteur bootstrap
    
    P_j_on_x <- create_dyadic_P_splines(x_grid, Z, j, degree)
    D_j_star = P_j_on_x%*%list_M_boot[[j_index]]%*%u_star_j #vecteur de taille nx 
    
    for (j_prime_index in j_index:length(I_hat)){
      j_prime = I_hat[j_prime_index]
      
      #compute u_star_j_prime
      u_hat_j_prime = matrix_all_u[j_prime_index,]
      u_star_j_prime = u_hat_j_prime*w_vector
      
      #compute D_j_prime sur les x
      P_j_prime_on_x <- create_dyadic_P_splines(x_grid, Z, j_prime, degree)
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
    
  return(max(vec_max))}

compute_J_hat <- function(theta_star, I_hat, x_grid, degree, Z, matrix_gamma, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x){
  j_index = 1
  J = I_hat[j_index]
  ratio = calcul_ratio(J, j_index, I_hat, x_grid, degree, Z, matrix_gamma, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x)
  while ((ratio > 10*theta_star && j_index<= length(I_hat))){
    j_index = j_index + 1 
    J = I_hat[j_index]
    ratio = calcul_ratio(J, j_index, I_hat, x_grid,degree, Z, matrix_gamma, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x)
  }
  
  return(J)
}

calcul_ratio <- function(J, J_index, I_hat, x_grid, degree, Z, matrix_gamma, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x){
  length_I_hat = length(I_hat)
  res_mat <- matrix(0, length_I_hat, length(x_grid))
  
  P_j_on_x <- create_dyadic_P_splines(x_grid, Z, J, degree) #on crée la base pour J donné
  if (J>1){
    g_hat_J_x <- P_j_on_x%*%matrix_gamma[J_index,][1:J]
  }else{
    g_hat_J_x <- matrix_gamma[J_index,][1:J]*P_j_on_x}
  
  sigma_J_vect <- matrix_sigma_all_j_all_x[J_index,]
  
  for (J_prime_index in J_index:length(I_hat)){
    J_prime = I_hat[J_prime_index]
    P_j_prime_on_x <- create_dyadic_P_splines(x_grid, Z, J_prime, degree)
    if (J>1){
      g_hat_J_prime_x <- P_j_prime_on_x%*%matrix_gamma[J_prime_index,][1:J_prime]
    }else{
      g_hat_J_prime_x <- matrix_gamma[J_prime_index,][1:J_prime]*P_j_prime_on_x}
    
    sigma_J_prime_vect <- matrix_sigma_all_j_all_x[J_prime_index,]
    for (x_index in 1:length(x_grid)){
      sigma_J_J_prime <- abs(matrix_sigma_all_j_all_x[J_index,x_index]^{2} + matrix_sigma_all_j_all_x[J_prime_index,x_index]^{2} - 2*matrix_sigma_all_j_j2_all_x[[x_index]][J_index, J_prime_index])
      denom = sqrt(sigma_J_J_prime)
      if (denom!=0){
        res_mat[J_prime_index, x_index] = abs(g_hat_J_prime_x[x_index] - g_hat_J_x[x_index])/(denom)
      }}
  }
  return(max(res_mat))
}

valid_dim_b_splines <- function(J, degree){ #attention modifier pour NS
  if (J-degree+1 <= 1){
    return(FALSE)
  }
  else{
    l = log(J - degree +1)/log(2)
    if (as.integer(l)!=l){
      return(FALSE)
    }else{return(TRUE)}}}

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
    return(J)}}


lepski_bootstrap <- function(n_boot,valid_dim, x_grid, W, Z, Y, degree){#attention : ici x_grid juste pour calculer les sup norm
  n = length(Z)
  J_max = compute_J_max(W, Z, Y, degree, valid_dim)
  I_hat = seq(as.integer(0.1*(log(J_max)^2)), as.integer(J_max), by = 1) 
  I_hat = sort(I_hat[sapply(I_hat,valid_dim, degree)]) #select only the valid dimensions
  length_I_hat = length(I_hat)
  
  # estimate the models (ie the gamma_J and M_boot for all J)
  matrix_gamma = matrix(0, nrow = length_I_hat, ncol = max(I_hat))
  list_M_boot <- vector("list", length_I_hat)
  for (j_index in 1:length_I_hat){
    j = I_hat[j_index]
    M_boot_j <- compute_M_bootstrap(j, W, Z, Y, degree)
    list_M_boot[[j_index]] <- M_boot_j
    gamma_J <- M_boot_j%*%Y
    gamma_J_zeros <- c(gamma_J, rep(0, max(I_hat) - j)) #add zeros for good dimension
    matrix_gamma[j_index,] = gamma_J_zeros #ATTENTION ! prendre que les J premiers termes dans les prochains algos pour pas prendre les zéros!!!
  }
  
  #draw the w_i
  W_boot_matrix <- matrix(rnorm(n_boot * n, mean = 0, sd = 1), nrow = n_boot, ncol = n)
  
  #compute the sigma_J(x) for all values of J and all values of x
  matrix_sigma_all_j_all_x <- matrix(0, nrow = length_I_hat, ncol = length(x_grid))
  matrix_all_u <- matrix(0, nrow = length_I_hat, ncol = n) #n number of data points (ie length(Z))
  for (j_index in 1:length_I_hat){
    j = I_hat[j_index]
    h_hat_j_on_Z = calcul_hat_g_on_Z(j, Z, matrix_gamma[j_index,1:j], a, b, Z, degree)
    u_hat_j = Y - h_hat_j_on_Z 
    matrix_all_u[j_index,] = u_hat_j
    
    P_j_on_x <- create_dyadic_P_splines(x_grid, Z, j, degree) #create for all x_grid
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
    P_j_on_x <- create_dyadic_P_splines(x_grid, Z, j, degree)
    
    for (j_prime_index in (j_index):length_I_hat){
      j_prime = I_hat[j_prime_index]
      
      #calcul de sigma_j,jprime(x) pour tous les x_grid et ensuite on met dans les matrices
      U_j_jp <- create_matrix_U(matrix_all_u[j_index,], matrix_all_u[j_prime_index,])
      P_j_prime_on_x <- create_dyadic_P_splines(x_grid, Z,j_prime,degree)
      
      middle_matrix <- list_M_boot[[j_index]]%*%U_j_jp%*%t(list_M_boot[[j_prime_index]])
      
      if (j == 1 && j_prime == 1){
        vec_sigma_j_jprime <- rep(0, length(x_grid))
        for (x_index in 1:length(x_grid)){
          vec_sigma_j_jprime[x_index] = middle_matrix*P_j_on_x[x_index]*P_j_prime_on_x[x_index]
          matrix_sigma_all_j_j2_all_x[[x_index]][j_index, j_prime_index] = vec_sigma_j_jprime[x_index]
        }
      }
      if (j == 1 && j_prime >1){
        vec_sigma_j_jprime <- rep(0, length(x_grid))
        for (x_index in 1:length(x_grid)){
          vec_sigma_j_jprime[x_index] = P_j_on_x[x_index]*(middle_matrix%*%P_j_prime_on_x[x_index,])
          matrix_sigma_all_j_j2_all_x[[x_index]][j_index, j_prime_index] = vec_sigma_j_jprime[x_index]
        }
      }
      if (j>1 && j_prime == 1){
        vec_sigma_j_jprime <- rep(0, length(x_grid))
        for (x_index in 1:length(x_grid)){
          vec_sigma_j_jprime[x_index] = (P_j_on_x[x_index,]%*%middle_matrix)*P_j_prime_on_x[x_index]
          matrix_sigma_all_j_j2_all_x[[x_index]][j_index, j_prime_index] = vec_sigma_j_jprime[x_index]
        }
      }
      if (j>1 && j_prime >1){
        vec_sigma_j_jprime <- rep(0, length(x_grid))
        for (x_index in 1:length(x_grid)){
          vec_sigma_j_jprime[x_index] = P_j_on_x[x_index,]%*%(middle_matrix)%*%P_j_prime_on_x[x_index,]
          matrix_sigma_all_j_j2_all_x[[x_index]][j_index, j_prime_index] = vec_sigma_j_jprime[x_index]
        }
      }
    }}
  
  
  #compute the T stat for each draw of w
  T_i = rep(0, n_boot)
  for (i in 1:n_boot){
    T_i[i] = T_stat_D(x_grid, Y, Z, W, degree, I_hat, W_boot_matrix[i,], list_M_boot, matrix_sigma_all_j_all_x, matrix_all_u, matrix_sigma_all_j_j2_all_x) 
  }
  alpha = min(0.5, 1/J_max)
  theta = quantile(T_i, probs = 1 - alpha)
  
  J_n_hat = I_hat[length(I_hat)-1] #on prend l'avant dernier pour être le plus petit strict
  
  rm(list_M_boot, W_boot_matrix, matrix_all_u)
  J_hat = compute_J_hat(theta,I_hat, x_grid, degree, Z, matrix_gamma, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x)
  
  rm(matrix_gamma,matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x)
  gc()
  J_tilde = min(J_hat, J_n_hat)
  
  
  return(J_tilde)
}

#x = seq(-2, 2, by = 0.1)
#for (i in 1: 100){
#  simul <- simulate_data_3(c(200, 0.5, 0.9), g_sim_3, 2)
#  test <- lepski_bootstrap(50, valid_dim_b_splines, x, simul$W, simul$Z, simul$Y, 3)
#  print(test)}


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

for (i in 1:100){
  print(i)
  simul <- simulate_data_3(c(200, 0.5, 0.9), g_sim_3, 2)
  test <- lepski_chen(10, simul$W, simul$Z, simul$Y, 3, valid_dim_b_splines)}

