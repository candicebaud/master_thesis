#### Packages ####
library(splines)
library(MASS)
library(caret)
library(expm)

#### Estimate the gamma in the series estimation ####
kernel_functions <- function(x,ker='normal',knorm="sq")
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

create_W <- function(W, h=1 ,ker='normal',knorm="sq",remove=FALSE)
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
estimation_gamma <- function(J, W, Z, Y, degree){ # J = (2^{l}-1) + degree
  n = length(Z)
  #compute Omega
  Omega <- create_W(W)
  
  #compute P
  P <- create_dyadic_P_splines(Z, Z, J, degree)
  #compute gamma
  gamma_step = t(P)%*%Omega%*%Y
  gamma_step_invert = solve(t(P)%*%Omega%*%P)
  gamma_hat = gamma_step_invert%*%gamma_step
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


#n = 200 #or 400
#rhoev = 0.5 #0.8, 0.8
#rhowz = 0.9 #0.9, 0.7

data_param_3 = c(1000, 0.5, 0.9)
simul_3 = simulate_data_3(data_param_3, g_sim_3, 2) #utilisent que le cas 2 et 3 dans le papier lapenta
W <- simul_3$W
Y <- simul_3$Y
Z <- simul_3$Z

#### Selection of J by cross-validation on M ####
calcul_M_g_hat <- function(gamma_hat, J, Z, Y, W, degree){
  Omega <- create_W(W)
  M = 0 
  a = min(Z)
  b = max(Z)
  n = length(Z)
  g_hat_on_Z <- rep(0, n)
  basis <- create_dyadic_P_splines(Z, Z, J, degree)
  g_hat_on_Z <- basis%*%gamma_hat
  for (i in 1:n){
    for (j in 1:n){
      M = M + (Y[i] - g_hat_on_Z[i])*(Y[j] - g_hat_on_Z[j])*Omega[i,j]
    }
  }
  return(M/(n**2))}

calcul_M_true_g <- function(g, J, Z, Y, W, case){
  M = 0
  n = length(Z)
  Omega <- create_W(W)
  for (i in 1:n){
    for (j in 1:n){
      M = M + (Y[i] - g(Z[i], case))*(Y[j] - g(Z[j], case))*Omega[i,j]
    }
    return(M/(n**2))
  }}

calcul_M_g_hat_test_sample <- function(g_hat_on_Z_test, Omega, n_test, Y_test){
  M = 0
  for (i in 1:n_test){
    for (j in 1:n_test){
      M = M + (Y_test[i] - g_hat_on_Z_test[i])*(Y_test[j] - g_hat_on_Z_test[j])*Omega[i,j]
    }
    return(M/(n_test**2))
  }}

sample_train_test <- function(Z, Y, W, p_train){
  # Create a single data frame with all variables
  df <- data.frame(Z = Z, Y = Y, W = W)
  
  # Create bins based on the X variable (you can modify the binning logic if needed)
  bins <- cut(df$Z, breaks = 3, labels = FALSE)
  
  # Sample indices for the training set using bins for stratification
  index_train <- createDataPartition(bins, p = p_train, list = FALSE)
  
  # Indices for the validation set (complementary to train indices)
  index_val <- setdiff(seq_len(nrow(df)), index_train)
  
  # Create train and validation sets based on indices
  train_set <- df[index_train, ]
  val_set   <- df[index_val, ]
  
  # Return the training and validation sets
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
    print(J)
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
  #print(c(M_values, J_opt))
  
  #evaluate hat_g_hat_J on a grid x
  gamma_hat <- estimation_gamma(J_opt,W,Z,Y,degree)
  
  g_hat_on_x_opt = rep(0, length(x_grid))
  basis <- create_dyadic_P_splines(x_grid, Z, J_opt, degree)
  g_hat_on_x_opt <- basis%*%gamma_hat

  #g_on_x = true_g(x_plot, case)
  #plot(x_plot, g_hat_on_x_opt, type = 'l', col = 'black')
  #lines(x_plot, g_on_x, type = 'l', col = 'green')
  
  return(list(J_opt = J_opt, gamma_hat_j_opt = gamma_hat, g_hat_on_x_opt = g_hat_on_x_opt, gamma_hat_list = gamma_hat_list, M_values = M_values))
  
}

#vect_J_to_test <- c()
#for (l in 1:3){
#  vect_J_to_test[l] = (2^{l}-1)+3
#}
#optimization_CV_M(Z,W,Y,vect_J_to_test,0.5,3, seq(-2, 2, by = 0.1)) #fonctionne


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
  #print(c(MSE_values, J_opt))
  gamma_hat <- estimation_gamma(J_opt,W,Z,Y,degree)
  
  g_hat_on_x_opt = rep(0, length(x_grid))
  basis <- create_dyadic_P_splines(x_grid, Z, J_opt, degree)
  g_hat_on_x_opt <- basis%*%gamma_hat
  
  #g_on_x = true_g(x_plot, case)
  #plot(x_plot, g_hat_on_x_opt, type = 'l', col = 'black')
  #lines(x_plot, g_on_x, type = 'l', col = 'green')
  
  return(list(J_opt = J_opt, gamma_hat_j_opt = gamma_hat, g_hat_on_x_opt = g_hat_on_x_opt, gamma_hat_list = gamma_hat_list, MSE_values = MSE_values))
  
}

#optimization_CV_MSE(Z,W,Y,vect_J_to_test,0.5,3, seq(-2, 2, by = 0.1)) #fonctionne


#### Lepski (paper bootstrap) ####

calcul_s_J <- function(J, W, Z, Y, degree){
  n = length(W)
  
  Omega <- create_W(W)
  Phi <- create_dyadic_P_splines(Z, Z, J, degree)
  tPhi_Phi <- t(Phi)%*%Phi
  
  inverse_tPhi_Phi <- solve(tPhi_Phi)
  sqrt_inverse_tPhi_Phi <- sqrtm(inverse_tPhi_Phi)
  res_mat <- n*sqrt_inverse_tPhi_Phi%*%t(Phi)%*%Omega%*%Phi%*%sqrt_inverse_tPhi_Phi
  
  min_eigen_val <- min(eigen(res_mat)$values)
  
  return(min_eigen_val)}

compute_M_bootstrap <- function(J, W, Z, Y, degree){
  n = length(Z)
  #compute Omega
  Omega <- create_W(W)
  
  #compute P
  P <- create_dyadic_P_splines(Z, Z, J, degree)
  
  #compute gamma
  M_boot_step = t(P)%*%Omega
  M_boot_step_invert = solve(t(P)%*%Omega%*%P)
  M_boot = M_boot_step_invert%*%M_boot_step
  
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
  
  diff_matrix_allx <- replicate(length(x_grid), matrix(0, length_I_hat, length_I_hat), simplify = FALSE)#matrice des diff de D_j
  matrix_ratio_allx <-  replicate(length(x_grid), matrix(0, length_I_hat, length_I_hat), simplify = FALSE)
  
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
    }
    
    #Stockage des résultats
    for (x_index in 1:length(x_grid)){
      diff_matrix_allx[[x_index]][j_index, j_prime_index] =abs(D_j_prime_star[x_index] - D_j_star[x_index]) 
      
      sigma_j_x = matrix_sigma_all_j_all_x[j_index, x_index]
      sigma_j_prime_x = matrix_sigma_all_j_all_x[j_prime_index, x_index]
      sigma_j_j_prime_x = matrix_sigma_all_j_j2_all_x[[x_index]][j_index, j_prime_index]
      sigma_val_x = sigma_j_x + sigma_j_prime_x - 2*sigma_j_j_prime_x
      
      if (sigma_val_x !=0){ #attention on a des sigma = 0 (bizarre) et donc pb
        matrix_ratio_allx[[x_index]][j_index, j_prime_index] = diff_matrix_allx[[x_index]][j_index,j_prime_index]/sqrt(sigma_val_x)
      }}  
    vec_max <- rep(0, length(x_grid))
    for (x_index in 1:length(x_grid)){
      vec_max[x_index] = max(matrix_ratio_allx[[x_index]])
    }
    
    return(max(vec_max))}}


compute_J_hat <- function(theta_star, I_hat, x_grid, degree, matrix_gamma, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x){
  j_index = 1
  J = I_hat[j_index]
  ratio = calcul_ratio(J, j_index, I_hat, x_grid, degree, matrix_gamma, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x)
  while (ratio > 1.1*theta_star){
    j_index = j_index + 1 
    J = I_hat[j_index]
    ratio = calcul_ratio(J, j_index, I_hat, x_grid,degree, matrix_gamma, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x)
  }
  
  return(J)
}


calcul_ratio <- function(J, J_index, I_hat, x_grid, degree, matrix_gamma, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x){
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

valid_dim <- function(J, degree){
  if (J-degree+1 <= 1){
    return(FALSE)
  }
  else{
    l = log(J - degree +1)/log(2)
    if (as.integer(l)!=l){
      return(FALSE)
    }else{return(TRUE)}}}


lepski_bootstrap <- function(n_boot,valid_dim,x_grid, W, Z, Y, degree){#attention : ici x_grid juste pour calculer les sup norm
  n = length(Z)
  #J_max = compute_J_max(W, Z, Y, degree)
  J_max = 18
  I_hat = seq(as.integer(0.1 * (log(J_max)^2))+1, as.integer(J_max), by = 1) #on commence au moins à 2 parce que 1 c'est pas une bonne solution
  I_hat = sort(I_hat[sapply(I_hat,valid_dim, degree)]) #select only the valid dimensions
  length_I_hat = length(I_hat)
  
  # estimate the models (ie the gamma_J and M_boot for all J)
  matrix_gamma = matrix(0, nrow = length_I_hat, ncol = max(I_hat))
  list_M_boot <- list()
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
  alpha = min(0.5, sqrt(log(J_max)/J_max))
  theta = quantile(T_i, probs = 1 - alpha)
  
  J_n_hat = I_hat[length(I_hat)-1] #on prend l'avant dernier pour être le plus petit strict
  
  J_hat = compute_J_hat(theta,I_hat, x_grid, degree, matrix_gamma, matrix_sigma_all_j_all_x, matrix_sigma_all_j_j2_all_x)
  
  J_tilde = min(J_hat, J_n_hat)
  
  return(J_tilde)
}

x = seq(min(Z), max(Z), by = 0.1)
lepski_bootstrap(100,valid_dim,x, W, Z, Y,3) #ok fonctionne avec splines


#### Lepski simpler version ####
compute_J_max <- function(W, Z, Y, degree){
  J_init = degree
  J_next = degree + 1
  n = length(Z)
  s_1 = 1/calcul_s_J(J_init, W, Z, Y, degree)
  s_2 = 1/calcul_s_J(J_next, W, Z, Y, degree)
  prev_ratio = s_1/sqrt(n) 
  new_ratio = s_2/sqrt(n)
  
  J = J_next
  if ((prev_ratio <= 10 && new_ratio > 10)) {
  } else {
    while (!(prev_ratio <= 10 && new_ratio > 10)) {
      print(c(prev_ratio, new_ratio))
      J = J + 1
      
      # Look for a valid J
      while (!valid_dim(J, degree)) {  
        J = J + 1}
      
      print(J)
      
      # Calculate s_hat_J and update ratios
      s_hat_J = calcul_s_J(J, W, Z, Y, degree)
      prev_ratio = new_ratio
      new_ratio = s_hat_J / sqrt(n)
    }
  }
  return(J)}


lepski_chen <- function(c_0,W,Z,Y,degree, valid_dim){
  n = length(Z)
  x_grid = seq(min(Z), max(Z), length.out = 10*n)
  
  #J_max = compute_J_max(W, Z, Y, degree)
  J_max = 50
  I_hat = seq(as.integer(0.1 * (log(J_max)^2))+1, as.integer(J_max), by = 1) #on commence au moins à 2 parce que 1 c'est pas une bonne solution
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
  J_opt = J
  return(J_opt)
}

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
  
  return(m_m*sum((g_hat_J_prime_x_grid - g_hat_J_x_grid)^{2})/n)
}

lepski_chen(1, W,Z,Y,3, valid_dim )

#### Monte Carlo d'abord sans sélection du paramètre J ####
#on choisit des valeurs de J d'abord et on regarde pour un nombre de samples ce qu'on obtient comme MSE par ex sur un grid


MC_fixed_J <- function(J, n_MC, degree, x_evaluation, g_0, case, data_param){
  list_gamma <- list() #liste des gamma obtenus
  list_g_hat_on_x <- list() #liste de l'estimation de g sur la grille
  list_W <- list()
  list_Y <- list()
  list_Z <- list()
  
  for (n in 1:n_MC){
    #generate data 
    simul <- simulate_data_3(data_param, g_0, case)
    W <- simul$W
    Y <- simul$Y
    Z <- simul$Z
    
    list_W[[n]] <- W
    list_Y[[n]] <- Y
    list_Z[[n]] <- Z
    
    #compute g_hat_J
    gamma_hat_J <- estimation_gamma(J, W, Z, Y, degree)
    list_gamma[[n]] <- gamma_hat_J
    
    #compute the function on the grid
    basis <- create_dyadic_P_splines(x_evaluation, Z, J, degree)
    g_hat_on_x <- basis%*%gamma_hat_J
    list_g_hat_on_x[[n]] <- g_hat_on_x
  }
  
  g_0_on_x <- g_0(x_evaluation, case)
  
  return(list(list_gamma = list_gamma, list_g_hat_on_x = list_g_hat_on_x,
              list_W = list_W, list_Y = list_Y, list_Z = list_Z, g_0_on_x = g_0_on_x))
  }

test <- MC_fixed_J(10, 100, 3, seq(-2, 2, by = 0.1), g_sim_3, 2, c(1000, 0.5, 0.9))

compute_perf <- function(res_MC, measure){
  n_MC = length(res_MC$list_gamma)
  MSE = 0 
  var = 0 
  bias = 0
  sup_norm = 0
  M = 0
  n_val <- length(res_MC$list_g_hat_on_x[[1]])
  if (measure == 'MSE'){
    for (i in 1:n_MC){
      MSE = MSE + sum((res_MC$list_g_hat_on_x[[i]] - res_MC$g_0_on_x)^{2})
    }
    MSE = MSE/(n_MC*n_val)
    return(MSE)
  }
  if (measure == 'Var'){
    avg <- rep(0, n_val)
    for (x in 1:n_val){
      for (n in 1:n_MC){
        avg[x] <- avg[x] + res_MC$list_g_hat_on_x[[n]][x]/n_MC
      }}
    for (n in 1:n_MC){
      var = var + sum((res_MC$list_g_hat_on_x[[n]] - avg)^{2})
    }
    var = var/(n_MC*n_val)
    return(var)
  }
  if (measure == 'bias'){
    n_val <- length(res_MC$list_g_hat_on_x[[1]])
    avg <- rep(0, n_val)
    for (x in 1:n_val){
      for (n in 1:n_MC){
        avg[x] <- avg[x] + res_MC$list_g_hat_on_x[[n]][x]
      }}
    
    for (n in 1:n_MC){
      bias = bias + sum((res_MC$g_0_on_x - avg/n_MC)^{2})
    }
    bias = bias/(n_MC*n_val)
    return(bias)
  }
  if (measure == 'supnorm'){
    sup_norm_vect <- rep(0, n_MC)
    for (n in 1:n_MC){
      sup_norm_vect[n] = max(abs(res_MC$list_g_hat_on_x[[n]] - res_MC$g_0_on_x))
    }
    sup_norm = max(sup_norm_vect)
    return(sup_norm)
  }
  if (measure == 'M'){
    M_vect <- rep(0, n_MC)
    for (n in 1:n_MC){
      Omega <- create_W(res_MC$list_W[[n]])
      M_vect[n] = calcul_M_g_hat_test_sample(res_MC$list_g_hat_on_x[[n]], Omega, n_val, res_MC$list_Y[[n]])
    }
    return (mean(M_vect))
  }
}


compute_perf(test, 'M') #calcul du critère M moyen

compute_perf(test, 'supnorm') #supnorm

compute_perf(test, 'Var')
compute_perf(test, 'MSE')
compute_perf(test, 'bias')

# ok on a bien MSE = var + bias -> à faire : créer grille avec résultats

library(ggplot2)
plot_mean_true <- function(res_MC, x_evaluation,J){
  n_MC = length(res_MC$list_gamma)
  n_val = length(res_MC$list_g_hat_on_x[[1]])
  
  #compute the avg
  avg <- rep(0, n_val)
  for (x in 1:n_val){
    for (n in 1:n_MC){
      avg[x] <- avg[x] + res_MC$list_g_hat_on_x[[n]][x]/n_MC
    }}
  
  data <- data.frame(x_evaluation, avg, g_0_on_x = res_MC$g_0_on_x)
  ggplot(data, aes(x = x_evaluation)) +
    geom_line(aes(y = g_0_on_x, colour = 'True Function')) + # Changed the label
    geom_line(aes(y = avg, colour = 'Estimate by MC')) + # Changed the label
    ylab("Function value") +
    xlab("Evaluation points") + 
    scale_colour_manual(
      name = "Functions", # Change this to your desired legend title
      values = c("True Function" = "blue", "Estimate by MC" = "red") # Customize colors
    ) +
    ggtitle(paste("True Function vs Estimation (n_MC =", n_MC, ", J =", J, ")"))
}

plot_mean_true(test, seq(-2, 2, by = 0.1), 10)


# à faire pour plusieurs paramètres, plusieurs tailles de data set... (modifier les titres etc s'il faut)







