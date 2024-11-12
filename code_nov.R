#### Packages ####
library(splines)
library(MASS)
library(caret)

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
  return(wmat/(n**2))
}

# splines
create_P_splines <- function(x, Z, J, degree){
  knots = quantile(Z, probs = seq(0, 1, length.out = J))
  return(bs(x, knots = knots, degree = degree, intercept = FALSE)[,1:length(knots)+degree-1])}

#trigo
p_j_trigo_0_1 <-function(j, J, z){
  a = 2*pi*j*z
  #print(a)
  if (j==1){
    return(1)
  }
  if (j%%2 == 0){
    cos_val = cos(a)
    return(sqrt(2)*cos_val)
  }
  else{
    sin_val = sin(a)
    return(sqrt(2)*sin_val)
  }
}

#hist
p_j_hist_0_1 <-function(j, J, z){
  l = length(z)
  if (l>1){
    res = rep(0, l)
    for (i in 1:l){
      if (z[i] >= (j - 1) / J && z[i] < j / J) {
        res[i] = sqrt(J)
      } else {
        res[i] = 0
      }
    }
    return(res)}
  else{
    if (z >= (j - 1) / J && z < j / J) {
      return(sqrt(J))
    } else{
      return(0)
    }
  }
}

# generalize for interval [a,b] : en fait ici on peut dire qu'on normalise Z entre 0 et 1
basis_a_b <-function(j, J, z, a, b, p_j_0_1){
  return(p_j_0_1(j,J,((z-a)/(b-a)))/sqrt(b-a))
}

# create P
create_P <- function(z_eval, J, Z, p_j, degree, spline_bool){
  n = length(Z)
  a = min(Z)
  b = max(Z)
  P = matrix(0, n, J)
  if (spline_bool == 0){
    for (j in 1:J){
      for (i in 1:n){
        P[i,j] = basis_a_b(j, J, Z[i], a, b, p_j)
      }
    }
  }
  else{
    P <- create_P_splines(z_eval, Z, J, degree)
  }
  return(P)
}

# compute gamma_hat
estimation_gamma <- function(J, W, Z, Y, p_j, bool_splines, degree){
  n = length(Z)
  a = min(Z)
  b = max(Z)
  #compute Omega
  Omega <- create_W(W)
  
  #compute P
  P <- create_P(Z, J, Z, p_j, degree, bool_splines)
  
  #compute gamma
  gamma_step = t(P)%*%Omega%*%Y
  gamma_step_invert = solve(t(P)%*%Omega%*%P)
  gamma_hat = gamma_step_invert%*%gamma_step
  
  return(gamma_hat)
}

# Compute g_hat
g_hat <- function(z, J, p_j, gamma_hat, a, b){
  p_j_vector = rep(0, J)
  for (j in 1:J){
    p_j_vector[j] = basis_a_b(j, J, z, a, b, p_j)
  }
  return(sum(gamma_hat*p_j_vector))
}


#### Simulation 1 from Chen ####
g_sim_1 <- function(z, case){
  return(sin(4*z)*log(z))
}

simulate_data_1 <- function(data_param, g_sim, case){
  n = data_param[1]
  corr = data_param[2]
  
  U_V = mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, corr, corr, 1),ncol = 2))
  
  U = U_V[,1]
  V = U_V[,2]
  X = rnorm(n, 0, 1)
  W = pnorm(X, 0, 1)
  D = rbinom(n, 1, 0.5)
  Z = pnorm(D*(X+V)+(1-D)*V)

  Y = g_sim(Z, none) + U
  
  return(list(U = U, V = V, X = X, W = W, D = D, Z = Z, Y = Y))
}

#data_param_1 = (n, corr)
n = 1250 #2500, 5000, 10 000
corr = 0.75
data_param_1 = c(n, corr)

simul_1 = simulate_data_1(data_param_1, g_sim_1, none) 


#### Simulation 2 from Chen ####
n = 1250 #2500, 5000, 10 000
Z = runif(n,0,1)
U = rnorm(n,0,1)
g_sim_2 <- function(z, case){
  return(sin(15*pi*z)*cos(z))
}

Y = g_sim_2(Z,none) + U

#attention ici pas d'instrument : à réfléchir




#### Simulation 3 #### 
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

data_param_3 = c(200, 0.5, 0.9)
simul_3 = simulate_data_3(data_param_3, g_sim_3, 2) #utilisent que le cas 2 et 3 dans le papier lapenta


#### Selection of J by cross-validation on M ####
calcul_M_g_hat <- function(gamma_hat, J, Z, Y, W, p_j_hat, bool_splines, degree){
  Omega <- create_W(W)
  M = 0 
  a = min(Z)
  b = max(Z)
  n = length(Z)
  g_hat_on_Z <- rep(0, n)
  if (bool_splines == 0){
    for (i in 1:n){
      g_hat_on_Z = g_hat(Z[i], J, p_j_hat, gamma_hat, a, b)
    }}
  else{
    basis <- create_P_splines(Z, Z, J, degree)
    g_hat_on_Z <- basis%*%gamma_hat
  }
  
  for (i in 1:n){
    for (j in 1:n){
      M = M + (Y[i] - g_hat_on_Z[i])*(Y[j] - g_hat_on_Z[j])*Omega[i,j]
    }
    return(M/(n**2))
  }}


calcul_M_true_g <- function(g, J, Z, Y, W){
  M = 0
  Omega <- create_W(W)
  for (i in 1:n){
    for (j in 1:n){
      M = M + (Y[i] - g(Z[i]))*(Y[j] - g(Z[j]))*Omega[i,j]
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


#we want to minimize M so we choose J such that M is minimal
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

optimization_CV_M <- function(Z, W, Y, vect_J_to_test, p_train, p_j, bool_splines, degree, x_plot, case, true_g){
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
  for (j in 1:length(vect_J_to_test)){
    #compute gamma_hat_train
    gamma_hat <- estimation_gamma(vect_J_to_test[j],W_train,Z_train,Y_train,p_j, bool_splines, degree)
    gamma_hat_list[[j]] <- gamma_hat
    
    # compute the prediction of g_hat_train on the test set 
    n_val <- length(Z_validation)
    if (bool_splines == 0){ 
      g_hat_on_Z_val = rep(0, n_val)
      for (i in 1:n_val){
        g_hat_on_Z_val[i] <- g_hat(Z_validation[i], vect_J_to_test[j], p_j, gamma_hat, a, b)
      }}
    else{
      basis <- create_P_splines(Z_validation, Z_train, vect_J_to_test[j], degree)
      g_hat_on_Z_val <-  basis%*%gamma_hat
    }
    
    # compute M and store it 
    Omega_val <- create_W(W_validation)
    M_values[j] <- calcul_M_g_hat_test_sample(g_hat_on_Z_val, Omega_val, n_val, Y_validation)
  }
  
  # return J corresponding to the smallest MSE
  J_opt = vect_J_to_test[which.min(M_values)]
  print(c(M_values, J_opt))
  gamma_hat <- estimation_gamma(J_opt,W,Z,Y,p_j, bool_splines, degree)
  
  g_hat_on_x_opt = rep(0, length(x_plot))
  if (bool_splines == 0){ 
    for (i in 1:length(x_plot)){
      g_hat_on_x_opt[i] <- g_hat(x_plot[i], J_opt, p_j, gamma_hat, a, b)
    }}
  else{
    basis <- create_P_splines(x_plot, Z, J_opt, degree)
    g_hat_on_x_opt <- basis%*%gamma_hat
  }
  
  g_on_x = true_g(x_plot, case)
  plot(x_plot, g_hat_on_x_opt, type = 'l', col = 'black')
  lines(x_plot, g_on_x, type = 'l', col = 'green')
  
  return(J_opt)
  
}


#simulation 1 
data_param_1 = c(1250, 0.75)
simul_1 <- simulate_data_1(data_param_1, g_sim_1, none)
vect_J_to_test = seq(3, 10, by = 1)
x = seq(min(Z), max(Z), by = 0.1)
optimization_CV_M(simul_1$Z, simul_1$W, simul_1$Y, vect_J_to_test, 0.6, p_j_trigo_0_1, 0, none, x, none, g_sim_1)

vect_J_to_test = seq(3, 10, by = 1)
x = seq(min(Z), max(Z), by = 0.1)
optimization_CV_M(simul_1$Z, simul_1$W, simul_1$Y, vect_J_to_test, 0.6, none, 1, 3, x, none, g_sim_1)


#simulation 3 : avec base trigo, pas l'impression que ça fonctionne de ouf
data_param_3 = c(200, 0.5, 0.9)
simul_3 = simulate_data_3(data_param_3, g_sim_3, 2) 

vect_J_to_test = seq(3, 20, by = 1)
x = seq(min(Z), max(Z), by = 0.1)
optimization_CV_M(simul_3$Z, simul_3$W, simul_3$Y, vect_J_to_test, 0.8, p_j_trigo_0_1, 0, none, x, 2, g_sim_3)


vect_J_to_test = seq(3, 20, by = 1)
x = seq(min(Z), max(Z), by = 0.1)
optimization_CV_M(simul_3$Z, simul_3$W, simul_3$Y, vect_J_to_test, 0.8, none, 1, 3, x, 2, g_sim_3)


#### Selection of J by cross-validation on the MSE ####
optimization_CV_MSE <- function(Z, W, Y, vect_J_to_test, p_train, p_j, bool_splines, degree, x_plot, case, true_g){
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
  for (j in 1:length(vect_J_to_test)){
    #compute gamma_hat_train
    gamma_hat <- estimation_gamma(vect_J_to_test[j],W_train,Z_train,Y_train,p_j, bool_splines, degree)
    gamma_hat_list[[j]] <- gamma_hat
    
    # compute the prediction of g_hat_train on the test set 
    n_val <- length(Z_validation)
    if (bool_splines == 0){ 
      g_hat_on_Z_val = rep(0, n_val)
      for (i in 1:n_val){
        g_hat_on_Z_val[i] <- g_hat(Z_validation[i], vect_J_to_test[j], p_j, gamma_hat, a, b)
      }}
    else{
      basis <- create_P_splines(Z_validation, Z_train, vect_J_to_test[j], degree)
      g_hat_on_Z_val <-  basis%*%gamma_hat
    }
    
    # compute MSE and store it
    MSE = sum((g_hat_on_Z_val - Z_validation)**2)
    MSE_values[j] = MSE
     }
  
  # return J corresponding to the smallest MSE
  J_opt = vect_J_to_test[which.min(MSE_values)]
  print(c(MSE_values, J_opt))
  gamma_hat <- estimation_gamma(J_opt,W,Z,Y,p_j, bool_splines, degree)
  
  g_hat_on_x_opt = rep(0, length(x_plot))
  if (bool_splines == 0){ 
    for (i in 1:length(x_plot)){
      g_hat_on_x_opt[i] <- g_hat(x_plot[i], J_opt, p_j, gamma_hat, a, b)
    }}
  else{
    basis <- create_P_splines(x_plot, Z, J_opt, degree)
    g_hat_on_x_opt <- basis%*%gamma_hat
  }
  
  g_on_x = true_g(x_plot, case)
  plot(x_plot, g_hat_on_x_opt, type = 'l', col = 'black')
  lines(x_plot, g_on_x, type = 'l', col = 'green')
  
  return(J_opt)
  
}

#simulation 1 
data_param_1 = c(1250, 0.75)
simul_1 <- simulate_data_1(data_param_1, g_sim_1, none)
vect_J_to_test = seq(3, 10, by = 1)
x = seq(min(Z), max(Z), by = 0.1)
optimization_CV_MSE(simul_1$Z, simul_1$W, simul_1$Y, vect_J_to_test, 0.6, p_j_trigo_0_1, 0, none, x, none, g_sim_1)

vect_J_to_test = seq(3, 10, by = 1)
x = seq(min(Z), max(Z), by = 0.1)
optimization_CV_MSE(simul_1$Z, simul_1$W, simul_1$Y, vect_J_to_test, 0.6, none, 1, 3, x, none, g_sim_1)


#simulation 3 
data_param_3 = c(200, 0.5, 0.9)
simul_3 = simulate_data_3(data_param_3, g_sim_3, 2) 

vect_J_to_test = seq(3, 20, by = 1)
x = seq(min(Z), max(Z), by = 0.1)
optimization_CV_MSE(simul_3$Z, simul_3$W, simul_3$Y, vect_J_to_test, 0.8, p_j_trigo_0_1, 0, none, x, 2, g_sim_3)


vect_J_to_test = seq(3, 20, by = 1)
x = seq(min(Z), max(Z), by = 0.1)
optimization_CV_MSE(simul_3$Z, simul_3$W, simul_3$Y, vect_J_to_test, 0.8, none, 1, 3, x, 2, g_sim_3)




#### Selection of J by Lepski 1 (Comte) #### 


#### Selection of J by Lepski 2 (Chen) #### 
lepski_chen <- function(c_0,W,Z,Y,p_j,bool_splines,degree, valid_dim){
  prev_ratio = 0
  new_ratio = 1000
  J = 1
  
  #find J_max
  while ((prev_ratio>10 || new_ratio<= 10)){
    s_hat_J = calcul_s_J(todo) #TODO
    prev_ratio = new_ratio
    new_ratio = s_hat_J*J*sqrt(log(J))/sqrt(n) 
    J = J+1
    while (!valid_dim(J)){# as long as the value of J is not valid we look for another one
      J = J+1
     }}

  J_max = J
  
  I_hat = seq(as.integer(0.1 * log(J_max)^2), as.integer(J_max) + 1, by = 1)
  I_hat = I_hat[sapply(I_hat, valid_dim)] #select only the valid dimensions
  
  index = 1
  J = I_hat[index]
  condition = FALSE
  while (J %in% I_hat && condition == FALSE){
    for (J_prime in I_hat){
      if (J_prime>J){
        #compute g_J and g_J'
        gamma_J <- estimation_gamma(J,W,Z,Y,p_j, bool_splines, degree)
        gamma_J_prime <- estimation_gamma(J_prime,W,Z,Y,p_j, bool_splines, degree)
        
        #compute V_hat(J) and V_hat(J')
        V_J = V_hat_J(J)
        V_J_prime = V_hat_J(J_prime)
        
        #compute |g_J - g_J'| (a)
        a = difference_norm(g_J, g_J_prime)
        
        #compute c_0(V_hat(J)+V_hat(J')) (b)
        b = c_0*(V_J + V_J_prime)
        
        if(a<b){
          condition = TRUE}
        else{
          condition = FALSE
        }  
        index = index + 1
        J = I_hat[index]
      }
      else{
        index = index + 1
        J = I_hat[index]
      }}}
  J_opt = J
  return(J_opt)
}



valid_dim_splines <- function(J){#checks if J is a valid dimension in Tau  
  
}

V_hat_J <- function(J, n, s_J_hat){
  return (sqrt(log(n)/n)/s_J_hat) #the one in our case 
}

difference_norm <- function(g_J, g_J_prime){
  # TO DO
}


#### Lepski with bootstrap ####
lepski_bootstrap <- function(n_boot,n,valid_dim,x_grid){
  prev_ratio = 0
  new_ratio = 1000
  J = 1
  
  #find J_max
  while ((prev_ratio>10 || new_ratio<= 10)){
    s_hat_J = calcul_s_J(todo) #TODO
    prev_ratio = new_ratio
    new_ratio = s_hat_J*J*sqrt(log(J))/sqrt(n) 
    J = J+1
    while (!valid_dim(J)){# as long as the value of J is not valid we look for another one
      J = J+1
    }}
  
  J_max = J
  
  I_hat = seq(as.integer(0.1 * log(J_max)^2), as.integer(J_max) + 1, by = 1)
  I_hat = I_hat[sapply(I_hat, valid_dim)] #select only the valid dimensions
  length_I_hat = length(I_hat)
  
  
  # estimate the models (ie the gamma_J and M_boot for all J)
  matrix_gamma = matrix(0, nrow = length_I_hat, ncol = length_I_hat)
  list_M_boot = list()
  for (j in 1:length_I_hat){
    M_boot_j <- compute_M_bootstrap(j, W, Z, Y, p_j, bool_splines, degree)
    list_M_boot[[j]] <- M_boot_j
    
    gamma_J <- M_boot_j%*%Y
    gamma_J_zeros <- c(gamma_J, rep(0, length_I_hat - j))
    matrix_gamma[j,] = gamma_J_zeros
  }
  
  #draw the w_i
  W_boot_matrix <- matrix(rnorm(n_boot * n, mean = 0, sd = 1), nrow = n_boot, ncol = n)
  
  T_i = rep(0, n_boot)
  for (i in 1:n_boot){
    T_i[i] = T_stat(x_grid, I_hat,W_boot_matrix[i,],matrix_gamma)
  }
  alpha = min(0.5, sqrt(log(J_max)/J_max))
  theta = quantile(T_i, probs = 1 - alpha)
  
  J_n_hat = I_hat[length(I_hat)-1] #on prend l'avant dernier pour être le plus petit strict
  
  J_hat = todo
    
  J_tilde = min(J_hat, J_n_hat)
  
  return(J_tilde)
}


compute_M_bootstrap <- function(J, W, Z, Y, p_j, bool_splines, degree){
  n = length(Z)
  a = min(Z)
  b = max(Z)
  #compute Omega
  Omega <- create_W(W)
  
  #compute P
  P <- create_P(Z, J, Z, p_j, degree, bool_splines)
  
  #compute gamma
  M_boot_step = t(P)%*%Omega
  M_boot_step_invert = solve(t(P)%*%Omega%*%P)
  M_boot = M_boot_step_invert%*%M_boot_step
  
  return(M_boot)
}


