#### Functions ####

# data generating function
g <- function(z){
  return (z*exp(-z**2))
}

# Fourier transform for the normal density
fourier_transform <- function(x){
  return(sqrt(2)*exp(0.5*x**2))
}

# Omega 
calcul_omega <- function(W,n){
  Omega = matrix(0, n, n)
  for (i in 1:n){
    for (j in 1:n){
      Omega[i,j] = fourier_transform(W[i]-W[j])/n**2
    }}
  return(Omega)}


# M
calcul_M <- function(g_test, gamma_hat, J, n, Z, Y, W, Omega, p_j, a, b){
  M = 0
  for (i in 1:n){
    for (j in 1:n){
      M = M + (Y[i] - g_test(Z[i], J, p_j, gamma_hat, Z, a,b))*(Y[j] - g_test(Z[j], J, p_j, gamma_hat, Z, a, b))*Omega[i,j]}
  }
  return(M)
}

create_P <- function(J, n, Z, p_j, a, b){
  a = min(Z)
  b = max(Z)
  P = matrix(0, J, n)
  for (j in 1:J){
    for (i in 1:n){
      P[j,i] = basis_a_b(j, J, Z[i], a, b, p_j)
    }
  }
  return(P)}

# create P for the splines
library(splines)
create_P_splines <- function(Z, J){
  return(bs(Z, degree = J, intercept = FALSE))
}

# Compute g_hat
g_hat <- function(z, J, p_j, gamma_hat, Z, a, b){
  p_j_vector = rep(0, J)
  for (j in 1:J){
    p_j_vector[j] = basis_a_b(j, J, z, a, b, p_j)
  }
  return(sum(gamma_hat*p_j_vector))
}



#### Bases ####
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


# generalize for interval [a,b]
basis_a_b <-function(j, J, z, a, b, p_j_0_1){
  return(p_j_0_1(j,J,((z-a)/(b-a)))/sqrt(b-a))
}


#### Gamma program ####
estimation_gamma <- function(n, J, W, Z, Y, p_j, bool_splines){
  a = min(Z)
  b = max(Z)
  #compute Omega
  Omega <- calcul_omega(W,n)
  
  #compute P
  if (bool_splines == 1){
    P <- t(create_P_splines(Z, J))
  }
  else{
    P <- create_P(J,n,Z,p_j,a,b)
  }
  
  #compute gamma
  gamma_step = P%*%Omega%*%Y
  gamma_step_invert = P%*%Omega%*%t(P)
  gamma_hat = solve(gamma_step_invert)%*%gamma_step
}


#### Plot the function program ####
plot_function_true_est <- function(x, gamma_hat, p_j_hat, true_g, Z, bool_splines){
  a = min(Z)
  b = max(Z)
  g_hat_on_x = rep(0, length(x))
  if (bool_splines == 0){ 
    for (i in 1:length(x)){
      g_hat_on_x[i] <- g_hat(x[i], J, p_j_hat, gamma_hat, Z, a, b)
    }}
  else{
    basis <- t(as.matrix(bs(x, degree = J, intercept = FALSE)))
    g_hat_on_x <- t(gamma_hat)%*%basis
  }
  g_on_x = g(x)
  plot(x, g_hat_on_x, type = 'l', col = 'black')
  lines(x, g_on_x, type = 'l', col = 'green')
}


#### Generate data ####
n = 2000
W = rnorm(n, 0, 1)
V = rnorm(n, 0, 1)
rhoev = 0.5
rhowz = 0.9
beta = sqrt((rhoev**2)/(1 - rhoev**2))
Z = (beta*W + V)/sqrt(1+beta**2)

eta = rnorm(n, 0, 1)
a = sqrt((rhowz**2)/(1 - rhowz**2))
eps = (a*V + eta)/sqrt(1+beta**2)

Y = g(Z) + eps



#### Estimation #### 
# trigonometric basis
J = 10
gamma_hat <- estimation_gamma(n,J,W,Z,Y,p_j_trigo_0_1, 0)

x = seq(-5, 5, by = 0.1)
plot_function_true_est(x, gamma_hat, p_j_trigo_0_1, g, Z, 0)


#histogram basis
J = 15
gamma_hat <- estimation_gamma(n,J,W,Z,Y,p_j_hist_0_1,0)

x = seq(-5, 5, by = 0.1)
plot_function_true_est(x, gamma_hat, p_j_hist_0_1, g, Z,0)


# spline basis 
J = 7
gamma_hat <- estimation_gamma(n,J,W,Z,Y, empty, 1)

x = seq(-5, 5, by = 0.1)
plot_function_true_est(x, gamma_hat, empty , g, Z, 1)



#### Optimization of J ####
library(caret)
library(splines)

optimization <- function(Z, W, Y, n, vect_J_to_test, p_train, p_j, bool_splines){
  gamma_hat_list <- vector("list", length(vect_J_to_test))
  a = min(Z)
  b = max(Z)
  MSE_values <- numeric(length(vect_J_to_test))
  
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
    gamma_hat <- estimation_gamma(n_train,j,W_train,Z_train,Y_train,p_j,bool_splines)
    gamma_hat_list[[j]] <- gamma_hat
    
  # compute the prediction of g_hat_train on the test set 
    n_val <- length(Z_validation)
    g_hat_on_Z_val = rep(0, n_val)
    if (bool_splines == 0){ 
      for (i in 1:n_val){
        g_hat_on_Z_val[i] <- g_hat(Z_validation[i], j, p_j, gamma_hat, Z_train, a, b)
      }}
    else{
      basis <- t(as.matrix(bs(Z_validation, degree = j, intercept = FALSE)))
      g_hat_on_Z_val <- t(gamma_hat)%*%basis
    }
    
    # compute M and store it 
    #Omega_val <- calcul_omega(W_validation,n_val)
    #M_values[j] <- calcul_M(g_hat_on_Z_val, gamma_hat, j, n_val, Z_validation, Y_validation, W_validation, Omega_val, p_j, a, b)
    
    MSE_values[j] = sum((g_hat_on_Z_val - Y_validation)**2)
    }
  
  # return J corresponding to the smallest MSE
  J_opt = vect_J_to_test[which.min(MSE_values)]
  print(c(MSE_values, J_opt))
  gamma_hat <- estimation_gamma(n,J_opt,W,Z,Y,p_j,bool_splines)
  
  estimation_points <- rep(0, n)
  if (bool_splines == 0){ 
    for (i in 1:length(x)){
      estimation_points[i] <- g_hat(Z[i], J_opt, p_j, gamma_hat, Z, a, b)
    }}
  else{
    basis <- t(as.matrix(bs(Z, degree = J_opt, intercept = FALSE)))
    estimation_points <- t(gamma_hat)%*%basis
  }
  
  plot(Z, estimation_points, col = 'black')
  points(Z, Y, col = 'green')
}


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


# Test 
vect_J_to_test = seq(2, 15, by = 1)
optimization(Z, W, Y, n, vect_J_to_test, 0.8, p_j_trigo_0_1, 0)


