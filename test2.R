#### Functions ####

# data generating function
g <- function(z){
  return (z**2)
}

# Fourier transform for the normal density
#fourier_transform <- function(x){
#  return(sqrt(2)*exp(-0.5*x**2))
#}

#Test avec w laplace density 
fourier_transform <- function(x){
  return (0.5*exp(-abs(x)))
}



# Omega 
calcul_omega <- function(W,n){
  Omega = matrix(0, n, n)
  for (i in 1:n){
    for (j in 1:n){
      Omega[i,j] = fourier_transform(W[i]-W[j])/n**2
    }}
  return(Omega)}


# M (inutile à priori)
calcul_M <- function(g_test, gamma_hat, J, n, Z, Y, W, Omega, p_j, a, b, bool_splines, degree, g_hat_on_Z_val){
  M = 0
  if (bool_splines == 0){
    for (i in 1:n){
      for (j in 1:n){
        M = M + (Y[i] - g_test(Z[i], J, p_j, gamma_hat, Z, a,b))*(Y[j] - g_test(Z[j], J, p_j, gamma_hat, Z, a, b))*Omega[i,j]}
    }}
  else{
    for (i in 1:n){
      for (j in 1:n){
        M = M + (Y[i] - g_hat_on_Z_val[i])*(Y[j] - g_hat_on_Z_val[j])*Omega[i,j]
      }
    }
  }
  
  return(M)
}

create_P <- function(J, n, Z, p_j){
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
create_P_splines <- function(x, Z, J, degree){
  knots = quantile(Z, probs = seq(0, 1, length.out = J))
  return(bs(x, knots = knots, degree = degree, intercept = FALSE)[,1:length(knots)+degree-1])}


#create_P_splines <- function(x, Z, J){
#  return(bs(x, degree = J, intercept = FALSE))
#}

#P1 <- create_P_splines(Z,Z,3)

# Splines sans noeuds avec juste les poly
#create_P_splines <- function(x, Z, J){
#  n = length(x)
#  P <- matrix(0, n, J)
#  for (i in 1:n){
#   for (j in 1:J){
#      P[i,j] = (1/factorial(j-1))*x[i]**(j-1)
#    }
#  }
#  return(P)
#}

#create_P_splines <- function(x, Z, J){
#  n = length(Z)
#  P <- matrix(0, n, n+2)
#  for (i in 1:n){
#    P[i,1] = 1
#    P[i,2] = x[i]
#    for (j in 1:n){
#        P[i,j+2] = abs(x[i]-Z[j])**J
#    }
#  }
#  return(P)
#}
#P2 <- create_P_splines(Z,Z,3)

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
estimation_gamma <- function(n, J, W, Z, Y, p_j, bool_splines, degree){
  a = min(Z)
  b = max(Z)
  #compute Omega
  Omega <- calcul_omega(W,n)
  
  #compute P
  if (bool_splines == 1){
    P <- t(create_P_splines(Z, Z, J, degree))
  }
  else{
    P <- create_P(J,n,Z,p_j)
  }
  
  #compute gamma
  gamma_step = P%*%Omega%*%Y
  gamma_step_invert = solve(P%*%Omega%*%t(P))
  gamma_hat = gamma_step_invert%*%gamma_step
  
  return(gamma_hat)
}


#### Plot the function program ####
plot_function_true_est <- function(x, gamma_hat, p_j_hat, true_g, Z, bool_splines, degree){
  a = min(Z)
  b = max(Z)
  g_hat_on_x = rep(0, length(x))
  if (bool_splines == 0){ 
    for (i in 1:length(x)){
      g_hat_on_x[i] <- g_hat(x[i], J, p_j_hat, gamma_hat, Z, a, b)
    }}
  else{
    basis <- create_P_splines(x, Z, J, degree)
    g_hat_on_x <- basis%*%gamma_hat
  }
  g_on_x = g(x)
  plot(x, g_hat_on_x, type = 'l', col = 'black')
  lines(x, g_on_x, type = 'l', col = 'green')
}


#### Generate data ####
n = 1000
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
Omega <- calcul_omega(W,n)


# trigonometric basis
J = 10
gamma_hat <- estimation_gamma(n,J,W,Z,Y,p_j_trigo_0_1, 0, none)

x = seq(-5, 5, by = 0.1)
plot_function_true_est(x, gamma_hat, p_j_trigo_0_1, g, Z, 0, none)


#histogram basis
J = 3
gamma_hat <- estimation_gamma(n,J,W,Z,Y,p_j_hist_0_1,0,none)

x = seq(-5, 5, by = 0.1)
plot_function_true_est(x, gamma_hat, p_j_hist_0_1, g, Z,0,none)


# spline basis 
J = 3
gamma_hat <- estimation_gamma(n,J,W,Z,Y, empty, 1, 3)

x = seq(-5, 5, by = 0.1)
plot_function_true_est(x, gamma_hat, empty , g, Z, 1, 3)



#### Optimization of J with cross validation on M ####
library(caret)
library(splines)

optimization_CV <- function(Z, W, Y, n, vect_J_to_test, p_train, p_j, bool_splines, degree){
  gamma_hat_list <- vector("list", length(vect_J_to_test))
  a = min(Z)
  b = max(Z)
  M_values <- numeric(length(vect_J_to_test))
  
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
    gamma_hat <- estimation_gamma(n_train,vect_J_to_test[j],W_train,Z_train,Y_train,p_j, bool_splines, degree)
    gamma_hat_list[[j]] <- gamma_hat
    
  # compute the prediction of g_hat_train on the test set 
    n_val <- length(Z_validation)
    if (bool_splines == 0){ 
      g_hat_on_Z_val = rep(0, n_val)
      for (i in 1:n_val){
        g_hat_on_Z_val[i] <- g_hat(Z_validation[i], vect_J_to_test[j], p_j, gamma_hat, Z_train, a, b)
      }}
    else{
      basis <- create_P_splines(Z_validation, Z_train, vect_J_to_test[j], degree)
      g_hat_on_Z_val <-  basis%*%gamma_hat
    }
    
    # compute M and store it 
    Omega_val <- calcul_omega(W_validation,n_val)
    M_values[j] <- calcul_M(g_hat, gamma_hat, vect_J_to_test[j], n_val, Z_validation, Y_validation, W_validation, Omega_val, p_j, a, b, bool_splines, degree, g_hat_on_Z_val)
    }
  
  # return J corresponding to the smallest MSE
  J_opt = vect_J_to_test[which.min(M_values)]
  print(c(M_values, J_opt))
  gamma_hat <- estimation_gamma(n,J_opt,W,Z,Y,p_j,bool_splines, degree)
  #gamma_hat <- gamma_hat_list[J_opt] à tester voir si c good pour pas recalcuelr
  
  estimation_points <- rep(0, n)
  if (bool_splines == 0){ 
    for (i in 1:length(Z)){
      estimation_points[i] <- g_hat(Z[i], J_opt, p_j, gamma_hat, Z, a, b)
    }}
  else{
    basis <- create_P_splines(Z, Z, J_opt, degree)
    estimation_points <-  basis%*%gamma_hat
  }
  
  #On plot tous les points, mais on pourrait plot soit sur des points de l'intervalle soit sur seulement les pts de validation
  plot(Z, estimation_points)
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
vect_J_to_test = seq(3, 20, by = 1)
optimization_CV(Z, W, Y, n, vect_J_to_test, 0.8, p_j_trigo_0_1, 0, none)


#### Optimization of J with Adaptive method (Fabienne Comte) ####
# à modifier en fonction de la borne sup que j'aurai trouvée et des propriétés des bases
optimization_adaptive <- function(Z, W, Y, n, vect_J_to_test, p_j, bool_splines, degree, xpenalty){
  a = min(Z)
  b = max(Z)
  
  estimation_points = rep(0,n)
  
  to_min = rep(0, length(vect_J_to_test))
  gamma_squared = rep(0, length(vect_J_to_test))
  penalty = rep(0, length(vect_J_to_test))
  
  for (j in 1:length(vect_J_to_test)){
    j_test <- vect_J_to_test[j]
    gamma_hat <- estimation_gamma(n,j_test,W,Z,Y,p_j, bool_splines, degree)
    gamma_squared[j] = sum(gamma_hat**2)
    penalty[j] = xpenalty * j_test /n
    to_min[j] = penalty[j] - gamma_squared[j]
  }
  
  J_opt = vect_J_to_test[which.min(to_min)]
  print(J_opt)
  print(to_min)
  gamma_hat <- estimation_gamma(n,J_opt,W,Z,Y,p_j, bool_splines, degree)
  if (bool_splines == 0){ 
    for (i in 1:length(Z)){
      estimation_points[i] <- g_hat(Z[i], J_opt, p_j, gamma_hat, Z, a, b)
    }}
  else{
    basis <- create_P_splines(Z, Z, J_opt, degree)
    estimation_points <-  basis%*%gamma_hat
  }
  plot(Z, estimation_points)
  points(Z, Y, col = 'green')
}

vect_J_to_test = seq(2, 20, by = 1)
#dépendant de la pénalisation et de la valeur de ||phi|| on a des formes différentes : algo à adapter + tard
optimization_adaptive(Z, W, Y, n, vect_J_to_test,p_j_trigo_0_1, 0, none, 2000)


