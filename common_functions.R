#### R script for Cross validation on M ####
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

compute_new_MC_selection <- function(MC){
  zero_indices <- which(sapply(MC$list_gamma, function(x) all(x == 0)))
  
  if (length(zero_indices)>0){
    filtered_gamma <- MC$list_gamma[-zero_indices]
    filtered_g_hat_on_x <- MC$list_g_hat_on_x[-zero_indices]
    filtered_W <- MC$list_W[-zero_indices]
    filtered_Y <- MC$list_Y[-zero_indices]
    filtered_Z <- MC$list_Z[-zero_indices]
    filtered_J_opt <- MC$list_J_opt[-zero_indices] #attention
    
    new_MC <- list(list_J_opt = filtered_J_opt, list_gamma = filtered_gamma, list_g_hat_on_x = filtered_g_hat_on_x, 
                   list_W = filtered_W, list_Y = filtered_Y, list_Z = filtered_Z, 
                   g_0_on_x = MC$g_0_on_x)
  }else{
    new_MC <- MC
  }
  return(new_MC)
}

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
    sup_norm = mean(sup_norm_vect)
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

calcul_M_g_hat_test_sample <- function(g_hat_on_Z_test, Omega, n_test, Y_test){
  M = 0
  for (i in 1:n_test){
    for (j in 1:n_test){
      M = M + (Y_test[i] - g_hat_on_Z_test[i])*(Y_test[j] - g_hat_on_Z_test[j])*Omega[i,j]
    }}
  return(M/(n_test**2))
}
