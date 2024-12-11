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

create_dyadic_P_splines_bs <- function(x, Z, J, degree){
  n_dyadic = log(J - degree)/log(2)
  #print(c("create_dyadic_P_splines_bs", n_dyadic))
  if (as.integer(n_dyadic)!=n_dyadic){
    print("Error, dimension is not valid")
  }
  else{
    a = min(Z)
    b = max(Z)
    n_pts = 2^{n_dyadic} +1 
    knots = create_dyadic_interval(n_dyadic, a, b)[2:2^{n_dyadic}]
    basis <- bs(x, knots = knots, degree = degree, 
                intercept = TRUE)
    rm(knots, a, b, n_pts)
    gc()
    return(basis) 
  }}

verify_P_is_working <- function(x, knots){
  result <- tryCatch({
    ns(x, knots = knots, Boundary.knots = c(-2,2), intercept = TRUE) # attempt
    TRUE       # If successful, return TRUE
  }, error = function(e) {
    FALSE      # If an error occurs, return FALSE
  })
  return(result)
}

create_dyadic_P_splines_ns <- function(x, Z, J, degree){
  n_dyadic = log(J-1)/log(2)
  #print(c("create_dyadic_P_splines_ns", n_dyadic))
  if (as.integer(n_dyadic)!=n_dyadic){
    print("Error, dimension is not valid")
  }
  else{
    a = min(Z)
    b = max(Z)
    n_pts = 2^{n_dyadic} - 1 
    knots = create_dyadic_interval(n_dyadic, a, b)[2:2^{n_dyadic}] #taille 2^{l} - 1
    bool = verify_P_is_working(x, knots)
    
    if (bool == TRUE){
      return(ns(x, knots = knots, Boundary.knots = c(-2,2), intercept = TRUE))}
    else{
      return(matrix(0, nrow = length(x), ncol = J))
    }
    
  }}

compute_M_bootstrap <- function(J, W, Z, Y, degree, create_P){
  n = length(Z)
  Omega <- create_W(W)
  P <- create_P(Z, Z, J, degree)
  
  #compute gamma
  M_boot_step = t(P)%*%Omega
  mat = t(P)%*%Omega%*%P
  
  bool = verify_solve_is_working(mat)
  
  if (bool == TRUE){
    M_boot_step_invert = solve(mat)
    M_boot = M_boot_step_invert%*%M_boot_step
    rm(mat, M_boot_step_invert)}
  else{
    #gamma_step_invert = solve(t(P)%*%Omega%*%P + 0.1*diag(1, J)) #première solution
    M_boot_step_invert = matrix(0, nrow = J, ncol = J) #2e solution : pour identifier les cas où c'est tout pourri pcq ça marche pas
    M_boot = M_boot_step_invert%*%M_boot_step
    rm(M_boot_step_invert)
  }
  
  rm(Omega,P)
  gc()
  return(M_boot)
}


valid_dim_b_splines <- function(J, degree){ #attention modifier pour NS
  if (J-degree <= 1){
    return(FALSE)
  }
  else{
    l = log(J - degree)/log(2)
    if (as.integer(l)!=l){
      return(FALSE)
    }else{return(TRUE)}}}

valid_dim_ns <- function(J, degree){
  if (J < degree){
    return(FALSE)
  }
  else{
    l = log(J-1)/log(2)
    if (as.integer(l)!=l){
      return(FALSE)
    }else{return(TRUE)}}}

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

estimation_gamma <- function(J, W, Z, Y, degree, create_P){ # J = (2^{l}-1) + degree
  n = length(Z)
  #compute Omega
  Omega <- create_W(W)
  
  #compute P
  P <- create_P(Z, Z, J, degree)
  
  if (sum(P) == 0){
    gamma_hat = rep(0,J)
  }
  else{
    #compute gamma
    gamma_step = t(P)%*%Omega%*%Y
    mat = t(P)%*%Omega%*%P
    
    bool = verify_solve_is_working(mat)
    
    if (bool == TRUE){
      gamma_step_invert = solve(mat)
      gamma_hat = gamma_step_invert%*%gamma_step
      rm(gamma_step_invert)}
    else{
      #gamma_step_invert = solve(t(P)%*%Omega%*%P + 0.1*diag(1, J)) #première solution
      gamma_hat = rep(0, J) #2e solution : pour identifier les cas où c'est tout pourri pcq ça marche pas
    }
    rm(gamma_step, mat)}
  
  rm(Omega, P)
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


#### Other ####
calcul_M_g_hat_test_sample <- function(g_hat_on_Z_test, Omega, n_test, Y_test){
  M = 0
  for (i in 1:n_test){
    for (j in 1:n_test){
      M = M + (Y_test[i] - g_hat_on_Z_test[i])*(Y_test[j] - g_hat_on_Z_test[j])*Omega[i,j]
    }}
  return(M/(n_test**2))
}

#### common lepski ####
calcul_s_J <- function(J, W, Z, Y, degree, create_P){
  n = length(W)
  
  Omega <- create_W(W)
  Phi <- create_P(Z, Z, J, degree)
  if (sum(Phi) == 0){
    return(999)
  }
  else{
    
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
    }}}


compute_J_max <- function(W, Z, Y, degree, bs_bool){
  if (bs_bool == 1){
    J_init = 5
    J_next = 7
    valid_dim = valid_dim_b_splines
    create_P = create_dyadic_P_splines_bs
  }
  else{
    J_init = 3
    J_next = 5
    valid_dim = valid_dim_ns
    create_P = create_dyadic_P_splines_ns
  }
  
  n = length(Z)
  s_1 = 1/calcul_s_J(J_init, W, Z, Y, degree, create_P) #à checker
  s_2 = 1/calcul_s_J(J_next, W, Z, Y, degree, create_P)
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
      s_hat_J = calcul_s_J(J, W, Z, Y, degree, create_P)
      prev_ratio = new_ratio
      new_ratio = s_hat_J / sqrt(n)
    }
    rm(J_init, J_next, n, s_1, s_2, prev_ratio, new_ratio, s_hat_J)
    gc()
    return(J)}}


