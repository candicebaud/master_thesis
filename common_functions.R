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
