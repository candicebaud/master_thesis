#### Kernel matrix ####
#kernel evaluation
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
  return(wmat)
}

#### Basis ####
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

# generalize for interval [a,b]
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

#### Estimation ####
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
g_hat <- function(z, J, p_j, gamma_hat, Z, a, b){
  p_j_vector = rep(0, J)
  for (j in 1:J){
    p_j_vector[j] = basis_a_b(j, J, z, a, b, p_j)
  }
  return(sum(gamma_hat*p_j_vector))
}

#### Plot ####

# plot the function
plot_function_true_est <- function(x, gamma_hat, p_j_hat, true_g, Z, bool_splines, degree, case){
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
  g_on_x = g(x, case)
  plot(x, g_hat_on_x, type = 'l', col = 'black')
  lines(x, g_on_x, type = 'l', col = 'green')
}
  

#### Simulation ####
g <- function(x,case)  
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

Y = g(Z, 2) + eps


#### Tests #### 
J = 20
gamma_hat <- estimation_gamma(J,W,Z,Y,p_j_trigo_0_1, 0, none)

x = seq(-5, 5, by = 0.1)
plot_function_true_est(x, gamma_hat, p_j_trigo_0_1, g, Z, 0, none,2)

  