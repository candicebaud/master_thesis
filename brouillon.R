#### Base de splines #### 
library(splines)
J = 100
x = seq(-1, 5, by = 0.1)
knots = rnorm(5, 0, 1)
b_spline_basis_custom_knots <- bs(x, degree = 5, knots = x, intercept = FALSE)

b_spline_basis_custom_knots <- ns(x, knots = knots, intercept = FALSE)








#### Base de poly #### 
legendre <- function(j,J,z){
  if (j == 0){
    return (1/sqrt(2))
  }
  elif (j<4){
    return (sqrt((2*j+1)/2)/(factorial(j)*2**j)*deriv_poly(j,z))
  }
}

deriv_poly<-function(j,z){
  if (j == 1){
    return (2*z*j*(x**2-1)**(j-1))
  }
  if (j == 2){
    return(((x**2-1)**(j-2))*(x**2*(2*j+4*j*(j-1))-2*j))
  }
  if (j==3){
    return(2*x*(x**2-1)**(j-2)*(x**2*(2*j+4*j*(j-1)*(j+1))-6*j(j-1)))
  }
}


p_j_poly_0_1 <- function(j, J, z){
  
  
}



#### Modification des fonctions pour pas avoir de bool####

#fonctions modifiées
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

g_hat <- function(z, J, p_j, gamma_hat, Z, a, b){
  p_j_vector = rep(0, J)
  for (j in 1:J){
    p_j_vector[j] = basis_a_b(j, J, z, a, b, p_j)
  }
  return(sum(gamma_hat*p_j_vector))
}

calcul_M <- function(g_test, gamma_hat, J, n, Z, Y, W, Omega, p_j, a, b){
  M = 0
  for (i in 1:n){
    for (j in 1:n){
      M = M + (Y[i] - g_test(Z[i], J, p_j, gamma_hat, Z, a,b))*(Y[j] - g_test(Z[j], J, p_j, gamma_hat, Z, a, b))*Omega[i,j]}
  }
  return(M)
}

estimation_gamma <- function(n, J, W, Z, Y, p_j){
  a = min(Z)
  b = max(Z)
  #compute Omega
  Omega <- calcul_omega(W,n)
  
  #compute P
  P <- create_P(J,n,Z,p_j,a,b)
  
  #compute gamma
  gamma_step = P%*%Omega%*%Y
  gamma_step_invert = P%*%Omega%*%t(P)
  gamma_hat = solve(gamma_step_invert)%*%gamma_step
}

plot_function_true_est <- function(x, gamma_hat, p_j_hat, true_g, Z){
  a = min(Z)
  b = max(Z)
  g_hat_on_x = rep(0, length(x))
  for (i in 1:length(x)){
    g_hat_on_x[i] <- g_hat(x[i], J, p_j_hat, gamma_hat, Z, a, b)
  }
  g_on_x = g(x)
  plot(x, g_hat_on_x, type = 'l', col = 'black')
  lines(x, g_on_x, type = 'l', col = 'green')
  
}


# fonctions non modifiées

fourier_transform <- function(x){
  return(sqrt(2)*exp(0.5*x**2))
}


calcul_omega <- function(W,n){
  Omega = matrix(0, n, n)
  for (i in 1:n){
    for (j in 1:n){
      Omega[i,j] = fourier_transform(W[i]-W[j])/n**2
    }}
  return(Omega)}

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

basis_a_b <-function(j, J, z, a, b, p_j_0_1){
  return(p_j_0_1(j,J,((z-a)/(b-a)))/sqrt(b-a))
}

#test 
g <- function(z){
  return (z**3)/sqrt(2)
}
n = 100
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

J = 3
gamma_hat <- estimation_gamma(n,J,W,Z,Y,p_j_trigo_0_1)

x = seq(-5, 5, by = 0.1)
plot_function_true_est(x, gamma_hat, p_j_trigo_0_1, g, Z)

#### OLD ####
#### Functions ####

# data generating function
g <- function(z){
  return (z**3)/sqrt(2)
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
calcul_M <- function(g_test, gamma_hat, J, n, Z, Y, W, Omega, p_j, bool_P){
  M = 0
  for (i in 1:n){
    for (j in 1:n){
      M = M + (Y[i] - g_test(Z[i], J, p_j, gamma_hat, Z, bool_P))*(Y[j] - g_test(Z[j], J, p_j, gamma_hat, Z, bool_P))*Omega[i,j]}
  }
  return(M)
}

create_P <- function(J, n, Z, p_j, bool){
  if (bool == 0){
    P = matrix(0, J, n)
    for (j in 1:J){
      for (i in 1:n){
        P[j,i] = p_j(j, J, Z[i])
      }
    }}
  else{
    a = min(Z)
    b = max(Z)
    P = matrix(0, J, n)
    for (j in 1:J){
      for (i in 1:n){
        P[j,i] = basis_a_b(j, J, Z[i], a, b, p_j)
      }
    }}
  return(P)}

# Compute g_hat
g_hat <- function(z, J, p_j, gamma_hat, Z, bool_P){
  a = min(Z)
  b = max(Z)
  p_j_vector = rep(0, J)
  if (bool_P == 0){
    for (j in 1:J){
      p_j_vector[j] = p_j(j,J,z)
    }
  }
  else{
    for (j in 1:J){
      p_j_vector[j] = basis_a_b(j, J, z, a, b, p_j)
    }
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
estimation_gamma <- function(n, J, W, Z, Y, p_j, bool_P){
  #compute Omega
  Omega <- calcul_omega(W,n)
  
  #compute P
  P <- create_P(J,n,Z,p_j,bool_P)
  
  #compute gamma
  gamma_step = P%*%Omega%*%Y
  gamma_step_invert = P%*%Omega%*%t(P)
  gamma_hat = solve(gamma_step_invert)%*%gamma_step
}

#### Plot the function program ####
plot_function_true_est <- function(x, gamma_hat, p_j_hat, true_g, Z, bool_P){
  g_hat_on_x = rep(0, length(x))
  for (i in 1:length(x)){
    g_hat_on_x[i] <- g_hat(x[i], J, p_j_hat, gamma_hat, Z, bool_P)
  }
  g_on_x = g(x)
  plot(x, g_hat_on_x, type = 'l', col = 'black')
  lines(x, g_on_x, type = 'l', col = 'green')
  
}


#### Generate data ####
n = 100
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
J = 15
gamma_hat <- estimation_gamma(n,J,W,Z,Y,p_j_trigo_0_1, 1)

x = seq(-5, 5, by = 0.1)
plot_function_true_est(x, gamma_hat, p_j_trigo_0_1, g, Z, 1)


#histogram basis
J = 5
gamma_hat <- estimation_gamma(n,J,W,Z,Y,p_j_hist_0_1, 1)

x = seq(-5, 5, by = 0.1)
plot_function_true_est(x, gamma_hat, p_j_hist_0_1, g, Z, 1)











