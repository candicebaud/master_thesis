#### Base de splines #### 
library(splines)
knots <- c(2, 5, 3)
J = 5
x = seq(-1, 5, by = 0.1)
b_spline_basis_custom_knots <- bs(x, knots = knots, degree = J)

b_spline_basis_custom_knots






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