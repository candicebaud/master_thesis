#### Lepski Chen ####
V_hat_J <- function(J, n, s_J_hat){
  return (sqrt(log(n)/n)/s_J_hat) #the one in our case 
}

lepski_chen <- function(I_hat, c_0, m_m, list_g_hat_J, W, Z, Y, degree, valid_dim, create_P){
  n = length(Z)
  
  length_I_hat = length(I_hat)
  
  J_index = 1
  J = I_hat[J_index]
  condition = FALSE
  while ((J_index<length_I_hat && condition == FALSE)){
    s_J <- calcul_s_J(J, W, Z, Y, degree, create_P)
    V_J = V_hat_J(J, n, s_J)
    J_prime_index = J_index + 1 
    J_prime = I_hat[J_prime_index]
    while ((J_prime_index < length_I_hat && condition == FALSE)){
      s_J_prime <- calcul_s_J(J_prime, W, Z, Y, degree, create_P)
      V_J_prime = V_hat_J(J, n, s_J_prime)
      a = m_m*sum((list_g_hat_J[[J_prime_index]] - list_g_hat_J[[J_index]])^{2})/n
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
  rm(condition, J_index)
  gc()
  J_opt = J
  return(J_opt)
}

