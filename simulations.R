#### Import code ####
source("code_b_splines_monte_carlo.R")

#### Simulation : degree = 3, J = 4, n = 200,rho_1 = 0.5, rho_2 = 0.9, case2 ####
# First, evaluate performances for different values of J 
data_param = c(200, 0.5, 0.9)
n_MC = 2000
J = 4
degree = 3
x_evaluation = seq(-2, 2, length.out = 100)
MC_J4_n200_deg3_rhocase1_gcase2 <- MC_fixed_J(J, n_MC, degree, x_evaluation, g_sim_3, 2, data_param)

perf_MC_J4_n200_deg3_rhocase1_gcase2 <- rep(0, 5)
perf_MC_J4_n200_deg3_rhocase1_gcase2[1] = compute_perf(MC_J4_n200_deg3_rhocase1_gcase2, 'M')
perf_MC_J4_n200_deg3_rhocase1_gcase2[2] = compute_perf(MC_J4_n200_deg3_rhocase1_gcase2, 'supnorm')
perf_MC_J4_n200_deg3_rhocase1_gcase2[3] = compute_perf(MC_J4_n200_deg3_rhocase1_gcase2, 'Var')
perf_MC_J4_n200_deg3_rhocase1_gcase2[4] = compute_perf(MC_J4_n200_deg3_rhocase1_gcase2, 'MSE')
perf_MC_J4_n200_deg3_rhocase1_gcase2[5] = compute_perf(MC_J4_n200_deg3_rhocase1_gcase2, 'bias')

save.image(file = "C:/Users/candi/Desktop/ETUDES/2025 - ENSAE 4A - EPFL3A/pdm/code/github/perf_MC_J4_n200_deg3_rhocase1_gcase2.RData")
