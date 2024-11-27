#### Source code ####
source("code_b_splines_monte_carlo.R")
library(splines)
library(MASS)
library(caret)
library(expm)
library(foreach)
library(doParallel)
library(ggplot2)

#### Import the data and plot, for fixed J ####

#plot_mean_true(res_MC, x_evaluation,J, degree, rhozw, rhouv, case) : plots the average estimations and the true function
#plot_allcurves_true(res_MC, x_evaluation,J, degree, rhozw, rhouv, case) : plots all the estimated curves, the true one and the average one


#### Plot the distribution of the selected values of J_opt ####
#to do 