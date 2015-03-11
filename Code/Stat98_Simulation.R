####################
# Import libraries #
####################
setwd( "C:/Users/Matteo/Desktop/Stat98_SimulationProject/Code")
library(truncnorm)
library(mice)
library("mi")
library("calibrate")

source("helper_functions.R")


###############################
# Run Simulation Combinations #
###############################

# coefficients used to generate missingness
coeff.miss <- c(3, 3, 3)

# missing probability
prob <- c(.05, .15, .25)

# scenarios:
# MCAR, MAR, MNAR
# complete, single, multiple
# prob = .15, .35, .25
run_simulation(num_iters=1000, missing_method="MCAR", coeff.miss=coeff.miss, prob[1], 
               impute_method="complete", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MCAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="complete", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MCAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="complete", level=.9, true_betas=c(10, .7, .8))


run_simulation(num_iters=1000, missing_method="MAR", coeff.miss=coeff.miss, prob[1], 
               impute_method="complete", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="complete", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="complete", level=.9, true_betas=c(10, .7, .8))

run_simulation(num_iters=1000, missing_method="MNAR", coeff.miss=coeff.miss, prob[1], 
               impute_method="complete", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MNAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="complete", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MNAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="complete", level=.9, true_betas=c(10, .7, .8))

run_simulation(num_iters=100, missing_method="MCAR", coeff.miss=coeff.miss, prob[1], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=100, missing_method="MCAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=100, missing_method="MCAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))

run_simulation(num_iters=100, missing_method="MAR", coeff.miss=coeff.miss, prob[1], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=100, missing_method="MAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=100, missing_method="MAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))

run_simulation(num_iters=100, missing_method="MNAR", coeff.miss=coeff.miss, prob[1], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=100, missing_method="MNAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=100, missing_method="MNAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))

run_simulation(num_iters=100, missing_method="MCAR", coeff.miss=coeff.miss, prob[1], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=100, missing_method="MCAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=100, missing_method="MCAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))

run_simulation(num_iters=100, missing_method="MAR", coeff.miss=coeff.miss, prob[1], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=100, missing_method="MAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=100, missing_method="MAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))

run_simulation(num_iters=100, missing_method="MNAR", coeff.miss=coeff.miss, prob[1], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=100, missing_method="MNAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=100, missing_method="MNAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))


# Run these to view the imputations
# stripplot(imp, pch = c(1, 20))
# bwplot(imp)
  



