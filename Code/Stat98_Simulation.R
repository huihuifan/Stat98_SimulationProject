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
# prob = .15, .35, .65
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

run_simulation(num_iters=1000, missing_method="MCAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MCAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MCAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))

run_simulation(num_iters=1000, missing_method="MAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))

run_simulation(num_iters=1000, missing_method="MNAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MNAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MNAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="single", level=.9, true_betas=c(10, .7, .8))

run_simulation(num_iters=1000, missing_method="MCAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MCAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MCAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))

run_simulation(num_iters=1000, missing_method="MAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))

run_simulation(num_iters=1000, missing_method="MNAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MNAR", coeff.miss=coeff.miss, prob[2], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))
run_simulation(num_iters=1000, missing_method="MNAR", coeff.miss=coeff.miss, prob[3], 
               impute_method="multiple", level=.9, true_betas=c(10, .7, .8))


dat <- gen_data()
d.w.m <- genMNAR(df = dat, prop = c(.75, .75, .75), beta.missing = c(1, 1, 1), 
                 vec.col = c(1, 2, 3))
d.f <- d.w.m[complete.cases(d.w.m), ]
#15
fit <- lm(logincome ~ age + edu, data = d.f)
confint(fit, "age")
confint(fit, "edu") # just stops working at a point


# METHOD 3: Multiple Imputation



















par(mfrow = c(1, 1))

head(data)
# MICE's method for extracting the missing data information 
ini<-mice(data,maxit=0)

# Observed values are blue and missing in red
# Display the missing data patern without clustering 
missing.pattern.plot(data,clustered=FALSE) 
# Display the missing data patern with clustering  
missing.pattern.plot(data) 
# Let us order the data
mp.plot(data,y.order=TRUE, x.order=TRUE)

# The influx 
# of a variable quantifies how well its missing data connect to the observed data on 
# other variables. The outflux of a variable quantifies how well its observed data connect 
# to the missing data on other variables. 
fx<-flux(data)
plot(fx$influx,fx$outflux)
textxy(fx$influx,fx$outflux,colnames(data))
x<-seq(0,1,1/101)
y<--x+1
lines(x,y)

fluxplot(data)

# Let's have a look at the flux 
# Outflux is equal to the number of variable pairs with Yj observed and Yk missing, 
# divided by the total number of incomplete data cells. Outflux is an indicator of 
# the potential usefulness of Yj for imputing other variables. 
# Influx is equal to the number of variable pairs (Yj , Yk) with Yj missing and Yk 
# observed, divided by the total number of observed data cells. Influx depends on the 
# proportion of missing data of the variable. Influx of a completely observed variable is 
# equal to 0, whereas for completely missing variables wehave influx = 1


# Run the Imputation


# Visualize Imputation
# Blue is observed
# Red is imputed
stripplot(___, pch = c(1, 20))

# Alternatively we could use a boxplot
bwplot(___)

# Also we can plot the density
densityplot(___)
# Note that we can use a comparison between the observed and the imputed density 
# as a naive check to see if the data is MAR or MNAR (Abayomi, Gelman and Levy 2008 JRSS)
# This is important as MICE ONLY works if the data is MAR. 




d.i <- mice(data.MAR, method = c("pmm", "pmm", "logreg"))
library(mice)








