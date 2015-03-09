###################
# Data Generation #
###################
n = 1000
p = .25

edu <- rbinom(n, 1, p)
epsilon <- rnorm(n, 0, 1)
age <- rep(NA, n)

# Loop through education and simulate from different normals for different
# education levels of the population to create dependency 
for (i in 1:n) {
  if (edu[i] == 0) {
    age[i] <- rnorm(1, 50, 30)
    
  }
  else {
    age[i] <- rnorm(1, 45, 20)    
  }
}

# set the true values of the beta parameters
beta <- c(10, .7, .8)

# generate income data using the covariates
logincome <- beta[1] + beta[2]*age + beta[3]*edu + epsilon

# bind all the covariates and income data into a dataframe
data <- data.frame(age, edu, epsilon, y)

########################
# Generate Missingness #
########################

######################
# Imputation Methods #
######################

# METHOD 1: Complete Case analysis
# Here, we want to remove the rows that have some missing values 
# using dataframe subsetting
complete_data <- data[complete.cases(data), ]


# METHOD 2: Single Imputation

# METHOD 3: Multiple Imputation
library("mice")
library("mi")
library("calibrate")
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


##########################
# Coverage Probabilities #
##########################

in.interval <- function(x, lo, hi) {
  abs(x-(hi+lo)/2) < (hi-lo)/2 
}

calc_covg_prob <- function(n, p, n_iters) {
  results <- rep(NA, n_iters)
  for (i in 1:n_iters) {
    binoms <- rbinom(n,1,p)
    mu <- mean(binoms)
    vari <- var(binoms)
    low <- mu - 1.96*vari*(1/sqrt(n)) 
    high <- mu + 1.96*vari*(1/sqrt(n))
    
    results[i] <- in.interval(p, low, high)
  }
  return(sum(results)/n_iters)
}

for (i in 1:length(n)) {
  for (j in 1:length(p)) {
    print(calc_covg_prob(n[i], p[j], 10000))
  }
}
