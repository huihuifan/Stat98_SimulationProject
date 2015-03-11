# gen.MCAR takes a dataframe n by m, a vector of probabilities (default: 0.5),
# a vector of integers that specifies for which columns we would like to generate
# missing data (default: all columns).
genMCAR <- function(df, vec.prob, vec.col) {
  
  NumRowDf <- nrow(df)
  NumColDf <- ncol(df)
  
  missing.matrix <- matrix(rbinom(n = NumRowDf*length(vec.col), 1, vec.prob), 
                           ncol = length(vec.col), byrow = T)

  DfMCAR <- df
  
  for (i in 1:length(vec.col)) {
    if (any(which(missing.matrix[,i] == 1))) {
      DfMCAR[which(missing.matrix[,i] == 1), vec.col[i]] <- NA
    }
  }
  
  return(DfMCAR)
}
d <- gen_data()
m <- genMCAR(d, c(.15,0.15,0.15), c(1,2,3))

###############################################################################

genMAR <- function(df, prop, beta.missing, vec.col){
  # genMAR takes a dataframe (n by m), a vector (1 by m) of proportions and
  # a vector (1 by m) indicating in which columns we should induce MAR.
  # it returns a dataframe containing missing values. For each column selected
  # on average we will have the respective proportion of missing values.
  # Assume first column of dataframe is response
  
  lCol <- length(vec.col)
  
  # We might consider changing this or factoring it out
  vec.prob <- matrix(NA, ncol = lCol, nrow = nrow(df))
  vec.beta0 <- rep(NA, lCol)
  
  for (i in 1:lCol){
    
    f <- function (beta0) { 
      prop[i] - mean(1/(1 + exp(-beta0 - 
                                  as.matrix(df[,vec.col[-i]])%*%beta.missing[-i])))
    }
    
    vec.beta0[i] <- uniroot(f, interval = c(-100000,100000))$root
    
    # I construct the probabilities fitting a logit model
    # taking out column i makes sure that the missingness is at random.
    # The way we computed beta0 makes sure that the expected prop. of 
    # missing values in column i is as desired.
    vec.prob[,i] <- 1/(1 + exp(-vec.beta0[i] - 
                                 as.matrix(df[,vec.col[-i]])%*%beta.missing[-i]))   
  }
  
  missing.matrix <- matrix(rbinom(n = nrow(df)*lCol, 1, c(vec.prob)), 
                           ncol = lCol)
  
  DfMAR <- df
  
  for (i in 1:length(vec.col)) {
    if (any(which(missing.matrix[,i] == 1))) {
      DfMAR[which(missing.matrix[,i] == 1), vec.col[i]] <- NA
    }
  }
  return(DfMAR)
}

#############################################################################
# Exactly the same as genMAR, with the exception that, for a given column
# and a give obs, we construct the probability of that unit to have a missing
# using the complete X matrix (i.e. without excluding that column, otherwise
# missing data is MAR)
##############################################################################

genMNAR <- function(df, prop, beta.missing, vec.col) {
  
  lCol <- length(vec.col)
  vec.prob <- matrix(NA, ncol = lCol, nrow = nrow(df))
  vec.beta0 <- rep(NA, lCol)
  
  for (i in 1:lCol){
    
    f <- function (beta0) { 
      prop[i] - mean(1/(1 + exp(-beta0 - 
                                  as.matrix(df[,vec.col])%*%beta.missing)))
    }
    
    # I hope that the interval makes the function converge
    vec.beta0[i] <- uniroot(f, interval = c(-100000,100000))$root
    vec.prob[,i] <- 1/(1 + exp(-vec.beta0[i] - 
                                 as.matrix(df[,vec.col])%*%beta.missing))   
  }
  
  missing.matrix <- matrix(rbinom(n = nrow(df)*lCol, 1, c(vec.prob)), 
                           ncol = lCol)
  
  missing <- rbinom(length(dat$logincome), 1, .4)
  missing[which(dat$age < 40 | dat$logincome < 50)] <- 0
  
  missing.matrix[, 1] <- missing
  
  DfMNAR <- df
  
  for (i in 1:length(vec.col)) {
    DfMNAR[which(missing.matrix[,i] == 1), vec.col[i]] <- NA
  }
  return(DfMNAR) 
  
}


###################################
# Generate Coverage Probabilities #
###################################

in.interval <- function(x, lo, hi) {
  # Tests if a value is inside an interval
  abs(x-(hi+lo)/2) < (hi-lo)/2 
}

gen_data <- function() {
  # Generates Census data following this model:
  # logincome = beta_0 + beta_1 * age + beta_2 * education + e_i
  
  n <- 100
  p.grad <- .25
  edu <- rbinom(n, 1, p.grad)
  epsilon <- rnorm(n, 0, 5)
  age <- rep(NA, n)
  
  # Draw from truncated normal
  age[which(edu == 0)] <- rtruncnorm(n - sum(edu), a = 10, b = 105, 
                                     mean = 50, sd = 10)
  age[which(edu == 1)] <- rtruncnorm(sum(edu), a = 10, b = 105, 
                                     mean = 45, sd = 10)

  # set the true values of the beta parameters
  beta <- c(10, .7, .8)
  
  # generate income data using the covariates
  logincome <- beta[1] + beta[2]*age + beta[3]*edu + epsilon
  
  # bind all the covariates and income data into a dataframe
  data <- data.frame(logincome, age, edu) 
  colnames(data) <- c("logincome", "age", "edu")
  
  return(data)
}

run_simulation <- function(num_iters=1000, missing_method="MCAR", 
                           coeff.miss="none", prob, impute_method="complete", 
                           level=.9, true_betas=c(10, .7, .8)) {
  
  # Calculates coverage probability for given number of iterations
  # Generates data, missingness, imputations, and computes confidence interval

  # set a different, deterministic seed each time a new dataset is generated
  seed_vec <- 1000:2000
  
  # generate empty vectors to hold the confidence intervals 
  age.ci.res <- rep(NA, num_iters)
  edu.ci.res <- rep(NA, num_iters)
  
  for (i in 1:num_iters) {
    set.seed(seed_vec[i])
    dat <- gen_data()
    # compute confidence intervals
    res <- coverage_probs(iters=1000, data=dat, missing_method=missing_method, coeff.miss=coeff.miss,
                                    prob=c(prob, prob, prob), impute_method=impute_method, level=level,
                                    true_betas=true_betas)
    age.ci.res[i] <- res[1]
    edu.ci.res[i] <- res[2]
  }
  
  # compute coverage probability 
  return(c(sum(age.ci.res)/num_iters, sum(edu.ci.res)/num_iters)) 
  
}

coverage_probs <- function(data, missing_method="MCAR",
                           coeff.miss="none", prob, impute_method="complete", 
                           level=.9, iters=1000, true_betas=c(10, .7, .8)) {
 
  # Returns True/False for if the confidence interval for beta for 
  # age and education contains the true value
  
  col.missing <- c(1:ncol(data))
  
  # given the missingness method, generate the missing data
  if (missing_method == "MCAR") {
    # d.w.m means "data with missingness"
    d.w.m <- genMCAR(df = data, vec.prob = prob, vec.col = col.missing)
  }
  else if (missing_method == "MAR") {
    d.w.m <- genMAR(df = data, prop = prob, beta.missing = coeff.miss, 
                                 vec.col = col.missing)
  }
  else {
    d.w.m <- genMNAR(df = data, prop = prob, beta.missing = coeff.miss, 
                   vec.col = col.missing)
  }
  
  # given the imputation method, generate the dataset for analysis
  if (impute_method == "complete") {
    # d.f means "data filled"
    d.f <- d.w.m[complete.cases(d.w.m), ]
  }
  else if (impute_method == "single") {
    
  }
  else if (impute_method == "multiple") {
    
    dat <- gen_data()
    missing_dat <- genMAR(dat, prop=c(.25, .25, .25), vec.col=c(1, 2, 3), beta.missing=coeff.miss)
    
    missing_dat$edu <- as.factor(missing_dat$edu)
    ini<-mice(missing_dat,maxit=0)
    pred<-quickpred(missing_dat)
    vis<- order(ini$nmis)[1:3]
    mth<-ini$method # mice knows they are factors
    
    imp <- mice(missing_dat, m=5, method = mth, predictorMatrix = pred,
                visitSequence = vis, maxit = 5)  
    
    fit <- with(imp, lm(logincome ~ age + edu)) 
    est <- pool(fit) # pools the values with rubin's combining rules
    
    val <- data.frame(summary(pool(fit)))    
    ci.age <- c(val$est[2] - 1.96*val$se[2], val$est[2] + 1.96*val$se[2])
    ci.edu <- c(val$est[3] - 1.96*val$se[3], val$est[3] + 1.96*val$se[3])
    
    
  }
    
  fit <- lm(logincome ~ age + edu, data=d.f)
  fit.ci.age <- confint(fit, "age", level=level)
  fit.ci.edu <- confint(fit, "edu", level=level)
  
  age_in <- in.interval(true_betas[2], fit.ci.age[1], fit.ci.age[2])
  edu_in <- in.interval(true_betas[3], fit.ci.edu[1], fit.ci.edu[2])
  
  return(c(age_in, edu_in))

}














