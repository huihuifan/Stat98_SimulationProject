genObs <- function (covariates.matrix, beta, error.term) {
  # gen.obs takes a covariate matrix n by m, 
  # a beta vector (including intercept) 1 by m and a error term
  # and returns a vector 1 by n of observations y = XB + e  
  obs <- covariates.matrix %*% beta + error.term
  return (obs)
}

# gen.MCAR takes a dataframe n by m, a vector of probabilities (default: 0.5),
# a vector of integers that specifies for which columns we would like to generate
# missing data (default: all columns).
genMCAR <- function(df, vec.prob, vec.col) {
  
  if (missing(vec.col) && !missing(vec.prob)) {
    stop("Specify columns in data frame where to induce MCAR")
  }
  
  else if (missing(vec.col) && missing(vec.prob)) {
    vec.col <- c(1:ncol(df))
    vec.prob <- rep(0.5, length(vec.col))
  }
  
  else if (missing(vec.prob) && !missing(vec.col)){
    vec.prob <- rep(0.5, length(vec.col))
  }
  
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

###############################################################################
genMAR <- function(df, prop, vec.col){
  # genMAR takes a dataframe (n by m), a vector (1 by m) of proportions and
  # a vector (1 by m) indicating in which column we should induce MAR.
  # it returns a dataframe containing missing values. For each column selected
  # on average we will have the respective proportion of missing values.
  
  lCol <- length(vec.col)
  betas <- 5 + 5*runif(lCol)
  vec.prob <- matrix(NA, ncol = lCol, nrow = nrow(df))
  vec.beta0 <- rep(NA, lCol)
  
  for (i in 1:lCol){
    
    f <- function (beta0) { 
      prop[i] - mean(1/(1 + exp(-beta0 - 
                                  as.matrix(df[,vec.col[-i]])%*%betas[-i])))
    }
    
    vec.beta0[i] <- uniroot(f, interval = c(-100000,100000))$root
    vec.prob[,i] <- 1/(1 + exp(-vec.beta0[i] - 
                           as.matrix(df[,vec.col[-i]])%*%betas[-i]))   
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
