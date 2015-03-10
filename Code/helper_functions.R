###############################################################################
# TOMORROW I WILL TRY TO OPTIMIZE MY CODE, FOR INSTANCE, I CAN FACTOR OUT
# A LOT OF STUFF, BUT FOR NOW THIS SHOUDL WORK OR AT LEAST I HOPE SO!!
###############################################################################

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
  
  # I should factor out the rest of the code
  missing.matrix <- matrix(rbinom(n = nrow(df)*lCol, 1, c(vec.prob)), 
                           ncol = lCol)
  
  DfMNAR <- df
  
  for (i in 1:length(vec.col)) {
    if (any(which(missing.matrix[,i] == 1))) {
      DfMNAR[which(missing.matrix[,i] == 1), vec.col[i]] <- NA
    }
  }
  return(DfMNAR)
}