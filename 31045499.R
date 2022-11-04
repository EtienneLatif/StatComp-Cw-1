##################
##################
### Question 1 ###
##################
##################

#######
## a ## 
#######

# Import the data set
eco.df <- read.table("ecoStudy.txt", header = TRUE)

#####
# i #
#####

# Check quantiles of density variable
summary(eco.df$density)

# The values asked for in the question are:
# Minimum      : 0.000
# 1st quartile : 1.413
# Median       : 2.744
# Mean         : 3.103
# 3rd quartile : 4.193
# Maximum      : 20.497

######
# ii #
######

# Find the unique values in the column
unique(eco.df$habitat)

# Find the number of observations of each value
dim(eco.df[eco.df$habitat=="A",]) # 100 observations of A
dim(eco.df[eco.df$habitat=="B",]) # 125 observations of B
dim(eco.df[eco.df$habitat=="C",]) # 75 observations of C

# There are 100 observations of type A, 125 observations of 
# type B and 75 observations of type C

#######
## b ##
#######

# Create the log density column
eco.df$logDensity <- log(1 + eco.df$density)

#####
# i #
#####

# Plot density by habitat type
boxplot(eco.df$density ~ eco.df$habitat, 
        main = "Density versus habitat type",
        xlab = "Habitat type",
        ylab = "Density (count per km^2)")

######
# ii #
######

# Plot log density by habitat type
boxplot(eco.df$logDensity ~ eco.df$habitat, 
        main = "Log density versus habitat type",
        xlab = "Habitat type",
        ylab = "Log density")

#######
# iii #
#######

# Display the two boxplots in a 2x1 lattice
par(mfrow = c(1, 2))
boxplot(eco.df$density ~ eco.df$habitat, 
        main = "Density versus habitat type",
        xlab = "Habitat type",
        ylab = "Density (count per km^2)")
boxplot(eco.df$logDensity ~ eco.df$habitat, 
        main = "Log density versus habitat type",
        xlab = "Habitat type",
        ylab = "Log density")

#######
## c ##
#######

# Function to return F-test statistic
f_stat <- function(y, factor){
  # y       : the continuous variable that we wish to test whether is the same
  #           across groups
  # factor  : the categorical variables which with different groups to test
  #           across
  # returns : the F-test statistic is returned
  
  # Find the levels of the factor
  levels <- unique(factor) 
  
  # Find the overall and factor sample means of y
  overall_mean <- mean(y) # Overall sample mean
  factor_means <- factor_means(y, factor, levels) # Sample means by factor.
  
  # Find the total and factor-level number of observations
  K   <- length(levels) # Number of levels of the factor
  N_i <- lapply(levels, function(x) length(y[factor == x])) # Number of
                                                            # observations per
                                                            # factor
  
  # Calculate the numerator
  numerator   <- f_stat_numerator(overall_mean, factor_means, K, N_i)
  
  # Calculate the denominator
  denominator <- f_stat_denominator(y, factor, levels, factor_means, K, N_i)
  
  return(numerator / denominator)
}

# Function to find means by levels of a factor
factor_means <- function(y, factor, levels){
  # y       : the continuous variable that we wish to find means of
  # factor  : the categorical variable with different groups
  # groups  : the groups of the factor
  # returns : a list of means

  factor_means <- rep(NA, length(levels)) # Create an empty list to store the 
                                          # means
  
  # Loop over the levels of the factor and add the means to the list
  i <- 1
  for (level in levels) {
    factor_means[[i]] <- mean(y[factor == level])
    i <- i + 1
  }
  
  return(factor_means)
}

# Function to calculate the numerator of the F-test statistic
f_stat_numerator <- function(overall_mean, factor_means, K, N_i){
  # overall_mean : the overall sample mean of the continuous variable
  # factor_means : the sample mean of the continuous variable by the factor
  # K            : the number of levels of the factor
  # N_i          : the number of observations by factor
  # returns      : the value of the numerator of the F-test statistic
  
  # Calculate the values of the factor-level means minus the overall mean
  # squared
  mean_diff_squared <- (factor_means - overall_mean) ^ 2
  
  # Multiply the mean difference squared values by the number of observations
  # for each factor and sum these values together.
  weighted_mean_diff <- 0
  for (i in 1:K){
    weighted_mean_diff <- weighted_mean_diff + 
                          (N_i[[i]] * mean_diff_squared[[i]])
  }
  
  # Divide the sum of the weighted mean difference squared values by the number
  # of levels of the factor minus one.
  numerator <- weighted_mean_diff / (K - 1)
  
  return(numerator)
}

# Function to calculate the denominator of the F-test statistic
f_stat_denominator <- function(y, factor, levels, factor_means, K, N_i){
  # y            : the continuous variable that we wish to test whether is the 
  #                same across groups
  # factor       : the categorical variables which with different groups to test
  #                across
  # levels       : the levels of the factor
  # factor_means : the sample mean of the continuous variable by the factor
  # K            : the number of levels of the factor
  # N_i          : the number of observations by factor
  # returns      : the value of the denominator of the F-test statistic
  
  denominator <- 0
  for (i in 1:K){
    for (j in 1:N_i[[i]]){
      denominator <- denominator + ((y[factor == levels[[i]]][j] - 
                                     factor_means[i]) ^ 2)
    }
  }
  
  return(denominator)
}

(((100 * (a_mean - overall_mean)^2) + (125 * (b_mean - overall_mean)^2) +
  (75 * (c_mean - overall_mean)^2)) / 2) / 1796.535

# Calculate the F-test statistic for the data
f_stat(eco.df$density, eco.df$habitat)

# The F-test statistic is 0.035

######
# ii #
######

# Find the p-value of the F test
pf(0.352901, 2, 297)

# The p-value of the F-test is 0.297


##################
##################
### Question 2 ###
##################
##################

#######
## a ##
#######

# Function to return the natural log of the joint density of iid Gammma mixture
# variables
calcLogMix2Gamma <- function(x, alpha, beta, p){
  # x       (numeric vector) : vector of inputs 
  # alpha   (numeric vector) : vector of alpha parameters
  # beta    (numeric vector) : vector of beta parameters
  # p       (numeric)        : probability parameter
  # returns (numeric)        : log joint density of the model
  
  # Get the information about validity of the inputs
  error_check <- verify_gamma_input(x, alpha, beta, p)
  
  # If inputs are invalid, terminate the function and display the reasons
  if(error_check[1]){
    stop(error_check[2])
  }
  
  # Create a vector of outputs from the Gamma mixture model pdf
  f_x <- mix2GammaPdf(x, alpha, beta, p)
  
  # Multiply the outputs together to get the joint density
  joint_density <- 1
  for (f_xi in f_x){
    joint_density <- joint_density * f_xi
  }
  
  # Take the natural logarithm of the joint density
  log_density <- log(joint_density)
  
  # Return the output
  return(log_density)
}

# Function to return the density of a single Gamma mixture variable
mix2GammaPdf <- function(x_i, alpha, beta, p){
  # x       (numeric vector) : input 
  # alpha   (numeric vector) : vector of alpha parameters
  # beta    (numeric vector) : vector of beta parameters
  # p       (numeric)        : probability parameter
  # returns (numeric)        : density of a single Gamma mixture model variable
  
  # Extract the individual parameters from their vectors
  alpha1 <- alpha[1]
  alpha2 <- alpha[2]
  beta1  <-  beta[1]
  beta2  <-  beta[2]

  # Find the results of the two Gamma pdfs
  gamma1 <- gammaPdf(x_i, alpha1, beta1)
  gamma2 <- gammaPdf(x_i, alpha2, beta2)
  
  # Find the output of the Gamma mixture model pdf
  output <- (p * gamma1) + ((1 - p) * gamma2)
  
  return(output)
}

# Function to return the density of a single Gamma variable
gammaPdf <- function(x_i, alpha, beta){
  # x       (numeric vector) : input
  # alpha   (numeric vector) : vector of alpha parameters
  # beta    (numeric vector) : vector of beta parameters
  # returns (numeric)        : density of a single Gamma variable
  
  # Calculate the output of a Gamma pdf
  output <- ((x_i ^ (alpha - 1)) * (exp(-1 * beta * x_i)) * (beta ^ alpha)) / 
    gamma(alpha)
  
  return(output)
}

#######
## b ##
#######

# Function to verify the validity of inputs to the Gamma mixture model function
verify_gamma_input <- function(x, alpha, beta, p){
  # x       (numeric vector) : input 
  # alpha   (numeric vector) : vector of alpha parameters
  # beta    (numeric vector) : vector of beta parameters
  # p       (numeric)        : probability parameter
  # returns (list)           : a list with the first entry as a boolean 
  #                            representing whether the inputs are valid (TRUE
  #                            indicates invalid) and the second entry a string
  #                            displaying information about the error
  
  # Default to TRUE and set to FALSE if no errors are found to avoid
  # repetitively reassigning value to TRUE if multiple errors found
  error_check <- c(TRUE)
  
  # Create a list to store messages about errors
  # These error messages will be concatenated into one string after all errors
  # are checked for and added to the error_check object
  error_check_msg <- list()
  
  # Check for correct data type of x
  if(!is.vector(x) | !is.numeric(x)){
    error_check_msg <- append(error_check_msg, "- x must be a numeric vector")
  } 
  
  # Check that all x values are positive
  else if(any(x <= 0)) {
    error_check_msg <- append(error_check_msg, 
                              "- x must contain only positive values")
  } 
  
  # Check for correct data type of alpha
  if(!is.vector(alpha) | !is.numeric(alpha) | length(alpha) != 2) {
    error_check_msg <- append(error_check_msg, 
                              "- alpha must be a numeric vector of length 2")
  }
  
  # Check that both alpha values are positive
  else if(any(alpha <= 0)) {
    error_check_msg <- append(error_check_msg, 
                              "- alpha must contain only positive values")
  }
  
  # Check for correct data type of beta
  if(!is.vector(beta) | !is.numeric(beta) | length(beta) != 2) {
    error_check_msg <- append(error_check_msg, 
                              "- beta must be a numeric vector of length 2")
  }
  
  # Check that both beta values are positive
  else if(any(beta <= 0)) {
    error_check_msg <- append(error_check_msg, 
                              "- beta must contain only positive values")
  }
  
  # Check for the type of p
  if(!is.numeric(p) | length(p) != 1) {
    error_check_msg <- append(error_check_msg, 
                              "- p must be a single numeric value")
  }
  
  # Check that p takes a value in [0,1]
  else if(p < 0 | p > 1){
    error_check_msg <- append(error_check_msg, 
                              "- p must be between 0 and 1 inclusive")
  }
  
  # If any error messages have been logged, concatenate these to a single string
  # and append to the return error_check list
  if(length(error_check_msg) > 0) {
    error_check_msg <- paste(c("Input is invalid because: ", error_check_msg),
                             collapse = " \n ")
    error_check <- append(error_check, error_check_msg)
  }
  
  # Else set the error_check boolean to FALSE as no errors were logged
  else {
    error_check[1] <- FALSE
  }
  
  return(error_check)
}

#######
## c ##
#######

# Assign the inputs and parameters
x     <- c(0.06, 2.32, 4.81, 0.02, 2.33, 2.18, 0.83, 2.45, 2.10, 3.27)
alpha <- c(0.5, 6)
beta  <- c(0.5, 2.5)
p     <- 0.35

# Perform the required calculation
calcLogMix2Gamma(x, alpha, beta, p)

# The natural logarithm of the joint density for the given values is -13.928


##################
##################
### Question 3 ###
##################
##################

#######
## a ##
#######

smoothCurve <- function(x, alpha, beta, p, q, k){
  # x       (numeric value/vector) : input value/vector
  # alpha   (numeric)              : alpha parameter
  # beta    (numeric)              : beta parameter
  # p       (numeric)              : p parameter
  # q       (numeric)              : q parameter
  # k       (numeric)              : k parameter
  # returns (numeric value/vector) : output value/vector of function g(.) 
  #                                  displayed in Eq. (5)
  
  # Get the information about validity of the inputs
  error_check <- verifySmoothInput(x, alpha, beta, p, q, k)
  
  # If inputs are invalid, terminate the function and display the reasons
  if(error_check[1]){
    stop(error_check[2])
  }
  
  # Calculate the input to the log(.) function in Eq. (5)
  numerator   <- (gamma(x + beta) ^ k) * (p ^ alpha) * (q ^ beta)
  denominator <- (gamma(alpha) ^ k) * (gamma(beta) ^ k)
  log_input   <- numerator / denominator
  
  # Calculate the input to the sin(.) function in Eq. (5)
  sin_input   <- (1 / (2 * alpha)) * log(log_input)
  
  # Calculate the output of the function in Eq. (5)
  output      <- sin(sin_input)
  
  return(output)
}

# Function to verify the validity of inputs to the smooth curve function
verifySmoothInput <- function(x, alpha, beta, p, q, k){
  # Default to TRUE and set to FALSE if no errors are found to avoid
  # repetitively reassigning value to TRUE if multiple errors found
  error_check <- c(TRUE)
  
  # Create a list to store messages about errors
  # These error messages will be concatenated into one string after all errors
  # are checked for and added to the error_check object
  error_check_msg <- list()
  
  # Check for correct data type of x
  if(!is.numeric(x)){
    error_check_msg <- append(error_check_msg, 
                              "- x must be a numeric value or numeric vector")
  }
  
  # Check for correct data type of alpha
  if(length(alpha) > 1 | !is.numeric(alpha)) {
    error_check_msg <- append(error_check_msg, 
                              "- alpha must be a single numeric value")
  }
  
  # Check that alpha is positive
  else if(alpha <= 0) {
    error_check_msg <- append(error_check_msg, 
                              "- alpha can only take positive values")
  }
  
  # Check for correct data type of beta
  if(length(beta) > 1 | !is.numeric(beta)) {
    error_check_msg <- append(error_check_msg, 
                              "- beta must be a single numeric value")
  }
  
  # Check that both beta values are positive
  else if(beta <= 0) {
    error_check_msg <- append(error_check_msg, 
                              "- beta can only take positive values")
  }
  
  if(length(p) > 1 | !is.numeric(p)) {
    error_check_msg <- append(error_check_msg, 
                              "- p must be a single numeric value")
  }
  
    if(length(q) > 1 | !is.numeric(q)) {
    error_check_msg <- append(error_check_msg, 
                              "- q must be a single numeric value")
  }

  if(length(k) > 1 | !is.numeric(k)) {
    error_check_msg <- append(error_check_msg, 
                              "- k must be a single numeric value")
  }
  
  # If any error messages have been logged, concatenate these to a single string
  # and append to the return error_check list
  if(length(error_check_msg) > 0) {
    error_check_msg <- paste(c("Input is invalid because: ", error_check_msg),
                             collapse = " \n ")
    error_check <- append(error_check, error_check_msg)
  }
  # Else set the error_check boolean to FALSE as no errors were logged
  else {
    error_check[1] <- FALSE
  }
  
  return(error_check)
}

#######
## b ##
#######

calcMidRiemannLoop <- function(xVec, alpha, beta, p, q, k){
  # xVec    (numeric vector) : input vector
  # alpha   (numeric)        : alpha parameter
  # beta    (numeric)        : beta parameter
  # p       (numeric)        : p parameter
  # q       (numeric)        : q parameter
  # k       (numeric)        : k parameter
  # returns
  
  # Get the number of subintervals
  N <- length(xVec) - 1
  
  # Initialise the sum at zero
  midReimannSum <- 0
  
  # Perform the sum
  for(i in 2:N){
    smoothCurveInput <- (x[i] + x[i-1]) / 2
    midReimannSum    <- midReimannSum + ((x[i] - x[i-1]) * 
                                           smoothCurve(smoothCurveInput))
  }
  
  return(midReimannSum)
}

calcMidRiemannLoop(c(1, 2, 3, 4), 1, 1, 1, 1 ,1)

