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

# Habitat can take values "A", "B" and "C"

# Find the number of observations of each value
dim(eco.df[eco.df$habitat=="A",])[1] # 100 observations of A
dim(eco.df[eco.df$habitat=="B",])[1] # 125 observations of B
dim(eco.df[eco.df$habitat=="C",])[1] # 75 observations of C

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
        cex.main = 1.5,
        xlab = "Habitat type",
        ylab = "Density (count per km^2)",
        cex.lab = 1.5,
        col = "#274472",
        pch = 19,
        cex = 1.25,
        outcol = "#41729F")

######
# ii #
######

# Plot log density by habitat type
boxplot(eco.df$logDensity ~ eco.df$habitat, 
        main = "Log density versus habitat type",
        cex.main = 1.5,
        xlab = "Habitat type",
        ylab = "Log density", # log variables are unitless
        cex.lab = 1.5,
        col = "#274472",
        pch = 19,
        cex = 1.25,
        outcol = "#41729F")

#######
# iii #
#######

# Display the two boxplots in a 2x1 lattice
par(mfrow = c(2, 1))
boxplot(eco.df$density ~ eco.df$habitat, 
        main = "Density versus habitat type",
        cex.main = 1.5,
        xlab = "Habitat type",
        ylab = "Density (count per km^2)",
        cex.lab = 1.5,
        col = "#274472",
        pch = 19,
        cex = 1.25,
        outcol = "#41729F")
boxplot(eco.df$logDensity ~ eco.df$habitat, 
        main = "Log density versus habitat type",
        cex.main = 1.5,
        xlab = "Habitat type",
        ylab = "Log density", # log variables are unitless
        cex.lab = 1.5,
        col = "#274472",
        pch = 19,
        cex = 1.25,
        outcol = "#41729F")
dev.off() # reset options for par()

class(eco.df$habitat)

#######
## c ##
#######

# Calculate the F-test statistic for a given continuous variable and factor
fStat <- function(y, factor){
  #
  # y and factor must be in the same data frame
  #
  # y      (numeric vector) : the continuous variable that we wish to test 
  #                           whether is the same across groups
  # factor (vector)         : the categorical variable with different 
  #                           groups to test across
  # returns (numeric)       : the F-test statistic
  
  # Find the number of observations
  overallCount   <- length(y)
  
  # Find the overall mean of the continuous variable
  overallMean    <- mean(y)
  
  # Find the groups of the factor
  groups         <- unique(factor)
  
  # Find the number of groups
  groupCount     <- length(groups)
  
  # Calculate the F-test statistic numerator
  fStatNumerator <- fStatNumerator(y, factor, groups, overallMean, groupCount)
  print(fStatNumerator)
  
  # Calculate the F-test statistic denominator
  fStatDenominator <- fStatDenominator(y, factor, groups, overallCount, groupCount)
  print(fStatDenominator)
  
  # Calculate the F-test statistic
  fStat <- fStatNumerator / fStatDenominator
  
  return(fStat)
}

fStatNumerator <- function(y, factor, groups, overallMean, groupCount){
  # y            (numeric vector)  : the continuous variable that we wish to 
  #                                  test whether is the same across groups
  # factor      (vector)           : the categorical variable with different 
  #                                  groups to test across
  # groups      (character vector) : the group names in the factor we are 
  #                                  testing against
  # overallMean (numeric)          : the overall mean of y
  # groupCount  (numeric)          : the number of groups of the factor
  
  # Initialise the sum at zero
  weightedSumMeanDifference <- 0
  
  # Loop over the groups to perform the summation
  for(group in groups){
    # Find the mean of y for the given group
    groupMean                 <- mean(y[factor == group])
    
    # Find the difference between the group mean and overall mean squared
    meanDifference            <- (groupMean - overallMean) ^ 2
    
    # Weight the squared difference in means by the number of observations of y
    # for the given group
    weightedMeanDifference    <- length(y[factor == group]) * meanDifference
    
    # Add the result to the weighted sum of mean differences
    weightedSumMeanDifference <- weightedSumMeanDifference + 
                                 weightedMeanDifference
  }
  
  # Divide the sum by the degrees of freedom
  fStatNumerator <- weightedSumMeanDifference / (groupCount - 1)
  
  return(fStatNumerator)
}

fStatDenominator <- function(y, factor, groups, overallCount, groupCount){
  
  # Store the total sum of the difference between each observation and the mean
  # of y for the group it belongs to squared
  # Initialise at 0
  sumMeanDifference <- 0
  
  for(group in groups){
    # Get the observations of y for the given group
    ySubset                <- y[factor == group]
    
    # Get the mean of y for the given group
    groupMean              <- mean(ySubset)
    
    # Get the difference between each observation in the subset and the group
    # mean squared
    groupMeanDifference    <- (ySubset - groupMean) ^ 2
    
    # Sum the squared differences
    sumGroupMeanDifference <- sum(groupMeanDifference)
    
    # Add to the total squared differences of observations and their respective
    # group mean
    sumMeanDifference      <- sumMeanDifference + sumGroupMeanDifference
  } 
  
  fStatDenominator <- sumMeanDifference / (overallCount - groupCount)
  
  return(fStatDenominator)
}


# Calculate the F-test statistic for the data
fStat(eco.df$density, eco.df$habitat)

# The F-test statistic is 10.481 (3 d.p.)

######
# ii #
######

# Find the p-value of the F test by 1 - P(F < 10.481)
1 - pf(q = 10.481, df1 = 2, df2 = 297)

# The p-value of the F-test is 3.996x10^-5 (3 d.p.)

#######
## c ##
#######





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
  error_check <- verifyGammaInput(x, alpha, beta, p)
  
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
verifyGammaInput <- function(x, alpha, beta, p){
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
x_0     <- c(0.06, 2.32, 4.81, 0.02, 2.33, 2.18, 0.83, 2.45, 2.10, 3.27)
alpha_0 <- c(0.5, 6)
beta_0  <- c(0.5, 2.5)
p_0     <- 0.35

# Perform the required calculation
calcLogMix2Gamma(x_0, alpha_0, beta_0, p_0)

# The natural logarithm of the joint density for the given values is -13.928
# (3 d.p.)


##################
##################
### Question 3 ###
##################
##################

#######
## a ##
#######

# Function to calculate the output of g(.) as displated in Eq. (5)
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

## NEEDS COMMENTING BELOW ##

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
  
  # Check for correct data type of p
  if(length(p) > 1 | !is.numeric(p)) {
    error_check_msg <- append(error_check_msg, 
                              "- p must be a single numeric value")
  }
  
  # Check for correct data type of q
  if(length(q) > 1 | !is.numeric(q)) {
    error_check_msg <- append(error_check_msg, 
                              "- q must be a single numeric value")
  }
  
  # Check for correct data type of k
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
  # returns (numeric)        : an approximation to the area between a function 
  #                            given by smoothCurve() and the x-axis using the 
  #                            midpoint Riemann sum
  
  # Initialise the sum at zero
  midRiemannSum <- 0
  
  # Perform the sum
  for(i in 2:length(xVec)){
    # Find the midpoint of the subinterval
    midPoint      <- (xVec[i] + xVec[i-1]) / 2

    # Add the width of the subinterval times the absolute value of the function
    # at the midpoint of the interval to the Riemann sum
    # Use absolute value to get positive area for subintervals falling below
    # the x-axis
    midRiemannSum <- midRiemannSum + ((xVec[i] - xVec[i-1]) * 
                                           abs(smoothCurve(midPoint, alpha, 
                                                           beta, p, q, k)))
  }
  
  return(midRiemannSum)
}

#######
## c ##
#######

calcMidRiemann <- function(xVec, alpha, beta, p, q, k){
  # xVec    (numeric vector) : input vector
  # alpha   (numeric)        : alpha parameter
  # beta    (numeric)        : beta parameter
  # p       (numeric)        : p parameter
  # q       (numeric)        : q parameter
  # k       (numeric)        : k parameter
  # returns (numeric)        : an approximation to the area between a function 
  #                            given by smoothCurve() and the x-axis using the 
  #                            midpoint Riemann sum
  
  # Get the information about validity of the inputs
  error_check <- verifyRiemannInput(xVec, alpha, beta, p, q, k)
  
  # If inputs are invalid, terminate the function and display the reasons
  if(error_check[1]){
    stop(error_check[2])
  }
  
  # Get the number of elements in xVec
  N <- length(xVec)
  
  # Create two vectors of left-points and right-points for the subintervals
  leftPoints <- xVec[-N]
  rightPoints <- xVec[-1]
  
  # Use the left and right points to calculate a vector of midpoints
  midPoints <- (leftPoints + rightPoints) / 2
  
  # Calculate the (absolute) height of the curve from the x-axis at each 
  # midpoint
  rectHeights <- abs(smoothCurve(midPoints, alpha, beta, p, q, k))
  
  # Use the left and right points to find the width of each subinterval
  rectWidths <- rightPoints - leftPoints
  
  # Calculate the area of the rectangles formed by each subinterval and the
  # height of the curve at the midpoint of the subinterval
  rectAreas <- rectHeights * rectWidths
  
  # Calculate the total approximate area between the curve and x-axis by summing 
  # the area of each rectangle
  midRiemannSum <- sum(rectAreas)
  
  return(midRiemannSum)
}

#######
## d ##
#######

verifyRiemannInput <- function(xVec, alpha, beta, p, q, k){
  # xVec    (numeric vector) : input vector
  # alpha   (numeric)        : alpha parameter
  # beta    (numeric)        : beta parameter
  # p       (numeric)        : p parameter
  # q       (numeric)        : q parameter
  # k       (numeric)        : k parameter
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
  if(length(xVec) < 2 | !is.numeric(xVec)){
    error_check_msg <- append(error_check_msg, 
                              "- xVec must be a numeric vector with length of 
                              at least 2")
  }
  
  if(!monotoneInceasing(xVec)) {
    error_check_msg <- append(error_check_msg, 
                              "- elements of xVec must be in strictly increasing
                              order")
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
  
  # Check for correct data type of p
  if(length(p) > 1 | !is.numeric(p)) {
    error_check_msg <- append(error_check_msg, 
                              "- p must be a single numeric value")
  }
  
  # Check for correct data type of q
  if(length(q) > 1 | !is.numeric(q)) {
    error_check_msg <- append(error_check_msg, 
                              "- q must be a single numeric value")
  }
  
  # Check for correct data type of k
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

# Check if elements of a numeric vector are monotone increasing (not monotone 
# non-decreasing)
monotoneInceasing <- function(vec){
  # vec     (numeric vector) : vector to check
  # returns (logical)        : boolean value that is TRUE if the elements of the
  #                            numeric vector are monotone increasing and FALSE
  #                            otherwise
  
  # Boolean to store whether elements of the vector are monotone increasing
  monotoneInceasing <- TRUE
  
  # Loop over each pair of elements in the vector
  for(i in 2:length(vec)){
    # Set the boolean to FALSE if the second element is not strictly greater
    # than the first element
    if(vec[i-1] >= vec[i]){
      monotoneInceasing <- FALSE
    }
  }
  
  return(monotoneInceasing)
}

#######
## e ##
#######

# Initialise the inputs
xVec_1  <- seq(from = 2, to = 8.5, by = 0.01)
alpha_1 <- 2.1
beta_1  <- 0.5
p_1     <- 3
q_1     <- 6
k_1     <- 2

# Calculate the approximate integral
calcMidRiemann(xVec_1, alpha_1, beta_1, p_1, q_1, k_1)

# The area is calculated to be 4.737 (3 d.p.)

#######
## f ##
#######

calcMidRiemannAreas <- function(xSeqList, alpha, beta, p, q, k){
  
  totalArea <- 0
  for(xVec in xSeqList){
    totalArea <- totalArea + calcMidRiemann(xVec, alpha, beta, p, q, k)
  }
    
  return(totalArea)
}

#######
## g ##
#######

# Initialise the inputs
xVec_2_1   <- seq(from = 3.5, to = 8, by = 0.01)
xVec_2_2   <- seq(from = 2, to = 4.3, by = 0.1)
xVec_2_3   <- seq(from = 1.07, to = 9.012, by = 0.001)
xSeqList_2 <- list(xVec_2_1, xVec_2_2, xVec_2_3)
alpha_2    <- 1.9
beta_2     <- 0.75
p_2        <- 2.15
q_2        <- 7
k_2        <- 1.65

# Calculate the combined area between the curve and x-axis over the intervals
calcMidRiemannAreas(xSeqList_2, alpha_2, beta_2, p_2, q_2, k_2)

# The total area is 10.847 (3 d.p.)


##################
##################
### Question 4 ###
##################
##################

# Calculates the natural logarithm of the probability of an observed sequence
# of weather states
weatherSeqProb <- function(weatherSeq, trProbs, initProbs){
  # weatherSeq (character vector) : vector representing the sequence of weather
  #                                 states with "s" = sunny, "c" = cloudy, "r" =
  #                                 rainy
  # trProbs    (numeric matrix)   : matrix of the transition probabilities
  # initsProbs (numeric vector)   : vector of initial probabilities
  # returns    (numeric)          : the natural logarithm of the probability of
  #                                 the sequence given the initial and 
  #                                 transition probabilities
  
  
  # Encode the weather sequence as numbers for indexing the vector and array
  weatherSeq <- weatherSeqEncoder(weatherSeq)
  
  # Find the probability of the initial state
  seqProb <- initProbs[weatherSeq[1]]
  
  # Multiply by the transition probabilities
  for(i in 2:length(weatherSeq)){
    seqProb <- seqProb * trProbs[weatherSeq[i-1], weatherSeq[i]]
  }
  
  # Take the natural logarithm of the probability of the sequence
  logSeqProb <- log(seqProb)
  
  return(logSeqProb)
}

# Encodes the weather sequence as numeric values which can be used to index the
# probability data structures
weatherSeqEncoder <- function(weatherSeq){
  # weatherSeq (character vector) : vector representing the sequence of weather
  #                                 states with "s" = sunny, "c" = cloudy, "r" =
  #                                 rainy
  # returns    (numeric vector)   : vector representing the sequence of weather
  #                                 states encoded as 1 = sunny, 2 = cloudy, 3 =
  #                                 rainy

  # Loop over the vector and modify values individually
  for(i in 1:length(weatherSeq)){
    if(weatherSeq[i] == "s"){
      weatherSeq[i] <- 1
    }
    else if(weatherSeq[i] == "c"){
      weatherSeq[i] <- 2
    }
    else if(weatherSeq[i] == "r"){
      weatherSeq[i] <- 3
    }
    else{
      stop("invalid value in weather sequence")
    }
  }
  
  return(as.numeric(weatherSeq))
}

#######
## b ##
#######

# Initialise the inputs
weatherSeq_0 <- c("c", "s", "c", "r", "s", "s")
trProbs_0r   <- c(0.5, 0.4, 0.1, 0.33, 0.35, 0.32, 0.3, 0.3, 0.4) # raw values
trProbs_0    <- matrix(data = trProbs_0r, nrow = 3, byrow = TRUE)
initProbs_0  <- c(0.45, 0.25, 0.3)

# Compute the natural logarithm of the probability of the sequence
weatherSeqProb(weatherSeq_0, trProbs_0, initProbs_0)

# The natural log of the probability of the sequence is -6.448

#######
## c ##
#######

# Calculate the natural logarithm of the probability of the observed sequences
# from a weather and jacket colour HMM
weatherColourProbs <- function(colourSeq, emitProbs, weatherSeq, trProbs, 
                               initProbs){
  # colourSeq  (character vector) : vector representing the sequence of jacket
  #                                 states with "B" = black and "W" = white
  # emitProbs  (numeric matrix)   : matrix of emission probabilities for jacket
  #                                 colour given weather state
  # weatherSeq (character vector) : vector representing the sequence of weather
  #                                 states with "s" = sunny, "c" = cloudy, "r" =
  #                                 rainy
  # trProbs    (numeric matrix)   : matrix of the transition probabilities
  # initsProbs (numeric vector)   : vector of initial probabilities
  # returns    (numeric)          : the natural logarithm of the probability of
  #                                 the sequences given the initial, transition 
  #                                 and emission probabilities
  
  # The probability of the sequence:
  # P(Y_1 | X_1)P(X_1) P(Y_2 | X_2)P(X_2 | X_1) ... P(Y_n | X_n)P(X_n | X_n-1) 
  # can be written as:
  # P(X_1) P(X_2 | X_1) ... P(X_n | X_n-1) P(Y_1 | X_1) ... P(Y_n | X_n)
  # such that we can reuse the previous code to calculate the probability of the
  # "hidden" weather sequence using the transition probabilities and multiply 
  # this result by the emission probabilities of the observed sequence of 
  # jacket colours given weather states to get the joint probability of both
  # sequences
  
  # Calculate the natural logarithm of the probability of the hidden states
  weatherSeqProb <- weatherSeqProb(weatherSeq, trProbs, initProbs)
  
  # Calculate the natural logarithm of the probability of the observed jacket 
  # sequence given the "hidden" weather states
  colourSeqProb <- colourSeqProb(colourSeq, emitProbs, weatherSeq)
  
  # Using log(XY) = log(X) + log(Y) we can obtain the natural logarithm of the 
  # joint probability of the sequences by adding the two natural logs of
  # probabilities already found
  weatherColourProbs <- weatherSeqProb + colourSeqProb
  
  return(weatherColourProbs)
}

# Calculate the natural logarithm of the probability of the observed sequenced
# of jacket colours given the "hidden" sequence of weather states
colourSeqProb <- function(colourSeq, emitProbs, weatherSeq){
  # colourSeq  (character vector) : vector representing the sequence of jacket
  #                                 states with "B" = black and "W" = white
  # emitProbs  (numeric matrix)   : matrix of emission probabilities for jacket
  #                                 colour given weather state
  # weatherSeq (character vector) : vector representing the sequence of weather
  #                                 states with "s" = sunny, "c" = cloudy, "r" =
  #                                 rainy
  # returns    (numeric)          : the probability of the observed sequence of
  #                                 jacket colours given the sequence of
  #                                 "hidden" weather states
  
  
  # Encode the sequences as numbers for indexing the emission probability matrix
  weatherSeq <- weatherSeqEncoder(weatherSeq)
  colourSeq  <- colourSeqEncoder(colourSeq)
  
  # Initialise the probability at 1
  colourSeqProb <- 1
  
  # Multiply by the emission probabilities
  for(i in 1:length(colourSeq)){
    colourSeqProb <- colourSeqProb * emitProbs[weatherSeq[i], colourSeq[i]]
  }
  
  # Take the natural logarithm of the probability of the jacket colour sequence
  logColourSeqProb <- log(colourSeqProb)
  
  return(logColourSeqProb)
}

colourSeqEncoder <- function(colourSeq){
  # colourSeq  (character vector) : vector representing the sequence of jacket
  #                                 states with "B" = black and "W" = white
  # returns    (numeric vector)   : vector representing the sequence of jacket
  #                                 states encoded as 1 = black and 2 = white
  
  # Loop over the vector and modify values individually
  for(i in 1:length(colourSeq)){
    if(colourSeq[i] == "B"){
      colourSeq[i] <- 1
    }
    else if(colourSeq[i] == "W"){
      colourSeq[i] <- 2
    }
    else{
      stop("invalid value in jacket colour sequence")
    }
  }
  
  return(as.numeric(colourSeq))
}

#######
## d ##
#######

# Initialise the inputs
colourSeq_1  <- c("B", "W", "W", "B", "B", "W", "W", "W")
emitProbs_1r <- c(0.2, 0.8, 0.55, 0.45, 0.9, 0.1) # raw values
emitProbs_1  <- matrix(data = emitProbs_1r, nrow = 3, byrow = TRUE)
weatherSeq_1 <- c("r", "s", "c", "r", "c", "r", "s", "s")
trProbs_1r   <- c(0.55, 0.25, 0.2, 0.25, 0.35, 0.4, 0.2, 0.15, 0.65)
trProbs_1    <- matrix(data = trProbs_1r, nrow = 3, byrow = TRUE)
initProbs_1  <- c(0.35, 0.45, 0.2)

# Calculate the natural logarithm of the joint probability of the observed
# sequences
weatherColourProbs(colourSeq_1, emitProbs_1, weatherSeq_1, trProbs_1, 
                   initProbs_1)

# The natural log of the joint probability of the observed sequences is -15.121
