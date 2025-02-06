# 1. Function: simulate
# Purpose: To generate a dataset with two variables, X and Y. Y is a binary random variable 
#          following a Bernoulli distribution, and X is a continuous random variable whose 
#          distribution depends on the value of Y.
#
# Arguments:
#   n - The number of samples to generate. Default value is 1000.
#   p - The probability of Y being 1. Default value is 0.5, meaning a 50% chance for each outcome of Y.
#
# Returns:
#   A data frame with two columns: 'x' (the continuous variable) and 'y' (the binary variable).
#   The 'x' values are generated conditionally on the value of 'y' according to:
#     - If Y = 0, X is uniformly distributed between -3 and 1.
#     - If Y = 1, X is uniformly distributed between -1 and 3.
simulate = function(n=1000, p = 0.5){
  y = rbinom(n, 1, p)  # Generate binary Y variable from a Bernoulli distribution (p = 0.5 by default)
  x = ifelse(y == 0, runif(n, -3, 1), runif(n, -1, 3))  # Generate X based on the value of Y:
  # If Y = 0, X ~ Unif(-3, 1), If Y = 1, X ~ Unif(-1, 3)
  
  # Return the dataset as a data frame with columns 'x' and 'y'
  data.frame(x = x, y = y)
}


# 2. Function: bayes_classifier
# Purpose: A simple Bayes classifier based on the problem setup where X is a continuous variable 
#          and Y is binary. The classifier assigns Y = 1 if X is between 1 and 3, and Y = 0 otherwise.
#
# Arguments:
#   x - A numeric vector representing the values of the predictor variable X.
#
# Returns:
#   A binary vector where each element corresponds to the predicted class for each element in 'x':
#     - 1 if X is between 1 and 3, 
#     - 0 otherwise.
bayes_classifier = function(x) {
  # Assign class 1 if 1 < x <= 3, else class 0
  ifelse(x > 1 & x <= 3, 1, 0)
}


# 3. Function: reg_fun
# Purpose: The regression function used to compute the probability that Y = 1 given X = x.
#          It computes the posterior probability based on a simple mixture model of two uniform
#          distributions (f1 and f0) corresponding to the two possible values of Y (0 and 1).
#
# Arguments:
#   x - A numeric vector representing the values of the predictor variable X.
#
# Returns:
#   A numeric vector with the computed posterior probabilities P(Y = 1 | X = x) for each value of X.
#   This is the likelihood ratio of two uniform distributions, which is used to determine the class of Y.
reg_fun = function(x) {
  # Calculate the regression function r(x) = P(Y=1 | X=x) using a mixture of two uniform densities
  dunif(x, -1, 3) / (dunif(x, -1, 3) + dunif(x, -3, 1))
}

# 4.


