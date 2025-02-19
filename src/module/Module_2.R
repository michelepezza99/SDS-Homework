# 4. Function: generate_data
# Purpose: The data generation function creates a dataset of n samples with k features, where each feature
#          is drawn from a specified distribution (normal, uniform, exponential, or a t-distribution).
#          It allows for flexible simulation of data with different underlying distributional properties,
#          and assigns a class label to each observation.
#
# Arguments:
#   n           - A positive integer specifying the number of samples to generate.
#   mean_vector - A numeric vector of length k specifying the mean for each feature.
#   cov_matrix  - A k x k covariance matrix; its diagonal elements provide the variances for the corresponding features.
#   label       - A value (or factor) that will be assigned as the class label to all generated samples.
#   dist        - A character vector of length k specifying the distribution type for each feature.
#                 Supported options are:
#                   "normal"      : Standard normal distribution.
#                   "uniform"     : Uniform distribution adjusted to match the specified variance.
#                   "exponential" : Exponential distribution (requires a positive mean).
#                   "t<df>"       : t-distribution with <df> degrees of freedom (e.g., "t10").
#
# Returns:
#   A data frame with n rows and (k + 1) columns, where:
#     - The first k columns contain the generated features.
#     - The last column, named 'label', contains the assigned class label for each observation.
#
# Details:
#   - If all entries in 'dist' are "normal", the function leverages 'mvrnorm' for efficient multivariate sampling.
#   - When a mixture of distributions is specified, each feature is generated independently:
#       * "normal": Uses rnorm with the corresponding mean and standard deviation derived from the covariance matrix.
#       * "uniform": Computes the range such that the variance equals (range^2)/12, then uses runif.
#       * "exponential": Uses rexp with a rate equal to 1/mean, ensuring the mean is positive.
#       * "t<df>": Extracts the degrees of freedom from the string, scales the t-distribution to match the variance,
#                 and then shifts it by the mean.
#   - The function stops with an error if an unsupported distribution is specified or if invalid parameters are provided.
generate_data <- function(n, mean_vector, cov_matrix, label, 
                          dist = rep("normal", length(mean_vector))) {
  k <- length(mean_vector)
  if(all(dist == "normal")) {
    data <- mvrnorm(n, mu = mean_vector, Sigma = cov_matrix)
  } else {
    data <- matrix(NA, nrow = n, ncol = k)
    for(i in 1:k) {
      if(dist[i] == "normal") {
        data[, i] <- rnorm(n, mean = mean_vector[i], sd = sqrt(cov_matrix[i, i]))
      } else if(dist[i] == "uniform") {
        sigma <- sqrt(cov_matrix[i, i])
        range <- sqrt(12) * sigma  # so that variance = range^2/12
        a <- mean_vector[i] - range/2
        b <- mean_vector[i] + range/2
        data[, i] <- runif(n, min = a, max = b)
      } else if(dist[i] == "exponential") {
        if(mean_vector[i] <= 0) stop("Exponential distribution requires a positive mean.")
        lambda <- 1 / mean_vector[i]
        data[, i] <- rexp(n, rate = lambda)
      } else if(grepl("^t[0-9]+$", dist[i])) {
        df <- as.numeric(sub("t", "", dist[i]))
        if(df <= 2) stop("Degrees of freedom must be greater than 2 for the variance to exist.")
        scale_factor <- sqrt(cov_matrix[i, i] * (df - 2) / df)
        data[, i] <- rt(n, df = df) * scale_factor + mean_vector[i]
      } else {
        stop(paste("Distribution", dist[i], "not supported."))
      }
    }
  }
  data <- as.data.frame(data)
  data$label <- label
  return(data)
}


# 5. Function: bivariate_density_plot
# Purpose: This function generates an interactive 3D scatter plot that visualizes the bivariate density
#          of data generated from two classes. It simulates data for class 0 and class 1 using the 
#          'generate_data' function, estimates a 2D kernel density over the feature space, and assigns
#          each data point a density value from the nearest grid point. Finally, it creates a 3D scatter
#          plot using plotly, where the x and y axes represent the two features and the z axis represents
#          the estimated density, with different colors indicating the class labels.
#
# Arguments:
#   n0            - A positive integer specifying the number of samples to generate for class 0.
#   n1            - A positive integer specifying the number of samples to generate for class 1.
#   k             - An integer indicating the number of features; typically 2 for bivariate density plotting.
#   mean_0        - A numeric vector of length k providing the means for each feature for class 0.
#   mean_1        - A numeric vector of length k providing the means for each feature for class 1.
#   cov_matrix    - A k x k covariance matrix, where the diagonal elements represent the variance for each feature.
#   dist_choice0  - A character vector of length k specifying the distribution type for each feature of class 0.
#   dist_choice1  - A character vector of length k specifying the distribution type for each feature of class 1.
#   grid_n        - (Optional) A positive integer determining the resolution of the grid for density estimation (default: 100).
#
# Returns:
#   An interactive plotly object containing a 3D scatter plot:
#     - x-axis: Feature X1.
#     - y-axis: Feature X2.
#     - z-axis: Density estimate derived from a 2D kernel density estimation.
#     - Colors: Data points are colored according to their class label (0 or 1).
#
# Details:
#   - The function uses the 'generate_data' function to simulate data for both classes.
#   - After generating the datasets, the data from both classes are combined into a single data frame.
#   - If the first k feature columns are not already named "X1" and "X2", they are renamed accordingly.
#   - The 2D kernel density estimation is performed using MASS::kde2d over a grid defined by grid_n.
#   - For each data point, the function identifies the closest grid point and assigns its density value.
#   - The resulting 3D scatter plot is rendered with plotly, facilitating interactive exploration of the density
#     across the feature space.
bivariate_density_plot <- function(n0, n1, k, 
                                   mean_0, mean_1, 
                                   cov_matrix, 
                                   dist_choice0, dist_choice1,
                                   grid_n = 100) {
  # Generate data for class 0 and class 1 using the assumed generate_data function
  data_0 <- generate_data(n0, mean_0, cov_matrix, label = 0, dist = dist_choice0)
  data_1 <- generate_data(n1, mean_1, cov_matrix, label = 1, dist = dist_choice1)
  
  # Combine data
  data_all <- rbind(data_0, data_1)
  
  # Ensure the feature columns are named X1 and X2 if not already
  if (!all(c("X1", "X2") %in% colnames(data_all))) {
    colnames(data_all)[1:k] <- paste0("X", 1:k)
  }
  
  # ----- Estimate 2D Density -----
  # Using MASS::kde2d to obtain a density estimate over a grid of size grid_n x grid_n
  density_estimate <- MASS::kde2d(data_all$X1, data_all$X2, n = grid_n)
  
  # Assign each data point the density value from the nearest grid point
  X3_density <- numeric(nrow(data_all))
  for (i in 1:nrow(data_all)) {
    idx_x <- which.min(abs(density_estimate$x - data_all$X1[i]))
    idx_y <- which.min(abs(density_estimate$y - data_all$X2[i]))
    X3_density[i] <- density_estimate$z[idx_x, idx_y]
  }
  data_all$X3 <- X3_density
  
  # ----- Create Interactive 3D Plot using plotly -----
  p <- plotly::plot_ly(data_all, 
                       x = ~X1, y = ~X2, z = ~X3, 
                       color = ~as.factor(label),
                       colors = c("red", "blue"),
                       type = "scatter3d", mode = "markers",
                       marker = list(size = 2)) %>%
    plotly::layout(title = "3D Scatter Plot with Density as Z",
                   scene = list(xaxis = list(title = "X1"),
                                yaxis = list(title = "X2"),
                                zaxis = list(title = "Density")))
  return(p)
}

# 5. Function: plot_density_glm
# Purpose: This function simulates data for two classes using the generate_data function, splits the combined dataset 
#          into training and testing subsets, fits a logistic regression model (GLM) on two selected features (X1 and X2),
#          and visualizes the predicted probabilities as overlaid density histograms for each class.
#
# Arguments:
#   mean0       - Numeric vector of length k specifying the means for each feature of Class 0.
#   mean1       - Numeric vector of length k specifying the means for each feature of Class 1.
#   cov_matrix  - A k x k covariance matrix used for generating features; its diagonal entries determine the variance of each feature.
#   dist_choice0- Character vector of length k indicating the distribution type for each feature for Class 0.
#   dist_choice1- Character vector of length k indicating the distribution type for each feature for Class 1.
#   n0          - Integer indicating the number of samples to generate for Class 0.
#   n1          - Integer indicating the number of samples to generate for Class 1.
#   k           - Integer representing the number of features in the dataset.
#   train_frac  - A numeric value between 0 and 1 specifying the fraction of the combined data to use for training.
#
# Returns:
#   A plotly figure object displaying overlaid density histograms of the predicted probabilities from the GLM for both classes:
#     - The histogram for "Class1" (red) shows the density of predicted probabilities for the positive class.
#     - The histogram for "Class0" (blue) shows the density of predicted probabilities for the negative class.
#
# Details:
#   - The function first generates synthetic data for both classes using the generate_data function.
#   - It assigns column names "X1" to "Xk" to the feature columns and converts the class labels into a factor with levels "Class0" and "Class1".
#   - The combined dataset is randomly split into training and testing subsets based on the train_frac parameter.
#   - A logistic regression model is then fitted using only the first two features (X1 and X2) from the training data.
#   - The model predicts probabilities on the testing set, which are then used to plot normalized density histograms using plotly.
#   - The resulting plot provides a visual assessment of the model's discrimination ability between the two classes.
#
plot_density_glm <- function(mean0, mean1, cov_matrix, dist_choice0, dist_choice1, n0, n1, k, train_frac) {
  # Generate synthetic data for each class
  data_0 <- generate_data(n0, mean0, cov_matrix, label = 0, dist = dist_choice0)
  data_1 <- generate_data(n1, mean1, cov_matrix, label = 1, dist = dist_choice1)
  data_full <- rbind(data_0, data_1)
  
  # Rename the first k columns to X1, X2, ..., Xk and convert the label to a factor with appropriate class names
  colnames(data_full)[1:k] <- paste0("X", 1:k)
  data_full$label <- factor(data_full$label, levels = c(0, 1), labels = c("Class0", "Class1"))
  
  # Split the data into training and testing sets based on the specified training fraction
  train_idx <- sample(seq_len(nrow(data_full)), size = round(train_frac * nrow(data_full)))
  train_data <- data_full[train_idx, ]
  test_data  <- data_full[-train_idx, ]
  
  # For the GLM, use only the first two features (X1 and X2) along with the class label
  train_glm <- train_data[, c("X1", "X2", "label")]
  test_glm  <- test_data[, c("X1", "X2", "label")]
  
  # Fit a logistic regression model on the training data and predict probabilities on the testing data
  model_glm <- glm(label ~ ., data = train_glm, family = binomial)
  test_glm$score_glm <- predict(model_glm, newdata = test_glm, type = "response")
  
  # Create an interactive density plot of the predicted probabilities using plotly
  fig <- plot_ly() %>%
    add_trace(x = test_glm$score_glm[test_glm$label == "Class1"], type = 'histogram',
              histnorm = 'density', name = 'Class1', opacity = 0.5, marker = list(color = 'red')) %>%
    add_trace(x = test_glm$score_glm[test_glm$label == "Class0"], type = 'histogram',
              histnorm = 'density', name = 'Class0', opacity = 0.5, marker = list(color = 'blue')) %>%
    layout(barmode = "overlay",
           title = "Predicted Probability (Score) Density Plot",
           xaxis = list(title = "Predicted Probability (Score)"),
           yaxis = list(title = "Density"))
  
  return(fig)
}


# 5. Function: run_power_analysis_return
# Purpose: This function performs a simulation-based power analysis by generating datasets for two classes,
#          training a GLM classifier, and evaluating the ability of two statistical tests (Mann–Whitney and KS test)
#          to detect differences in the predicted probabilities (scores) between the classes.
#
# Process:
#   - For each of M simulation iterations:
#       1. Generate data for Class 0 and Class 1 using the generate_data function, with specified means,
#          covariance matrix, and distribution choices.
#       2. Combine the datasets and assign column names for the features (X1, X2, ..., Xk) and the class label.
#       3. Randomly split the data into a training set and a test set based on the train_frac parameter.
#       4. Train a logistic regression (GLM) model on the training data using the features.
#       5. Obtain predicted probabilities on the test data.
#       6. Separate the predicted scores for the two classes and apply:
#              - A Mann–Whitney (Wilcoxon rank-sum) test.
#              - A Kolmogorov–Smirnov (KS) test.
#       7. Count a rejection if the p-value from a test is below the significance level alpha.
#       8. Also count a combined rejection when both tests reject the null hypothesis.
#
# Arguments:
#   exp_name     - A descriptive name for the experiment (character string).
#   mean0        - A numeric vector specifying the means for the features for Class 0.
#   mean1        - A numeric vector specifying the means for the features for Class 1.
#   cov_matrix   - A covariance matrix used in data generation (the diagonal provides variances for each feature).
#   dist_choice0 - A character vector specifying the distribution for each feature in Class 0.
#   dist_choice1 - A character vector specifying the distribution for each feature in Class 1.
#   n0           - The number of samples to generate for Class 0 (default is 200).
#   n1           - The number of samples to generate for Class 1 (default is 200).
#   k            - The number of features to consider in the analysis (default is 2).
#   train_frac   - The fraction of the full dataset to use for training (default is 0.7).
#   M            - The number of simulation iterations to perform (default is 100).
#   alpha        - The significance level for the hypothesis tests (default is 0.05).
#
# Returns:
#   A list containing:
#     - glm_mw:       The estimated power (proportion of rejections) using the Mann–Whitney test.
#     - glm_ks:       The estimated power (proportion of rejections) using the KS test.
#     - glm_combined: The estimated power when both tests reject the null hypothesis.
#
# Note:
#   This power analysis evaluates the sensitivity of the GLM-based classifier scores in differentiating
#   between the two classes through non-parametric testing.

run_power_analysis_return <- function(exp_name, mean0, mean1, cov_matrix, dist_choice0, dist_choice1,
                                      n0 = 200, n1 = 200, k = 2, train_frac = 0.7, M = 100, alpha = 0.05) {
  glm_mw_rejections       <- 0
  glm_ks_rejections       <- 0
  glm_combined_rejections <- 0
  
  for(m in 1:M) {
    # Generate data for Class 0 and Class 1 using the provided parameters
    data_0 <- generate_data(n0, mean0, cov_matrix, label = 0, dist = dist_choice0)
    data_1 <- generate_data(n1, mean1, cov_matrix, label = 1, dist = dist_choice1)
    
    # Merge the datasets and assign feature names
    data_full <- rbind(data_0, data_1)
    colnames(data_full)[1:k] <- paste0("X", 1:k)
    data_full$label <- factor(data_full$label, levels = c(0,1), labels = c("Class0", "Class1"))
    
    # Hold-Out Split: Randomly partition the data into training and test sets
    train_idx <- sample(seq_len(nrow(data_full)), size = round(train_frac * nrow(data_full)))
    train_data <- data_full[train_idx, ]
    test_data  <- data_full[-train_idx, ]
    
    # GLM Classification: Train the logistic regression model on the training data
    train_glm <- train_data[, c("X1", "X2", "label")]
    test_glm  <- test_data[, c("X1", "X2", "label")]
    model_glm <- glm(label ~ ., data = train_glm, family = binomial)
    
    # Predict probabilities on the test data
    test_glm$score_glm <- predict(model_glm, newdata = test_glm, type = "response")
    
    # Separate the predicted scores by class label
    S_plus_glm  <- test_glm$score_glm[test_glm$label == "Class1"]
    S_minus_glm <- test_glm$score_glm[test_glm$label == "Class0"]
    
    # Perform statistical tests on the predicted scores
    mw_test_glm <- wilcox.test(S_plus_glm, S_minus_glm, alternative = "two.sided")
    ks_test_glm <- ks.test(S_plus_glm, S_minus_glm)
    
    # Count rejections based on the significance level alpha
    if(mw_test_glm$p.value < alpha) {
      glm_mw_rejections <- glm_mw_rejections + 1
    }
    if(ks_test_glm$p.value < alpha) {
      glm_ks_rejections <- glm_ks_rejections + 1
    }
    if(mw_test_glm$p.value < alpha & ks_test_glm$p.value < alpha) {
      glm_combined_rejections <- glm_combined_rejections + 1
    }
  }
  
  # Compute power estimates for each test by averaging over the M simulations
  power_glm_mw       <- glm_mw_rejections / M
  power_glm_ks       <- glm_ks_rejections / M
  power_glm_combined <- glm_combined_rejections / M
  
  # Return a list with the power estimates for each testing procedure
  return(list(glm_mw = power_glm_mw,
              glm_ks = power_glm_ks,
              glm_combined = power_glm_combined))
}


# 6. Function: run_full_experiment
# Purpose: The run_experiment function orchestrates a simulation experiment to evaluate the statistical power
#          of different tests (MW, KS, and a combined test) across varying sample sizes for Class 1. For each
#          specified sample size (n1), it performs a power analysis by invoking an external function 
#          (run_power_analysis_return), aggregates the results, and then produces both interactive visualizations
#          and a data table summarizing the outcomes.
#
# Arguments:
#   exp_name      - A character string naming the experiment.
#   mean0         - A numeric vector specifying the mean values for the features in Class 0.
#   mean1         - A numeric vector specifying the mean values for the features in Class 1.
#   cov_matrix    - A covariance matrix for the features.
#   dist_choice0  - A character vector indicating the distribution choices for Class 0 features.
#   dist_choice1  - A character vector indicating the distribution choices for Class 1 features.
#   n0_fixed      - An integer specifying the fixed sample size for Class 0.
#   n1_values     - A numeric vector containing different sample sizes for Class 1 over which the experiment is run.
#   k             - An integer representing the number of features.
#   train_frac    - A numeric value (between 0 and 1) indicating the fraction of data to be used for training.
#   M_sim         - An integer representing the number of simulation iterations.
#   alpha         - A numeric value indicating the significance level for the hypothesis tests.
#
# Returns:
#   A list with the following components:
#     - power_results: A data frame summarizing the estimated power for the MW test, KS test, and combined test
#                      across different n1 values. It also includes a decision message based on the combined power.
#     - dt_table     : An interactive DT table displaying the power analysis results.
#     - fig_power    : A Plotly figure visualizing the relationship between n1 (sample size for Class 1) and the
#                      estimated power for each test.
#     - fig_density  : A Plotly figure displaying the predicted probability density plot based on the generalized
#                      linear model (GLM) for a selected n1 value.
#
# Details:
#   - The function iterates over each sample size in n1_values. For each iteration, it:
#       1. Prints progress information to the console.
#       2. Calls run_power_analysis_return to compute the power estimates for the three tests.
#       3. Stores the results in a data frame and makes a decision based on whether the combined power exceeds
#          a predefined threshold (0.05 in this example).
#       4. Outputs a message summarizing the decision for the current sample size.
#   - After processing all sample sizes, it generates:
#       * A Plotly plot (fig_power) to illustrate the estimated power against n1.
#       * A density plot (fig_density) of the predicted probabilities using plot_density_glm.
#       * An interactive DT table (dt_table) summarizing the power analysis results.
#
run_full_experiment <- function(exp_name,
                                mean0, mean1, cov_matrix,
                                dist_choice0, dist_choice1,
                                n0_fixed, n1_values,
                                k, train_frac,
                                M_sim, alpha=0.05) {
  # Create a data frame to store the power results
  power_results <- data.frame(n1 = n1_values,
                              glm_mw = NA,
                              glm_ks = NA,
                              glm_combined = NA,
                              stringsAsFactors = FALSE)
  
  # Loop over each n1 value and run the simulation
  for (i in seq_along(n1_values)) {
    cat(">> Processing n1 =", n1_values[i], "\n")
    
    # Call the power analysis function (assumes run_power_analysis_return is defined)
    res <- run_power_analysis_return(exp_name, mean0, mean1, cov_matrix, 
                                     dist_choice0, dist_choice1,
                                     n0 = n0_fixed, n1 = n1_values[i], k = k, 
                                     train_frac = train_frac, M = M_sim, alpha = alpha)
    
    # Store the results
    power_results[i, 2:4] <- unlist(res)
    
    # Insert the decision into the data frame based on combined power threshold (here 0.05 as an example)
    power_results$Decision[i] <- ifelse(res$glm_combined > 0.05,
                                        "Reject H0: Classes likely different",
                                        "Fail to reject H0: Insufficient evidence")
    
    # Also print a message to the console summarizing the decision
    if (res$glm_combined > 0.05) {
      cat("   Combined Power =", round(res$glm_combined, 3), 
          "- Conclusion: Reject H0. The classes are likely different.\n\n")
    } else {
      cat("   Combined Power =", round(res$glm_combined, 3), 
          "- Conclusion: Fail to reject H0. Insufficient evidence of difference.\n\n")
    }
  }
  
  cat("========== Completed", exp_name, "==========\n\n")
  
  # Generate the power analysis Plotly plot
  cat("Generating power analysis Plotly plot for", exp_name, "...\n")
  fig_power <- plot_ly(power_results, x = ~n1) %>%
    add_trace(y = ~glm_mw, type = 'scatter', mode = 'lines+markers', name = 'MW Power') %>%
    add_trace(y = ~glm_ks, type = 'scatter', mode = 'lines+markers', name = 'KS Power') %>%
    add_trace(y = ~glm_combined, type = 'scatter', mode = 'lines+markers', name = 'Combined Power') %>%
    layout(title = paste("Power vs n1 for", exp_name),
           xaxis = list(title = "n1 (Sample Size for Class 1)"),
           yaxis = list(title = "Estimated Power", range = c(0, 1)))
  
  # Generate the predicted probability density plot using a helper function
  cat("Generating predicted probability density plot for", exp_name, "...\n")
  fig_density <- plot_density_glm(mean0, mean1, cov_matrix, dist_choice0, dist_choice1,
                                  n0_fixed, n1_values[5], k, train_frac)
  
  # Create an interactive DT table from power_results
  cat("Generating DT table for", exp_name, "...\n")
  dt_table <- DT::datatable(power_results, 
                            options = list(pageLength = 5, autoWidth = TRUE),
                            caption = htmltools::tags$caption(
                              style = 'caption-side: top; text-align: center; color: #2E86AB; font-size: 16px;',
                              paste("Power Analysis Results for", exp_name)
                            ))
  
  # Return the results and plots as a list
  return(list(
    power_results = power_results,
    dt_table = dt_table,
    fig_power  = fig_power,
    fig_density = fig_density
  ))
}

# 6. Function: run_full_experiment
# Purpose: This unified function merges the functionalities of both run_experiment and run_power_experiment.
#          It conducts a simulation-based power analysis by iterating over various sample sizes for Class 1.
#          For each specified n1 value, it performs M_sim simulation iterations that:
#            - Generate data for Class 0 and Class 1 using a flexible data-generation function.
#            - Combine and split the data into training and testing sets.
#            - Fit a logistic regression model (GLM) and compute predicted probabilities.
#            - Perform the Mann–Whitney (Wilcoxon rank-sum) and Kolmogorov–Smirnov tests on the GLM scores.
#          The function then aggregates the rejection counts to estimate the statistical power for each test,
#          makes a decision based on a combined power threshold, and finally produces interactive visualizations
#          (a Plotly power plot and a predicted probability density plot) and a DT table summarizing the results.
#
# Arguments:
#   exp_name      - A character string naming the experiment.
#   mean0         - A numeric vector of feature means for Class 0.
#   mean1         - A numeric vector of feature means for Class 1.
#   cov_matrix    - A covariance matrix for the features.
#   dist_choice0  - A character vector indicating the distribution type for each feature in Class 0.
#   dist_choice1  - A character vector indicating the distribution type for each feature in Class 1.
#   n0_fixed      - An integer specifying the fixed sample size for Class 0.
#   n1_values     - A numeric vector of different sample sizes for Class 1.
#   k             - An integer representing the number of features.
#   train_frac    - A numeric value (0 to 1) specifying the fraction of data used for training.
#   M_sim         - An integer denoting the number of simulation iterations per n1 value.
#   alpha         - A numeric significance level used for the hypothesis tests.
#
# Returns:
#   A list containing:
#     - power_results: A data frame summarizing the estimated power for the Mann–Whitney test (MW Power),
#                      the Kolmogorov–Smirnov test (KS Power), and the combined test.
#     - dt_table     : An interactive DT table displaying the power analysis results.
#     - fig_power    : A Plotly figure showing power versus n1 for each test.
#     - fig_density  : A Plotly density plot of the predicted probabilities based on the GLM,
#                      generated using one example n1 value.
#
# Note: This function assumes that the helper functions 'generate_data' and 'plot_density_glm' are defined.
run_full_experiment <- function(exp_name,
                                mean0, mean1, cov_matrix,
                                dist_choice0, dist_choice1,
                                n0_fixed, n1_values,
                                k, train_frac,
                                M_sim, alpha) {
  # Create a data frame to store power analysis results for each n1 value
  power_results <- data.frame(n1 = n1_values,
                              glm_mw = NA,
                              glm_ks = NA,
                              glm_combined = NA,
                              Decision = NA,
                              stringsAsFactors = FALSE)
  
  # Loop over each sample size for Class 1
  for (i in seq_along(n1_values)) {
    cat(">> Processing n1 =", n1_values[i], "\n")
    
    # Initialize counters for rejections
    mw_rejections       <- 0
    ks_rejections       <- 0
    combined_rejections <- 0
    
    # Run M_sim simulation iterations for the current n1 value
    for (m in 1:M_sim) {
      # Generate data for Class 0 and Class 1
      data0 <- generate_data(n0_fixed, mean0, cov_matrix, label = 0, dist = dist_choice0)
      data1 <- generate_data(n1_values[i], mean1, cov_matrix, label = 1, dist = dist_choice1)
      data_full <- rbind(data0, data1)
      
      # Rename first k columns to X1, X2, …, Xk and set factor labels for the class
      colnames(data_full)[1:k] <- paste0("X", 1:k)
      data_full$label <- factor(data_full$label, levels = c(0, 1), labels = c("Class0", "Class1"))
      
      # Randomly split into training and testing sets based on train_frac
      train_idx <- sample(seq_len(nrow(data_full)), size = round(train_frac * nrow(data_full)))
      train_data <- data_full[train_idx, ]
      test_data  <- data_full[-train_idx, ]
      
      # Fit a logistic regression (GLM) model on the training data
      model_glm <- glm(label ~ ., data = train_data[, c(paste0("X", 1:k), "label")], family = binomial)
      
      # Predict probabilities on the test data
      test_data$score_glm <- predict(model_glm, newdata = test_data, type = "response")
      
      # Separate the predicted scores by class
      scores_class1 <- test_data$score_glm[test_data$label == "Class1"]
      scores_class0 <- test_data$score_glm[test_data$label == "Class0"]
      
      # Apply the Mann–Whitney (Wilcoxon rank-sum) and Kolmogorov–Smirnov tests
      mw_test <- wilcox.test(scores_class1, scores_class0, alternative = "two.sided")
      ks_test <- ks.test(scores_class1, scores_class0)
      
      # Increment rejection counters if p-values are below the significance level
      if (mw_test$p.value < alpha) mw_rejections <- mw_rejections + 1
      if (ks_test$p.value < alpha) ks_rejections <- ks_rejections + 1
      if (mw_test$p.value < alpha && ks_test$p.value < alpha) {
        combined_rejections <- combined_rejections + 1
      }
    } # End simulation iterations for current n1 value
    
    # Compute power estimates as the proportion of rejections
    power_results$glm_mw[i]       <- mw_rejections / M_sim
    power_results$glm_ks[i]       <- ks_rejections / M_sim
    power_results$glm_combined[i] <- combined_rejections / M_sim
    
    # Make a decision based on the combined power (threshold here is 0.05)
    power_results$Decision[i] <- ifelse(power_results$glm_combined[i] > 0.05,
                                        "Reject H0: Classes likely different",
                                        "Fail to reject H0: Insufficient evidence")
    
    cat("\n----------------------------------------------------------\n")
    cat("Completed Power Analysis for", exp_name, "\n")
    cat("GLM Test Powers:\n")
    cat("   Mann–Whitney Power:", round(power_results$glm_mw[i], 3), "\n")
    cat("   KS Test Power    :", round(power_results$glm_ks[i], 3), "\n")
    cat("----------------------------------------------------------\n\n")
    # Print a summary message for the current n1 value
    if (power_results$glm_combined[i] > 0.05) {
      cat("   Combined Power =", round(power_results$glm_combined[i], 3), 
          "- Conclusion: Reject H0. The classes are likely different.\n\n")
    } else {
      cat("   Combined Power =", round(power_results$glm_combined[i], 3), 
          "- Conclusion: Fail to reject H0. Insufficient evidence of difference.\n\n")
    }
  } # End loop over n1 values
  
  cat("========== Completed", exp_name, " ==========\n\n")
  
  # Generate the power analysis Plotly plot
  fig_power <- plot_ly(power_results, x = ~n1) %>%
    add_trace(y = ~glm_mw, type = 'scatter', mode = 'lines+markers', name = 'MW Power') %>%
    add_trace(y = ~glm_ks, type = 'scatter', mode = 'lines+markers', name = 'KS Power') %>%
    add_trace(y = ~glm_combined, type = 'scatter', mode = 'lines+markers', name = 'Combined Power') %>%
    layout(title = paste("Power vs n1 for", exp_name),
           xaxis = list(title = "n1 (Sample Size for Class 1)"),
           yaxis = list(title = "Estimated Power", range = c(0, 1)))
  
  # --- Generate a density plot for the GLM-predicted probabilities ---
  # For the density plot, choose one n1 value (e.g., the 5th element if available, otherwise the first)
  density_index <- ifelse(length(n1_values) >= 5, 5, 1)
  n1_density <- n1_values[density_index]
  
  # Generate synthetic data for the density plot
  data0 <- generate_data(n0_fixed, mean0, cov_matrix, label = 0, dist = dist_choice0)
  data1 <- generate_data(n1_density, mean1, cov_matrix, label = 1, dist = dist_choice1)
  data_full <- rbind(data0, data1)
  colnames(data_full)[1:k] <- paste0("X", 1:k)
  data_full$label <- factor(data_full$label, levels = c(0, 1), labels = c("Class0", "Class1"))
  
  # Split the data into training and testing sets for the density plot
  train_idx <- sample(seq_len(nrow(data_full)), size = round(train_frac * nrow(data_full)))
  train_data <- data_full[train_idx, ]
  test_data  <- data_full[-train_idx, ]
  
  # Fit the GLM on the training set and predict probabilities on the test set
  model_glm <- glm(label ~ ., data = train_data[, c(paste0("X", 1:k), "label")], family = binomial)
  test_data$score_glm <- predict(model_glm, newdata = test_data, type = "response")
  
  # Create an interactive density plot of the predicted probabilities using Plotly
  fig_density <- plot_ly() %>%
    add_trace(x = test_data$score_glm[test_data$label == "Class1"],
              type = 'histogram', histnorm = 'density',
              name = 'Class1', opacity = 0.5, marker = list(color = 'red')) %>%
    add_trace(x = test_data$score_glm[test_data$label == "Class0"],
              type = 'histogram', histnorm = 'density',
              name = 'Class0', opacity = 0.5, marker = list(color = 'blue')) %>%
    layout(barmode = "overlay",
           title = "Predicted Probability (Score) Density Plot",
           xaxis = list(title = "Predicted Probability (Score)"),
           yaxis = list(title = "Density"))
  
  # Return a list containing the power results, DT table, and both plots
  return(list(
    power_results = power_results,
    #dt_table = dt_table,
    fig_power = fig_power,
    fig_density = fig_density
  ))
}