# 7. Function: summary_function
# Purpose: This function calculates summary statistics for a numeric time series.
#
# Arguments:
#   time_series - A numeric vector representing the time series data.
#
# Returns:
#   A named vector containing the following summary statistics:
#     - mean: The mean of the time series.
#     - sd: The standard deviation of the time series.
#     - slope: The mean of the differences between consecutive time points.
#     - sd_slope: The standard deviation of the differences between consecutive time points.
#     - min: The minimum value in the time series.
#     - max: The maximum value in the time series.
summary_function <- function(time_series) {
  x <- as.numeric(time_series)  # Ensure numeric input
  
  # Calculate the global features
  results <- c(
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    slope = mean(diff(x), na.rm = TRUE),
    sd_slope = sd(diff(x), na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
  
  return(results)
}

# 8. Function: summary_dataset
# Purpose: This function computes summary statistics for speed and altitude features in a dataset and combines them.
#
# Arguments:
#   df - A dataframe containing speed and altitude features, along with a 'HR_Zone' label column.
#
# Returns:
#   A dataframe containing the computed summary statistics for speed and altitude, with rows having missing values removed.
summary_dataset <- function(df){
  
  # Custom Function for hw dataset
  
  speed_stat <- t(apply(df[,2:61], 1, summary_function))
  altitude_stat <- t(apply(df[,62:121], 1, summary_function))
  
  # Add the "v" and "a" to the names of the columns to make them recognizable
  colnames(speed_stat) <- paste("v", colnames(speed_stat), sep = "_")
  colnames(altitude_stat) <- paste("a", colnames(altitude_stat), sep = "_")
  
  # Combine the statistics of speed and altitude 
  df_stats <- cbind(speed_stat, altitude_stat)
  
  # Add the labels
  df_stats$HR_Zone <- df$HR_Zone
  
  #Remove the rows with missing values
  df_stats <- df_stats[complete.cases(df_stats),]
  df_stats <- as.data.frame(df_stats)
  
  return(df_stats)
}

# 9. Function: find_uncorrelated_features
# Purpose: This function identifies and prints pairs of features in a dataset that are highly correlated.
#
# Arguments:
#   features_df - A dataframe containing the features to analyze.
#   correlation_threshold - A numeric value specifying the threshold for considering correlations as high (default is 0.8).
#
# Returns:
#   None. The function prints the highly correlated feature pairs and their correlation values.
#
find_uncorrelated_features <- function(features_df, correlation_threshold = 0.8) {
  # Compute correlation matrix
  cor_matrix <- cor(features_df)
  
  # Extract correlated pairs
  cor_pairs <- which(abs(cor_matrix) > correlation_threshold, arr.ind = TRUE)
  
  # Remove duplicates (only keep upper triangle of the correlation matrix)
  cor_pairs <- cor_pairs[cor_pairs[, 1] < cor_pairs[, 2], ]
  
  # Print correlated pairs with their correlation values
  cat("Highly correlated variable pairs (correlation > ", correlation_threshold, "):\n")
  for (i in 1:nrow(cor_pairs)) {
    var1 <- colnames(features_df)[cor_pairs[i, 1]]
    var2 <- colnames(features_df)[cor_pairs[i, 2]]
    cor_value <- cor_matrix[cor_pairs[i, 1], cor_pairs[i, 2]]
    cat(var1, "<->", var2, ":", round(cor_value, 4), "\n")
  }
}

# 10. Function: hold_out_evaluation_RF
# Purpose: This function performs hold-out evaluation of a Random Forest model and computes various performance metrics.
#
# Arguments:
#   df - A dataframe containing the dataset with a 'label' column indicating class membership.
#   ntree - An integer specifying the number of trees to grow in the Random Forest model (default is 500).
#   seed - An integer to set the seed for reproducibility (default is 123).
#
# Returns:
#   A list containing the following performance metrics:
#     - Accuracy: The accuracy of the model on the test set.
#     - ROC_AUC: The area under the ROC curve.
#     - Precision: The precision of the model for the positive class.
#     - Recall: The recall of the model for the positive class.
#     - F1_Score: The F1 score of the model for the positive class.
hold_out_evaluation_RF <- function(df, ntree = 500, seed = 123){
  
  if (!"label" %in% colnames(df)) {
    stop("Error: The dataframe doesn't have a label column.")
  }
  
  df$label <- factor(df$label, ordered = FALSE)
  
  set.seed(seed)
  
  # Split into training and test sets
  train_index <- sample(1:nrow(df), 0.7*nrow(df))
  train_data <- as.data.frame(df[train_index,])
  test_data <- as.data.frame(df[-train_index,])
  
  # Train Random Forest model
  model <- randomForest(label ~ ., data = train_data, ntree = ntree)
  
  # Predict probability for ROC curve
  scores <- predict(model, newdata = test_data, type = "prob")
  
  # Predict class labels
  label_score <- predict(model, newdata = test_data, type = "class")
  
  # ROC AUC
  roc_result <- roc(test_data$label, scores[,2])
  roc_auc <- auc(roc_result)
  
  # Accuracy
  acc <- mean(test_data$label == label_score) * 100
  
  # Confusion Matrix
  cm <- table(Predicted = label_score, Actual = test_data$label)
  
  # Precision, Recall, F1 Score (See the label 1 as Positive Label)
  precision <- cm[2,2] / (cm[2,2] + cm[1,2])
  recall <- cm[2,2] / (cm[2,2] + cm[2,1])
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  # Return the metrics
  return(list(
    Accuracy = acc,
    ROC_AUC = roc_auc,
    Precision = precision,
    Recall = recall,
    F1_Score = f1_score))
}

# 11. Function: perm_test_pvalue
# Purpose: This function performs a permutation test to compute the empirical p-value for comparing
#          the distributions of scores obtained from a Random Forest classifier. It evaluates whether
#          the distributions of scores for two classes are significantly different.
#
# Arguments:
#   df - A dataframe containing the dataset with a 'label' column indicating class membership.
#   ntree - An integer specifying the number of trees to grow in the Random Forest model (default is 500).
#   seed - An integer to set the seed for reproducibility (default is 123).
#   file_name - A string specifying the file name to save the results of the permutation test (default is "perm_test.RData").
#
# Returns:
#   A list containing:
#     - zero_scores: The predicted scores for instances labeled as 0.
#     - one_scores: The predicted scores for instances labeled as 1.
#     - test_statistic: The test statistic from the chosen statistical test.
#     - p_value: The empirical p-value computed from the permutation test.
#     - null_distribution: The null distribution of test statistics obtained from permutations.
perm_test_pvalue <- function(df, ntree = 500, seed = 123, file_name = "permutations/perm_test.RData"){
  
  if (!"label" %in% colnames(df)) {
    stop("Error: The dataframe doesn't have a label column.")
  }
  
  df$label <- factor(df$label, ordered = FALSE)
  
  set.seed(seed)
  
  model <- randomForest(label ~ ., data = df, ntree = ntree)
  
  scores <- predict(model, newdata = df, type = "prob")
  
  # Divide the score for each zone
  zero_scores <- scores[,1][df$label == 0]
  one_scores <- scores[,2][df$label == 1]
  
  test_pvalue <- wilcox.test(one_scores, zero_scores)
  
  test_statistic <- test_pvalue$statistic
  
  num_permutations <- 1000
  null_distribution <- numeric(num_permutations)
  
  # Set up parallel processing
  cl <- makeCluster(detectCores() - 1) # Leave one core free
  registerDoParallel(cl)
    
  null_distribution <- foreach(i = seq_len(num_permutations), .combine = 'c', .packages = c('randomForest')) %dopar% {
    
    y_permuted <- sample(df$label)
    
    permuted_model <- randomForest(y_permuted ~ ., data = df, ntree = ntree)
    
    permuted_scores <- predict(permuted_model, newdata = df, type = "prob")
    
    
    stat <- wilcox.test(permuted_scores[,1][y_permuted == 0], permuted_scores[,2][y_permuted == 1])$statistic
    
    stat
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # Compute empiric p.value
  p_value <- mean(null_distribution >= test_statistic)
  
  # Save offline in order to not repeat the permutation test and to allow retrieve results in a second moment
  save(zero_scores, one_scores, test_statistic, p_value, null_distribution, file = file_name)
  
  return(list(zero_scores, one_scores,
              test_statistic, p_value, null_distribution))
  
}