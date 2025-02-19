# Statistical Classification and Hypothesis Testing: Applied Methods and Analysis

## Overview

This repository contains the homework assignments for the course Stat4DS, completed by Andrea De Vincenzo, Domenico Azzarito, and Michele Pezza. The assignments cover various statistical and machine learning concepts, including classification, hypothesis testing, and data analysis.

## Part 1: Bayes Classification Rule

### Objective

In this part, we derive and understand the Bayes classification rule, denoted as η∗(x), in a binary classification scenario. The goal is to minimize classification errors by assigning each observation to the class with the higher posterior probability.

### Key Concepts

- **Bayes Classification Rule**: Derived using conditional probabilities and the regression function.
- **Simulation**: Generated a dataset of size n=1000 from the joint data model and evaluated the performance of the Bayes classifier.
- **Comparison**: Compared the Bayes classifier with a logistic regression classifier.

### Results

- **Bayes Classifier**: Achieved an accuracy of 73.3% on the simulated dataset.
- **Logistic Regression**: Achieved an accuracy of 75.8% on the same dataset.
- **Repeated Sampling**: Conducted repeated sampling to evaluate the generalization performance of both classifiers.

## Part 2: Friedman's Procedure for Two-Sample Testing

### Objective

Summarize and implement the methodology described by Jerome H. Friedman for multivariate goodness-of-fit and two-sample testing using binary classifiers.

### Key Concepts

- **Pseudo-code**: Provided a high-level overview of Friedman's procedure.
- **Classifiers**: Used generalized linear models (GLM) for classification.
- **Experiments**: Conducted experiments to evaluate the power of the testing procedure under different scenarios.

### Results

- **Power Analysis**: Compared the performance of Kolmogorov-Smirnov (KS) and Mann-Whitney-Wilcoxon (MWW) tests.
- **Distributional Differences**: Identified scenarios where both tests consistently detect significant differences.
- **Heavy-Tailed Distributions**: Observed the challenges in detecting differences in heavy-tailed distributions.

## Part 3: Application to Heart Rate (HR) Data

### Objective

Apply Friedman's procedure to a dataset of heart rate (HR) zones in running, using embedded features and a Random Forest classifier.

### Key Concepts

- **Feature Embedding**: Reduced dimensionality by applying statistical functionals to time series data.
- **Classification**: Used a Random Forest classifier to capture non-linear relationships.
- **Hypothesis Testing**: Performed permutation tests to evaluate the significance of distributional differences.

### Results

- **Feature Selection**: Identified highly correlated features and evaluated their impact on classification performance.
- **Classification Performance**: Achieved stable performance with both 100 and 1000 trees in the Random Forest classifier.
- **Hypothesis Testing**: Strongly rejected the null hypothesis, indicating significant differences in feature distributions across HR zones.

## Repository Structure

- **src/**:
  - Includes R Markdown file of the homework.
  - Contains CSS file for styling the rendered HTML.
  - Include the HR dataset in RData
  - **module/**: Contains three R files with custom functions, each corresponding to a part of the homework.
  - **permutations/**: Contains RData files with the results of the permutation tests to avoid recomputation.
- **Homework-SDS.html**: The main HTML file containing plots, results, and analysis.

## Contributors

- **Andrea De Vincenzo**: [GitHub Profile](https://github.com/andreadv2000)
- **Domenico Azzarito**: [GitHub Profile](https://github.com/azzadom)
- **Michele Pezza**: [GitHub Profile](https://github.com/michelepezza99)

## References

- Brutti, P. (2025). Stat4DS | THE Homework | Part 1, 2+3.
- Friedman, J. H. (2003). On Multivariate Goodness–of–Fit and Two–Sample Testing.

## License

This project is licensed under the MIT License.
