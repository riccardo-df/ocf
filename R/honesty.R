##' Honest Sample Split (Internal Use)
##'
##' Randomly spits the sample into a training sample and an honest sample.
##'
##' @param data \code{data.frame} or \code{matrix} to be split. The outcome must be located in the first column.
##' @param honesty.fraction Fraction of honest sample.
##' 
##' @details 
##' \code{class_honest_split} looks for balanced splits, i.e., splits such as all the outcome's classes are represented
##' in both the training and the honest sample. After 100 trials, the program throws an error.
##'
##' @return 
##' List with elements:
##'   \item{\code{train_sample}}{Training sample.}
##'   \item{\code{honest_sample}}{Honest sample.}
class_honest_split <- function(data, honesty.fraction = 0.5) { # Inspired by https://github.com/okasag/orf/blob/master/orf/R/honest_funs.R
  ## Handling inputs.
  n <- nrow(data)
  y <- data[, 1]
  
  size <- floor(n * honesty.fraction)
  
  condition <- TRUE
  i <- 0
  
  ## Try new splits until we have a balanced split or we tried too many times.
  while (condition) {
    train_sample_idx <- sample(1:n, size = size, replace = FALSE)
  
    i <- i + 1
    
    if(all(unique(y[train_sample_idx]) %in% unique(y[-train_sample_idx])) | i == 100) break
  }
  
  if(!all(unique(y[train_sample_idx]) %in% unique(y[-train_sample_idx]))) stop("At least one of the classes of Y contains too few observations, preventing balanced honest splits. Consider recoding your outcome into less categories or setting honesty = FALSE.", call. = FALSE)
  
  # Reordering to avoid rownames-related clashes.
  train_sample_idx <- as.integer(sort(train_sample_idx, decreasing = FALSE))
  
  ## Handling output.
  train_sample <- data[train_sample_idx, ]
  honest_sample <- data[-train_sample_idx, ]
  
  ## Output.
  return(list("train_sample" = train_sample,
              "honest_sample" = honest_sample))
}


##' Honest In-Sample Predictions (Internal Use)
##'
##' Computes honest in-sample predictions for a \code{morf.forest} object relative to the m-th class.
##'
##' @param forest \code{morf.forest} object.
##' @param train_sample Training sample.
##' @param honest_sample Honest sample. 
##' @param m The class for which predictions are desired.
##' 
##' @details 
##' \code{forest} must have been grown using only the training sample. \code{honest_fitted} replaces the leaf estimates 
##' using the outcome from the honest sample (using the prediction method of \code{\link{morf}}).
##' 
##' @return 
##' In-sample honest predictions.
honest_fitted <- function(forest, train_sample, honest_sample, m) { # Inspired by https://github.com/okasag/orf/blob/master/orf/R/honest_funs.R
  ## Handling inputs.
  # Getting terminal nodes for the training and honest sample.
  train_leaves  <- predict(forest, train_sample, type = "terminalNodes")$predictions
  honest_leaves <- predict(forest, honest_sample, type = "terminalNodes")$predictions
  
  # Unique leaves for each tree.
  unique_leaves_train <- apply(train_leaves, 2, unique)
  unique_leaves_honest <- apply(honest_leaves, 2, unique)
  
  # Honest outcome indicators for the m-th class.
  honest_y_m <- ifelse(honest_sample[, 1] <= m, 1, 0)
  honest_y_m_1 <- ifelse(honest_sample[, 1] <= m - 1, 1, 0)
  
  ## Computing honest predictions.
  honest_fitted_values <- as.matrix(honest_fitted_cpp(unique_leaves_honest, honest_y_m, honest_y_m_1, honest_leaves, train_leaves))
  
  # Combining in dataset (first honest rownames, then train rownames).
  rownames(honest_fitted_values) <- c(rownames(honest_sample), rownames(train_sample))
  forest_fitted_values <- honest_fitted_values[order(as.numeric(row.names(honest_fitted_values))), ]
  
  ## Output
  return(as.numeric(forest_fitted_values)) 
}


##' Honest Out-of-Sample Predictions
##'
##' Computes honest out-of-sample predictions for a \code{morf.forest} object relative to the m-th class.
##'
##' @param forest \code{morf.forest} object.
##' @param test_sample Test sample.
##' @param honest_sample Honest sample. 
##' @param m The class for which predictions are desired.
##'
##' @details 
##' \code{honest_predictions} replaces the leaf estimates of \code{forest} using the outcome from the associated 
##' honest sample (using the prediction method of \code{\link{morf}}). The honest sample must not have been used
##' to build the trees.
##'
##' @return 
##' Out-of-sample honest predictions.
honest_predictions <- function(forest, honest_sample, test_sample, m) {
  ## Handling inputs.
  # Getting terminal nodes for the honest and the test sample.
  honest_leaves <- predict(forest, honest_sample, type = "terminalNodes")$predictions
  test_leaves <- predict(forest, test_sample, type = "terminalNodes")$predictions
  
  # Unique leaves for each tree.
  unique_leaves_honest <- apply(honest_leaves, 2, unique)

  # Honest outcome indicators for the m-th class.
  honest_y_m <- ifelse(honest_sample[, 1] <= m, 1, 0)
  honest_y_m_1 <- ifelse(honest_sample[, 1] <= m - 1, 1, 0)
  
  ## Computing honest predictions.
  honest_predictions <- honest_predictions_cpp(unique_leaves_honest, honest_y_m, honest_y_m_1, honest_leaves, test_leaves)

  ## Output.
  return(as.numeric(honest_predictions))
}
