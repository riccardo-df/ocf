#' Honest Sample Split
#'
#' Randomly spits the sample into a training sample and an honest sample.
#'
#' @param data \code{data.frame} or \code{matrix} to be split. The outcome must be located in the first column.
#' @param honesty.fraction Fraction of honest sample.
#' 
#' @return 
#' List with elements:
#'   \item{\code{train_sample}}{Training sample.}
#'   \item{\code{honest_sample}}{Honest sample.}
#'   
#' @keywords internal
#' 
#' @details 
#' \code{class_honest_split} looks for balanced splits, i.e., splits such as all the outcome's classes are represented
#' in both the training and the honest sample. After 100 trials, the program throws an error.
class_honest_split <- function(data, honesty.fraction = 0.5) { # Inspired by https://github.com/okasag/orf/blob/master/orf/R/honest_funs.R
  ## Handling inputs.
  n <- nrow(data)
  y <- data[, 1]
  size <- floor(n * honesty.fraction)
  counter <- 0
  
  ## Sample split.
  while (TRUE) {
    train_sample_idx <- sample(1:n, size = size, replace = FALSE)
    counter <- counter + 1
    if(all(unique(y[train_sample_idx]) %in% unique(y[-train_sample_idx]))) break
    if(counter == 100) stop("Cannot find balanced splits. Maybe one of the classes contains too few observations. Consider recoding your outcome into less categories or setting 'honesty = FALSE'.", call. = FALSE)
  }

  ## Handle output.
  train_sample_idx <- as.integer(sort(train_sample_idx, decreasing = FALSE)) # Reordering to avoid rownames-related clashes.
  train_sample <- data[train_sample_idx, ]
  honest_sample <- data[-train_sample_idx, ]
  
  ## Output.
  return(list("train_sample" = train_sample,
              "honest_sample" = honest_sample))
}


#' Honest In-Sample Predictions
#'
#' Computes honest in-sample predictions for an \code{morf.forest} object.
#'
#' @param forest An \code{morf.forest} object.
#' @param train_sample Training sample.
#' @param honest_sample Honest sample. 
#' @param y_m_honest Indicator variable, whether the outcome is smaller than or equal to the m-th class.
#' @param y_m_1_honest Indicator variable, whether the outcome is smaller than or equal to the (m-1)-th class.
#'
#' @return 
#' In-sample honest predictions.
#' 
#' @keywords internal
#'
#' @details 
#' \code{forest} must have been grown using only the training sample. \code{honest_fitted} replaces the leaf estimates 
#' using the outcome from the honest sample (using the prediction method of \code{\link{morf}}).
honest_fitted <- function(forest, train_sample, honest_sample, y_m_honest, y_m_1_honest) { # Inspired by https://github.com/okasag/orf/blob/master/orf/R/honest_funs.R
  ## Handling inputs.
  # Getting terminal nodes for the training and honest sample.
  train_leaves  <- predict(forest, train_sample, type = "terminalNodes")$predictions
  honest_leaves <- predict(forest, honest_sample, type = "terminalNodes")$predictions
  
  # Unique leaves for each tree.
  unique_leaves_train <- apply(train_leaves, 2, unique)
  unique_leaves_honest <- apply(honest_leaves, 2, unique)
  
  ## Compute honest predictions. Notice that the output matrix stores first honest and then train units (row-wise).
  honest_fitted_values <- as.matrix(honest_fitted_cpp(unique_leaves_honest, y_m_honest, y_m_1_honest, honest_leaves, train_leaves))
  
  ## Combine.
  rownames(honest_fitted_values) <- c(rownames(honest_sample), rownames(train_sample))
  honest_fitted_values <- honest_fitted_values[order(as.numeric(row.names(honest_fitted_values))), ]
  
  ## Output
  return(as.numeric(honest_fitted_values)) 
}


#' Honest Out-of-Sample Predictions
#'
#' Computes honest out-of-sample predictions for an \code{morf.forest} object.
#'
#' @param forest \code{morf.forest} object.
#' @param test_sample Test sample.
#' @param honest_sample Honest sample. 
#' @param y_m_honest Indicator variable, whether the outcome is smaller than or equal to the m-th class.
#' @param y_m_1_honest Indicator variable, whether the outcome is smaller than or equal to the (m-1)-th class.
#'
#' @return 
#' Out-of-sample honest predictions.
#' 
#' @keywords internal
#'
#' @details 
#' \code{honest_predictions} replaces the leaf estimates of \code{forest} using the outcome from the associated 
#' honest sample (using the prediction method of \code{\link{morf}}). The honest sample must not have been used
#' to build the trees.
honest_predictions <- function(forest, honest_sample, test_sample, y_m_honest, y_m_1_honest) {
  ## Handling inputs.
  # Getting terminal nodes for the honest and the test sample.
  honest_leaves <- predict(forest, honest_sample, type = "terminalNodes")$predictions
  test_leaves <- predict(forest, test_sample, type = "terminalNodes")$predictions
  
  # Unique leaves for each tree.
  unique_leaves_honest <- apply(honest_leaves, 2, unique)
  
  ## Computing honest predictions.
  honest_predictions <- honest_predictions_cpp(unique_leaves_honest, y_m_honest, y_m_1_honest, honest_leaves, test_leaves)

  ## Output.
  return(as.numeric(honest_predictions))
}
