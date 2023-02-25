#' Forest In-Sample Honest Weights
#'
#' Computes forest in-sample honest weights for an \code{morf.forest} object.
#'
#' @param forest An \code{morf.forest} object.
#' @param train_sample Training sample.
#' @param honest_sample Honest sample. 
#' 
#' @return 
#' Matrix of in-sample honest weights.
#' 
#' @keywords internal
#' 
#' @details 
#' \code{forest} must have been grown using only the training sample. 
#' 
#' @importFrom stats ave
forest_weights_fitted <- function(forest, honest_sample, train_sample) { # Taken from https://github.com/okasag/orf/blob/master/orf/R/weight_funs.R
  ## Handling inputs.
  # Get terminal nodes for the honest sample.
  leaf_IDs_train <- predict(forest, train_sample, type = "terminalNodes")$predictions
  leaf_IDs_train <- lapply(seq_along(leaf_IDs_train[1, ]), function(i) leaf_IDs_train[, i])
  
  leaf_IDs_honest <- predict(forest, honest_sample, type = "terminalNodes")$predictions
  leaf_IDs_honest <- lapply(seq_along(leaf_IDs_honest[1, ]), function(i) leaf_IDs_honest[, i])
  
  # Compute leaf size using honest units.
  leaf_size_honest <- lapply(leaf_IDs_honest, function(x) ave(x, x, FUN = length))
  
  ## Compute weights for the whole sample. Notice that the output matrix stores first honest and then train units (row-wise).
  forest_weights <- forest_weights_fitted_cpp(leaf_IDs_train, leaf_IDs_honest, leaf_size_honest)
  
  ## Order according to original sample (ascending rownames).
  rownames(forest_weights) <- c(rownames(honest_sample), rownames(train_sample))
  forest_weights <- as.matrix(forest_weights[order(as.numeric(row.names(forest_weights))), ])
  
  ## Output.
  return(forest_weights)
}


#' Forest Out-of-Sample Weights
#'
#' Computes forest out-of-sample honest weights for an \code{morf.forest} object.
#'
#' @param forest An \code{morf.forest} object.
#' @param test_sample Test sample.
#' @param honest_sample Honest sample. 
#'
#' @return 
#' Matrix of out-of-sample honest weights.
#' 
#' @keywords internal
#'
#' @details 
#' \code{forest} must have been grown using only the training sample. 
#' 
#' @importFrom stats ave
predict_forest_weights <- function(forest, honest_sample, test_sample) { # Taken from https://github.com/okasag/orf/blob/master/orf/R/weight_funs.R
  ## Get terminal nodes for the honest sample.
  leaf_IDs_honest <- predict(forest, honest_sample, type = "terminalNodes")$predictions
  leaf_IDs_honest <- lapply(seq_along(leaf_IDs_honest[1, ]), function(i) leaf_IDs_honest[, i])
  
  leaf_IDs_test <- predict(forest, test_sample, type = "terminalNodes")$predictions
  leaf_IDs_test <- lapply(seq_along(leaf_IDs_test[1, ]), function(i) leaf_IDs_test[, i])
  
  ## Compute leaf size using honest units.
  leaf_size_honest <- lapply(leaf_IDs_honest, function(x) ave(x, x, FUN = length))
  
  ## Compute weights for the test sample.
  forest_weights <- forest_weights_predicted_cpp(leaf_IDs_test, leaf_IDs_honest, leaf_size_honest, 0)
  
  ## Output.
  return(forest_weights)
}
