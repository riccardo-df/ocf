#' Mean Squared Error
#'
#' Computes the mean squared error for evaluating the accuracy of ordered probability predictions.
#'
#' @param predictions Matrix of predictions (\code{n.samples} x \code{M}).
#' @param y Either the observed outcome vector or a matrix of true probabilities (\code{n.samples} x \code{M}).
#' @param use.true If \code{FALSE}, then the program assumes that \code{y} stores the observed outcome vector, otherwise it treats \code{y} as a matrix of true probabilities.
#'
#' @details 
#' If \code{use.true = FALSE}, the mean squared error is computed as follows:
#' 
#' \deqn{\frac{1}{n} \sum_{i = 1}^n \sum_{m = 1}^M (1 (Y_i = m) - \hat{p}_m (x))^2}
#' 
#' otherwise:
#' 
#' \deqn{\frac{1}{n} \sum_{i = 1}^n \sum_{m = 1}^M (p_m (x) - \hat{p}_m (x))^2}
#' 
#' where:
#' 
#' \deqn{p_m (x) = P(Y_i = m \, | \, X_i = x)}
#' 
#' The second formula is useful for simulation studies.
#'
#' @return The mean squared error of the method.
#' 
#' @export
mean_squared_error <- function(y, predictions, use.true = FALSE) { # Taken from https://github.com/okasag/orf/blob/master/orf/R/evals.R
  ## Handling inputs and checks.
  if (!(use.true %in% c(TRUE, FALSE))) stop("'use.true' must be logical.", call. = FALSE)
  if (use.true == TRUE & is.null(dim(y)) | use.true == FALSE & !is.null(dim(y))) stop("Combination of 'use.true' and 'y' is non-sensical.")
  n <- dim(predictions)[1]
  n.classes <- dim(predictions)[2]
  
  ## Generating matrix for comparison. Either a matrix of indicator variables, or the true probability matrix.
  if (use.true == FALSE) {
    indicator_mat <- matrix(0, nrow = n, ncol = n.classes)
    for (i in seq_len(n)) {
      indicator_mat[i, y[i]] <- 1
    }
  } else if (use.true == TRUE) {
    indicator_mat <- y
  }

  ## MSE.
  mse <- mean(rowSums(apply((indicator_mat - predictions), 2, function(x) {x^2})))
  
  ## Output.
  return(mse)
}


#' Mean Ranked Probability Score
#'
#' Computes the mean ranked probability score for evaluating the accuracy of ordered probability predictions.
#'
#' @param predictions Matrix of predictions (\code{n.samples} x \code{M}).
#' @param y Either the observed outcome vector or a matrix of true probabilities (\code{n.samples} x \code{M}). 
#' @param use.true If \code{FALSE} (the default), then the program assumes that \code{y} stores the observed outcome vector, otherwise it treats \code{y} as a matrix of true probabilities.
#'
#' @details 
#' If \code{use.true = FALSE}, the mean ranked probability score is computed as follows:
#' 
#' \deqn{\frac{1}{n} \sum_{i = 1}^n \frac{1}{M - 1} \sum_{m = 1}^M (1 (Y_i \leq m) - \hat{p}_m^* (x))^2}
#' 
#' otherwise:
#' 
#' \deqn{\frac{1}{n} \sum_{i = 1}^n \frac{1}{M - 1} \sum_{m = 1}^M (p_m^* (x) - \hat{p}_m^* (x))^2}
#' 
#' where:
#' 
#' \deqn{p_m^* (x) = P(Y_i \leq m | X_i = x)}
#' 
#' The second formula is useful for simulation studies.
#'
#' @return The mean ranked probability score of the method.
#'
#' @export
mean_ranked_score <- function(y, predictions, use.true = FALSE){ # Taken from https://github.com/okasag/orf/blob/master/orf/R/evals.R
  ## Handling inputs and checks.
  if (!(use.true %in% c(TRUE, FALSE))) stop("'use.true' must be logical.", call. = FALSE)
  if (use.true == TRUE & is.null(dim(y)) | use.true == FALSE & !is.null(dim(y))) stop("Combination of 'use.true' and 'y' is non-sensical.")
  n <- dim(predictions)[1]
  n.classes <- dim(predictions)[2]

  ## Generating matrix for comparison. Either a matrix of indicator variables, or the true probability matrix.
  if (use.true == FALSE) {
    indicator_mat <- matrix(0, nrow = n, ncol = n.classes)
    for (i in seq_len(n)) {
      indicator_mat[i, y[i]] <- 1
    }
  } else if (use.true == TRUE) {
    indicator_mat <- y
  }
  
  ## RPS.
  # Pre-allocating memory.
  rps <- numeric(n)
  cum <- numeric(n)
  
  # Computing the inner sum of squares. At the i-th iteration, we sum up to the m-th class. 
  for (m in seq_len(n.classes)) {
    cum <- cum + (rowSums(matrix(predictions[, 1:m], ncol = m)) - rowSums(matrix(indicator_mat[, 1:m], ncol = m)))^2
  }
  
  # Computing RPS for each observation and averaging.
  rps <- (1 / (n.classes - 1)) * cum
  mrps <- mean(rps)

  ## Output.
  return(mrps)
}
