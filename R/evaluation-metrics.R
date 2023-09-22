#' Accuracy Measures for Ordered Probability Predictions
#'
#' Accuracy measures for evaluating ordered probability predictions.
#'
#' @param y Either the observed outcome vector or a matrix of true probabilities.
#' @param predictions Predictions.
#' @param use.true If \code{TRUE}, then the program treats \code{y} as a matrix of true probabilities.
#' 
#' @return The MSE, the MAE, the RPS, or the CE of the method.
#' 
#' @examples 
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#' 
#' data <- generate_ordered_data(1000)
#' sample <- data$sample
#' Y <- sample$Y
#' X <- sample[, -1]
#' 
#' ## Training-test split.
#' train_idx <- sample(seq_len(length(Y)), floor(length(Y) * 0.5))
#' 
#' Y_tr <- Y[train_idx]
#' X_tr <- X[train_idx, ]
#' 
#' Y_test <- Y[-train_idx]
#' X_test <- X[-train_idx, ]
#' 
#' ## Fit ocf on training sample.
#' forests <- ocf(Y_tr, X_tr)
#' 
#' ## Accuracy measures on test sample.
#' predictions <- predict(forests, X_test)
#' 
#' mean_squared_error(Y_test, predictions$probabilities)
#' mean_ranked_score(Y_test, predictions$probabilities)
#' classification_error(Y_test, predictions$classification)}
#' 
#' @md
#' @details 
#' ## MSE, MAE, and RPS
#' When calling one of \code{\link{mean_squared_error}}, \code{\link{mean_absolute_error}}, or \code{\link{mean_ranked_score}}, 
#' \code{predictions} must be a matrix of predicted  class probabilities, with as many rows as observations in \code{y} and as 
#' many columns as classes of \code{y}.\cr
#' 
#' If \code{use.true == FALSE}, the mean squared error (MSE), the mean absolute error (MAE), and the mean ranked probability score 
#' (RPS) are computed as follows:
#' 
#' \deqn{MSE = \frac{1}{n} \sum_{i = 1}^n \sum_{m = 1}^M (1 (Y_i = m) - \hat{p}_m (x))^2}
#' 
#' \deqn{MAE = \frac{1}{n} \sum_{i = 1}^n \sum_{m = 1}^M |1 (Y_i = m) - \hat{p}_m (x)|}
#' 
#' \deqn{RPS = \frac{1}{n} \sum_{i = 1}^n \frac{1}{M - 1} \sum_{m = 1}^M (1 (Y_i \leq m) - \hat{p}_m^* (x))^2}
#' 
#' If \code{use.true == TRUE}, the MSE, the MAE, and the RPS are computed as follows (useful for simulation studies):
#' 
#' \deqn{MSE = \frac{1}{n} \sum_{i = 1}^n \sum_{m = 1}^M (p_m (x) - \hat{p}_m (x))^2}
#' 
#' \deqn{MSE = \frac{1}{n} \sum_{i = 1}^n \sum_{m = 1}^M |p_m (x) - \hat{p}_m (x)|}
#' 
#' \deqn{RPS = \frac{1}{n} \sum_{i = 1}^n \frac{1}{M - 1} \sum_{m = 1}^M (p_m^* (x) - \hat{p}_m^* (x))^2}
#' 
#' where:
#' 
#' \deqn{p_m (x) = P(Y_i = m | X_i = x)}
#' 
#' \deqn{p_m^* (x) = P(Y_i \leq m | X_i = x)}
#' 
#' ## Classification error
#' When calling \code{\link{classification_error}}, \code{predictions} must be a vector of predicted class labels.\cr
#'  
#' Classification error (CE) is computed as follows:
#' 
#' \deqn{CE = \frac{1}{n} \sum_{i = 1}^n 1 (Y_i \neq \hat{Y}_i)}
#' 
#' where Y_i are the observed class labels. 
#' 
#' @author Riccardo Di Francesco
#' 
#' @seealso \code{\link{mean_ranked_score}}
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


#' Mean Absolute Error
#'
#' @rdname mean_squared_error
#'
#' @export
mean_absolute_error <- function(y, predictions, use.true = FALSE) {
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
  
  ## MAE.
  mae <- mean(rowSums(apply((indicator_mat - predictions), 2, function(x) {abs(x)})))
  
  ## Output.
  return(mae)
}


#' Mean Ranked Probability Score
#'
#' @rdname mean_squared_error
#'
#' @export
mean_ranked_score <- function(y, predictions, use.true = FALSE) { # Taken from https://github.com/okasag/orf/blob/master/orf/R/evals.R
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


#' Classification Error
#'
#' @rdname mean_squared_error
#'
#' @export
classification_error <- function(y, predictions) { 
  ## Classification error.
  ce <- mean(y != predictions)
  
  ## Output.
  return(ce)
}
