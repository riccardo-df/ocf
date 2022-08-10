##' Marginal Effects for Modified Ordered Random Forests
##'
##' Estimates the marginal effects for a \code{morf} object.
##'
##' @param object \code{morf} object.
##' @param eval Evaluation point for marginal effects. Either \code{"mean"}, \code{"atmean"} or \code{"atmedian"}.
##' @param bandwitdh How many standard deviations \code{x_up} and \code{x_down} differ from \code{x}.
##' @param data Data set of class \code{data.frame} to estimate marginal effects. It must contain at least the same covariates used to train the forests. If \code{data} is \code{NULL}, marginal effects are estimated on the training sample.
##' 
##' @details
##' \code{marginal_effects} provides a non-parametric estimation of the marginal effects for a modified ordered random
##' forest.\cr
##' 
##' If the k-th covariate is continuous, its marginal effect relative to the m-th class is defined as: 
##' 
##' \deqn{ME_{i, k}^m ( x ) = \frac{\partial P( Y_i = m \, | \, X_{i,k} = x_k, X_{i,-k} = x_{-k})}{\partial x_k}}
##' 
##' Otherwise, if k-th covariate is discrete:
##' 
##' \deqn{ME_{i, k}^m ( x ) = P( Y_i = m \, | \, X_{i,k} = \lceil x_k \rceil, X_{i,-k} = x_{-k}) - 
##'                      P( Y_i = m \, | \, X_{i,k} = \lfloor x_k \rfloor, X_{i,-k} = x_{-k})}
##' 
##' The program assumes that covariates with more than 10 unique values are continuous. Otherwise, covariates are assumed
##' to be categorical or binary.
##' 
##' @return 
##' Object of class \code{morf.marginal} with elements:
##'   \item{\code{evaluation.type}}{Where the marginal effects are evaluated.}
##'   \item{\code{bandwitdh}}{The bandwitdh parameter.}
##'   \item{\code{n.classes}}{Number of classes.}
##'   \item{\code{n.samples}}{Number of samples.}
##'   \item{\code{n.trees}}{Number of trees in each forest.}
##'   \item{\code{marginal.effects}}{Matrix of marginal effects.}
##'
##' @importFrom stats median sd
##'
##' @seealso \code{\link{morf}}.
##' 
##' @author Riccardo Di Francesco
##'
##' @export
marginal_effects <- function(object, data = NULL, eval = "mean", bandwitdh = 0.1) { # Inspired by https://github.com/okasag/orf/blob/master/orf/R/margins.R
  ## Handling inputs and checks.
  if (!inherits(object, "morf")) stop("Invalid class of input object.", call. = FALSE) 
  n.classes <- object$n.classes
  
  # Data.
  if (is.null(data)) data <- object$data
  if (sum(!(object$forest.1$covariate.names %in% colnames(data))) > 0) stop("One or more covariates not found in 'data'.", call. = FALSE)
  if (length(colnames(data)) != length(object$forest.1$covariate.names) || 
      any(colnames(data) != object$forest.1$covariate.names)) data <- data[, object$forest.1$covariate.names, drop = FALSE]
  
  X <- data
  
  independent.variable.names <- colnames(X)
  if (length(independent.variable.names) < 1) stop("No covariates found.", call. = FALSE)

  # Covariates type.
  X_unique <- apply(data, 2, function(x) length(unique(x)))
  
  X_continuous <- which(X_unique > 10) 
  X_dummy <- which(X_unique == 2) 
  X_categorical <- which(X_unique > 2 & X_unique <= 10)

  if (any(X_unique == 1) | any(X_unique == 0)) stop("Some of the covariates are constant. This makes no sense for evaluating marginal effects.", call. = FALSE)
  
  ## Computing marginal effects.
  # Evaluation points.
  if (!(eval %in% c("mean", "atmean", "atmedian"))) stop("Invalid value for 'eval'.", call. = FALSE)
  
  if (eval == "atmean") {
    evaluation_points <- as.data.frame(t(colMeans(X)))
  } else if (eval == "atmedian") {
    evaluation_points <- as.data.frame(t(apply(X, 2, median)))
  } else if (eval == "mean") {
    evaluation_points <- X 
  }
  
  # Generating x_up and x_down.
  standard_deviations <- apply(X, 2, sd)

  X_up <- evaluation_points + bandwitdh * standard_deviations
  X_down <- evaluation_points - bandwitdh * standard_deviations
  
  # Checking whether X_up and X_down are in the support of X.
  n_rows <- nrow(evaluation_points)

  X_max <- matrix(rep(apply(X, 2, max), times = 1, each = n_rows), nrow = n_rows)
  X_min <- matrix(rep(apply(X, 2, min), times = 1, each = n_rows), nrow = n_rows)

  X_up <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
  X_up <- (X_up > X_min) * X_up + (X_up <= X_min) * (X_min + bandwitdh * standard_deviations)

  X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
  X_down <- (X_down < X_max) * X_down + (X_down >= X_max) * (X_max - bandwitdh * standard_deviations)
  
  # Checking whether X_up and X_down are equal.
  while (any(X_up == X_down)) {
    X_up   <- (X_up > X_down) * X_up   + (X_up == X_down) * (X_up + (bandwitdh + 0.1) * standard_deviations)
    X_down <- (X_up > X_down) * X_down + (X_up == X_down) * (X_down - (bandwitdh + 0.1) * standard_deviations)

    X_up <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
    X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
  }
  
  # Generating one data frame for covariate. The k-th data frame stores the "shifted" k-th covariate and the "non-shifted" 
  # values for the other covariates.
  X_up_data <- list()
  X_down_data <- list()
  
  for (k in seq_len(ncol(X))) {
    shifted_var_up <- X_up[, k]
    shifted_var_down <- X_down[, k]
    
    original_covariates <- evaluation_points[, -k]
    
    X_up_data[[length(X_up_data) + 1]] <- data.frame(shifted_var_up, original_covariates)
    X_down_data[[length(X_down_data) + 1]] <- data.frame(shifted_var_down, original_covariates)
    
    colnames(X_up_data[[k]])[1] <- names(X_up)[k]
    colnames(X_down_data[[k]])[1] <- names(X_down)[k]
  }
  
  # Correcting according to covariate's type. "Shifted" covariate is always in the first column.
  for (k in X_categorical) {
    X_up_data[[k]][1] <- ceiling(X_up_data[[k]][1])
    X_down_data[[k]][1] <- ifelse(ceiling(X_down_data[[k]][1]) == ceiling(X_up_data[[k]][1]),
                                  floor(X_down_data[[k]][1]),
                                  ceiling(X_down_data[[k]][1]))[[1]]
  }
  
  for (k in X_dummy) {
    X_up_data[[k]][1] <- max(X[, k])
    X_down_data[[k]][1] <- min(X[, k])
  }
  
  # Difference in conditional class probabilities. The i-th element of the list stores the difference in predicted 
  # conditional class probabilities when the i-th covariate is "shifted".
  numerators <- list()

  for (i in seq_len(length(X_up_data))) {
    predictions_up <- colMeans(predict(object, X_up_data[[i]])$predictions)
    predictions_down <- colMeans(predict(object, X_down_data[[i]])$predictions)
    
    numerators[[i]] <- predictions_up - predictions_down
  }
  
  # Scaling factor. Due to estimation design, it equals the formula below. Notice we need to set the scaling factor
  # of discrete covariates to 1.
  denominators <- 2 * bandwitdh * standard_deviations
  
  for (i in (union(X_categorical, X_dummy))) {
    denominators[[i]] <- 1
  }
  
  # Ratio.
  marginal_effects <- matrix(NA, ncol = n.classes, nrow = length(denominators))
  colnames(marginal_effects) <- sapply(seq_len(n.classes), function(x) paste("class", x, sep = "."))
  rownames(marginal_effects) <- colnames(X)
  
  for(class in seq_len(ncol(X))) {
    marginal_effects[class, ] <- numerators[[class]] / denominators[class]
  }
  
  ## Handling output.
  results <- list()
  
  results$evaluation.type <- eval
  results$bandwitdh <- bandwitdh
  results$n.classes <- n.classes
  results$n.samples <- nrow(data)
  results$n.trees <- object$n.trees
  results$marginal.effects <- round(marginal_effects, 3)
  
  class(results) <- "morf.marginal"
  
  ## Output.
  return(results)
}


















