#' Marginal Effects for Modified Ordered Random Forests
#'
#' Nonparametric estimation of marginal effects using a \code{morf} object.
#'
#' @param object \code{morf} object.
#' @param eval Evaluation point for marginal effects. Either \code{"mean"}, \code{"atmean"} or \code{"atmedian"}.
#' @param bandwitdh How many standard deviations \code{x_up} and \code{x_down} differ from \code{x}.
#' @param data Data set of class \code{data.frame} to estimate marginal effects. It must contain at least the same covariates used to train the forests. If \code{NULL}, marginal effects are estimated on \code{object$full_data}.
#' @param inference Whether to conduct weight-based inference. It considerably slows down the program. Not advisable if \code{eval = "mean"}.
#' 
#' @return 
#' Object of class \code{morf.marginal} with elements:
#'   \item{\code{evaluation.type}}{Where the marginal effects are evaluated.}
#'   \item{\code{bandwitdh}}{The bandwitdh parameter.}
#'   \item{\code{n.classes}}{Number of classes.}
#'   \item{\code{n.samples}}{Number of samples.}
#'   \item{\code{n.trees}}{Number of trees in each forest.}
#'   \item{\code{marginal.effects}}{Matrix of marginal effects.}
#'   \item{\code{standard.errors}}{Matrix of standard errors, one for each marginal effect.}
#'   \item{\code{p.values}}{Matrix of p-values testing the null hypotheses that the marginal effect is zero.}
#'   \item{\code{ci.upper}}{Matrix of upper intervals.}
#'   \item{\code{ci.lower}}{Matrix of lower intervals.}
#' 
#' @details
#' The marginal effect of the j-th covariate on \eqn{p_m \left( \cdot \right)} is defined according to the continuous or 
#' discrete nature of the covariate: 
#'  
#' \deqn{p_{m, j}^\prime \left( x \right) := \frac{\partial P \left( Y_i = m \, | \, X_i = x \right)}{\partial x_{j}}}
#' 
#' \deqn{p_{m, j}^\prime \left( x \right) := P \left( Y_i = m \, | \, X_i = \lceil x_j \rceil \right) - P \left( Y_i = m \, | \, X_i = \lfloor x_j \rfloor \right)}
#'     
#' where \eqn{x_j} is the j-th element of the vector \eqn{x} and \eqn{\lceil x_j \rceil} and \eqn{\lfloor x_j \rfloor}
#' correspond to \eqn{x} with its j-th element rounded up and down to the closest integer.\cr  
#'                      
#' If \code{eval = "mean"}, then the mean marginal effects are estimated:
#' 
#' \deqn{MME := \frac{1}{n} \sum_{i = 1}^n p_{m, j}^\prime \left( X_i \right)}
#' 
#' If \code{eval = "atmean"}, then the marginal effects evaluated at the mean covariates are estimated:
#' 
#' \deqn{MEM := p_{m, j}^\prime ( \bar{X}_i )}  
#' 
#' \code{eval = "atmedian"} is similar to \code{eval = "atmean"}, where the median of the covariates is used as 
#' evaluation point rather than the mean.\cr
#' 
#' The program assumes that covariates with more than ten unique values are continuous. Otherwise, covariates are assumed
#' to be categorical or binary.
#'
#' @importFrom stats median sd pnorm
#'
#' @seealso \code{\link{morf}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item S Athey, J Tibshirani, S Wager (2019). Generalized random forests. The Annals of Statistics. \doi{10.1214/18-AOS1709}.
#'   \item Lechner, M., & Okasa, G. (2019). Random forest estimation of the ordered choice model. arXiv preprint arXiv:1907.02436. \doi{10.48550/arXiv.1907.02436}.
#'   \item Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. \doi{10.18637/jss.v077.i01}.
#' }
#'
#' @export
marginal_effects <- function(object, data = NULL, eval = "atmean", bandwitdh = 0.001, inference = FALSE) { # Inspired by https://github.com/okasag/orf/blob/master/orf/R/margins.R
  ## Handling inputs and checks.
  if (!inherits(object, "morf")) stop("Invalid class of input object.", call. = FALSE) 
  if (inference & !object$honesty) stop("Invalid inference if forests are not honest. Please feed in a morf object estimated with honesty = TRUE.", call. = FALSE)
  n_honest <- dim(object$honest_data)[1]
  y.classes <- object$classes
  n.classes <- object$n.classes
  
  ## Data.
  # Handle and check.
  if (is.null(data)) data <- object$full_data
  if (sum(!(object$forest.1$covariate.names %in% colnames(data))) > 0) stop("One or more covariates not found in 'data'.", call. = FALSE)
  if (length(colnames(data)) != length(object$forest.1$covariate.names) || 
      any(colnames(data) != object[[1]]$covariate.names)) data <- data[, object$forest.1$covariate.names, drop = FALSE]
  
  X <- data
  independent.variable.names <- colnames(X)
  if (length(independent.variable.names) < 1) stop("No covariates found.", call. = FALSE)

  # Save the covariates' types.
  X_unique <- apply(data, 2, function(x) length(unique(x)))
  X_continuous <- which(X_unique > 10) 
  X_dummy <- which(X_unique == 2) 
  X_categorical <- which(X_unique > 2 & X_unique <= 10)
  if (any(X_unique == 1) | any(X_unique == 0)) stop("Some of the covariates are constant. This makes no sense for evaluating marginal effects.", call. = FALSE)
  
  ## Evaluation points.
  if (!(eval %in% c("mean", "atmean", "atmedian"))) stop("Invalid value for 'eval'.", call. = FALSE)
  if (eval == "atmean") {
    evaluation_points <- as.data.frame(t(colMeans(X)))
  } else if (eval == "atmedian") {
    evaluation_points <- as.data.frame(t(apply(X, 2, median)))
  } else if (eval == "mean") {
    evaluation_points <- X 
  }
  
  ## Shifting each prediction point up and down.
  standard_deviations <- apply(X, 2, sd)
  X_up <- evaluation_points + bandwitdh * standard_deviations
  X_down <- evaluation_points - bandwitdh * standard_deviations
  
  # Enforce X_up and X_down in the support of X.
  n_rows <- nrow(evaluation_points)
  X_max <- matrix(rep(apply(X, 2, max), times = 1, each = n_rows), nrow = n_rows)
  X_min <- matrix(rep(apply(X, 2, min), times = 1, each = n_rows), nrow = n_rows)
  X_up <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
  X_up <- (X_up > X_min) * X_up + (X_up <= X_min) * (X_min + bandwitdh * standard_deviations)
  X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
  X_down <- (X_down < X_max) * X_down + (X_down >= X_max) * (X_max - bandwitdh * standard_deviations)
  
  # Check whether X_up and X_down are equal.
  while (any(X_up == X_down)) {
    X_up <- (X_up > X_down) * X_up + (X_up == X_down) * (X_up + (bandwitdh + 0.1) * standard_deviations)
    X_down <- (X_up > X_down) * X_down + (X_up == X_down) * (X_down - (bandwitdh + 0.1) * standard_deviations)

    X_up <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
    X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
  }
  
  ## Generating data for estimation. Each data set shifts the j-th covariate and leaves the others untouched.
  X_up_data <- list() 
  X_down_data <- list()
  
  for (j in seq_len(ncol(X))) {
    shifted_var_up <- X_up[, j, drop = FALSE] 
    shifted_var_down <- X_down[, j, drop = FALSE] 
    original_covariates <- evaluation_points[, -j, drop = FALSE] 
    
    X_up_data[[j]] <- data.frame(shifted_var_up, original_covariates)
    colnames(X_up_data[[j]])[1] <- names(X_up)[j]
    X_down_data[[j]] <- data.frame(shifted_var_down, original_covariates)
    colnames(X_down_data[[j]])[1] <- names(X_down)[j]
  }
  
  ## Correcting for discrete covariates. The shifted covariate is always in the first column.
  for (j in X_categorical) {
    X_up_data[[j]][1] <- ceiling(X_up_data[[j]][1])
    X_down_data[[j]][1] <- ifelse(ceiling(X_down_data[[j]][1]) == ceiling(X_up_data[[j]][1]), floor(X_down_data[[j]][1]), ceiling(X_down_data[[j]][1]))[[1]]
  }
  
  for (j in X_dummy) {
    X_up_data[[j]][1] <- max(X[, j])
    X_down_data[[j]][1] <- min(X[, j])
  }
  
  ## Difference in conditional class probabilities. 
  if (inference) {
    # Data for honest estimation. 
    honest_sample <- object$honest_data
    honest_outcomes <- list()
    counter <- 1
    
    for (m in y.classes) {
      honest_outcomes[[counter]] <- data.frame("y_m_honest" = ifelse(object$honest_data$y_honest <= m, 1, 0), "y_m_1_honest" = ifelse(object$honest_data$y_honest <= m -1, 1, 0))
      counter <- counter + 1
    }

    # Storing forests in a separate list. Forests are always the first n.classes elements of a morf object.
    forests <- list()
    
    for (m in seq_len(n.classes)) {
      forests[[m]] <- object[[m]]
    }
    
    # Extracting weights. List of lists: outer list concerns forests, inner lists concern shifted data.
    weights_up <- lapply(forests, function(x) {lapply(X_up_data, function(y) predict_forest_weights(x, honest_sample, test_sample = y))}) 
    weights_down <- lapply(forests, function(x) {lapply(X_down_data, function(y) predict_forest_weights(x, honest_sample, test_sample = y))}) 
    
    # Using weights for prediction. The j-th iteration uses data set with the j-th covariate shifted. Notice the normalization step.
    numerators <- list()
    
    for (j in seq_len(length(X_up_data))) {
      predictions_up <- mapply(function(x, y) {x[[j]] %*% (y$y_m_honest - y$y_m_1_honest)}, weights_up, honest_outcomes)
      predictions_down <- mapply(function(x, y) {x[[j]] %*% (y$y_m_honest - y$y_m_1_honest)}, weights_down, honest_outcomes)
      
      predictions_up <- matrix(predictions_up / rowSums(matrix(predictions_up, ncol = n.classes)), ncol = n.classes)
      predictions_down <- matrix(predictions_down / rowSums(matrix(predictions_down, ncol = n.classes)), ncol = n.classes)
      
      numerators[[j]] <- predictions_up - predictions_down
    }
  } else { 
    # Predictions shifting the j-th covariate. Normalization step done within predict.morf.
    predictions_up <- lapply(X_up_data, function(x) {predict(object, x)$predictions})
    predictions_down <- lapply(X_down_data, function(x) {predict(object, x)$predictions})
    numerators <- mapply(function(x, y) {x - y}, predictions_up, predictions_down, SIMPLIFY = FALSE)
  }
  
  ## Approximating the marginal change in the covariates. Enforce this to one for discrete covariates.
  denominators <- 2 * bandwitdh * standard_deviations
  for (i in (union(X_categorical, X_dummy))) {
    denominators[[i]] <- 1
  }
  
  ## Marginal effects. First, estimate them for each prediction point (we have only one if atmean/atmedian). Then, 
  ## average for each class with colMeans (this does not affect results if atmean/atmedian).
  marginal_effects <- mapply(function(x, y) {x / y}, numerators, denominators, SIMPLIFY = FALSE)
  marginal_effects <- matrix(unlist(lapply(marginal_effects, function(x) {colMeans(x)}), use.names = FALSE), ncol = n.classes, byrow = TRUE)
  colnames(marginal_effects) <- paste("P'(Y=", y.classes, ")", sep = "")
  rownames(marginal_effects) <- colnames(X)
  
  ## Variance of marginal effects.
  if (inference) {
    # Pre-allocating memory.
    variances <- matrix(NA, ncol = n.classes, nrow = length(denominators))
    colnames(variances) <- paste("P'(Y=", y.classes, ")", sep = "")
    rownames(variances) <- colnames(X)
    
    # Constants.
    denominators_squared <- 1 / denominators^2
    sample_correction <- n_honest / (n_honest - 1)
    
    # Building the variance. Notice that colMeans picks the average weight of each honest unit if mean and does 
    # nothing if atmean/atmedian.
    for (j in seq_len(length(X_up_data))) {
      weights_difference <- mapply(function(x, y) {colMeans(x[[j]] - y[[j]])}, weights_up, weights_down, SIMPLIFY = FALSE)
      products <- mapply(function(x, y) {x * (y$y_m_honest - y$y_m_1_honest)}, weights_difference, honest_outcomes, SIMPLIFY = FALSE)
      sum_squares <- lapply(products, function(x) {sum((x - mean(x))^2)})
      variances[j, ] <- unlist(lapply(sum_squares, function(x) {sample_correction * denominators_squared[j] * x}), use.names = FALSE)
    }
    
    standard_errors <- sqrt(variances)
    t_values <- marginal_effects / standard_errors
    t_values[is.infinite(t_values)] <- 0
    p_values <- 2 * pnorm(-abs(t_values))
    
    ci_upper <- marginal_effects + 1.96 * standard_errors 
    ci_lower <- marginal_effects - 1.96 * standard_errors
  }
  
  ## Handling output.
  results <- list()
  
  results$evaluation.type <- eval
  results$bandwitdh <- bandwitdh
  results$n.classes <- n.classes
  results$n.samples <- nrow(X)
  results$n.trees <- object$n.trees
  results$honesty <- object$honesty
  results$marginal.effects <- marginal_effects
  results$standard.errors <- if (inference) standard_errors else list()
  results$p.values <- if (inference) p_values else list()
  results$ci.upper <- if (inference) ci_upper else list()
  results$ci.lower <- if (inference) ci_lower else list()
  
  class(results) <- "morf.marginal"
  
  ## Output.
  return(results)
}