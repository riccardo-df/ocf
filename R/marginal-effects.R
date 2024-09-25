#' Marginal Effects for Ordered Correlation Forest
#'
#' Nonparametric estimation of marginal effects using an \code{\link{ocf}} object.
#'
#' @param object An \code{\link{ocf}} object.
#' @param eval Evaluation point for marginal effects. Either \code{"mean"}, \code{"atmean"} or \code{"atmedian"}.
#' @param which_covariates Character vector storing the names of the covariates for which marginal effect estimation is desired. If empty (the default), marginal effects are estimated for all covariates.
#' @param bandwitdh How many standard deviations \code{x_up} and \code{x_down} differ from \code{x}.
#' @param data Data set of class \code{data.frame} to estimate marginal effects. It must contain at least the same covariates used to train the forests. If \code{NULL}, marginal effects are estimated on \code{object$full_data}.
#' @param inference Whether to extract weights and compute standard errors. The weights extraction considerably slows down the program.
#' 
#' @return 
#' Object of class \code{ocf.marginal}.
#' 
#' @examples 
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#' 
#' data <- generate_ordered_data(100)
#' sample <- data$sample
#' Y <- sample$Y
#' X <- sample[, -1]
#' 
#' ## Fit ocf.
#' forests <- ocf(Y, X)
#' 
#' ## Marginal effects at the mean.
#' me <- marginal_effects(forests, eval = "atmean")
#' 
#' print(me)
#' print(me, latex = TRUE)
#' 
#' ## Compute standard errors. This requires honest forests.
#' honest_forests <- ocf(Y, X, honesty = TRUE)
#' 
#' honest_me <- marginal_effects(honest_forests, eval = "atmean", inference = TRUE)
#' 
#' print(honest_me, latex = TRUE)}
#' 
#' @details
#' \code{\link{marginal_effects}} can estimate mean marginal effects, marginal effects at the mean, or marginal effects at the
#' median, according to the \code{eval} argument.\cr 
#' 
#' The routine assumes that covariates with more than ten unique values are continuous. Otherwise, covariates are assumed to 
#' be discrete.\cr  
#'
#' @importFrom stats median sd pnorm
#' 
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2023). Ordered Correlation Forest. arXiv preprint \href{https://arxiv.org/abs/2309.08755}{arXiv:2309.08755}.
#' }
#'
#' @seealso \code{\link{ocf}}
#'
#' @export
marginal_effects <- function(object, data = NULL, which_covariates = c(), 
                             eval = "atmean", bandwitdh = 0.1, inference = FALSE) { # Inspired by https://github.com/okasag/orf/blob/master/orf/R/margins.R
  ## 1.) Handling inputs and checks.
  if (!inherits(object, "ocf")) stop("Invalid 'object'.", call. = FALSE) 
  if (inference & !object$tuning.info$honesty) stop("Inference requires forests to be honest. Please feed in a ocf object estimated with 'honesty = TRUE'.", call. = FALSE)
  
  n_honest <- dim(object$honest_data)[1]
  y.classes <- sort(unique(object$full_data[, 1]))
  n.classes <- length(y.classes)

  ## 2.) Data.
  # 2a.) Handle and check.
  if (is.null(data)) data <- object$full_data
  if (sum(!(object$forests.info$forest.1$covariate.names %in% colnames(data))) > 0) stop("One or more covariates not found in 'data'.", call. = FALSE)
  if (length(colnames(data)) != length(object$forests.info$forest.1$covariate.names)) data <- data[, object$forests.info$forest.1$covariate.names, drop = FALSE]
  
  X <- data
  if (length(colnames(X)) < 1) stop("No covariates found. Maybe 'X' is missing colnames?", call. = FALSE)
  if (!is.null(which_covariates) & !(any(which_covariates %in% colnames(X)))) stop("One or more of 'which_covariates' has not been found in 'data'.", call. = FALSE)

  # 2b.) Save the covariates' types.
  X_unique <- apply(data, 2, function(x) length(unique(x)))
  X_continuous <- which(X_unique > 10) 
  X_dummy <- which(X_unique == 2) 
  X_categorical <- which(X_unique > 2 & X_unique <= 10)
  if (any(X_unique == 1) | any(X_unique == 0)) stop("Some of the covariates are constant. This makes no sense for evaluating marginal effects.", call. = FALSE)
  
  ## 3.) Evaluation points.
  # 3a.) Extract evaluation points.
  if (!(eval %in% c("mean", "atmean", "atmedian"))) stop("Invalid value for 'eval'.", call. = FALSE)
  if (eval == "atmean") {
    evaluation_points <- as.data.frame(t(colMeans(X)))
  } else if (eval == "atmedian") {
    evaluation_points <- as.data.frame(t(apply(X, 2, median)))
  } else if (eval == "mean") {
    evaluation_points <- X 
  }
  
  # 3b.) Shift each point a little bit up and down.
  standard_deviations <- apply(X, 2, sd)
  X_up <- evaluation_points + bandwitdh * standard_deviations
  X_down <- evaluation_points - bandwitdh * standard_deviations
  
  # 3c.) Enforce X_up and X_down in the support of X.
  n_rows <- nrow(evaluation_points)
  X_max <- matrix(rep(apply(X, 2, max), times = 1, each = n_rows), nrow = n_rows)
  X_min <- matrix(rep(apply(X, 2, min), times = 1, each = n_rows), nrow = n_rows)
  X_up <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
  X_up <- (X_up > X_min) * X_up + (X_up <= X_min) * (X_min + bandwitdh * standard_deviations)
  X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
  X_down <- (X_down < X_max) * X_down + (X_down >= X_max) * (X_max - bandwitdh * standard_deviations)
  
  # 3d.) Ensure that X_up and X_down are not equal.
  while (any(X_up == X_down)) {
    X_up <- (X_up > X_down) * X_up + (X_up == X_down) * (X_up + (bandwitdh + 0.1) * standard_deviations)
    X_down <- (X_up > X_down) * X_down + (X_up == X_down) * (X_down - (bandwitdh + 0.1) * standard_deviations)

    X_up <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
    X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
  }
  
  ## 4.) Generate data for estimation. 
  # 4a.) Lists of data sets. The j-th data set shifts the j-th covariate and leaves the others untouched.
  X_up_data <- list() 
  X_down_data <- list()
  
  for (j in seq_len(ncol(X))) {
    shifted_var_up <- X_up[, j, drop = FALSE] 
    shifted_var_down <- X_down[, j, drop = FALSE] 
    original_covariates <- evaluation_points[, -j, drop = FALSE] 
    
    X_up_data[[j]] <- data.frame(shifted_var_up, original_covariates)
    colnames(X_up_data[[j]])[1] <- names(X_up)[j]
    colnames(X_up_data[[j]])[-1] <- names(X_up)[-j]
    X_down_data[[j]] <- data.frame(shifted_var_down, original_covariates)
    colnames(X_down_data[[j]])[1] <- names(X_down)[j]
    colnames(X_down_data[[j]])[-1] <- names(X_down)[-j]
  }
  
  # 4b.) Correct for discrete covariates. The shifted covariate is always in the first column.
  for (j in X_categorical) {
    X_up_data[[j]][1] <- ceiling(X_up_data[[j]][1])
    X_down_data[[j]][1] <- ifelse(ceiling(X_down_data[[j]][1]) == ceiling(X_up_data[[j]][1]), floor(X_down_data[[j]][1]), ceiling(X_down_data[[j]][1]))[[1]]
  }
  
  for (j in X_dummy) {
    X_up_data[[j]][1] <- max(X[, j])
    X_down_data[[j]][1] <- min(X[, j])
  }
  
  # 4c.) Drop covariates not of interest.
  if (!is.null(which_covariates)) {
    idx_interest <- which(colnames(X) %in% which_covariates) 
    X_up_data <- sapply(idx_interest, function(x) {X_up_data[[x]]}, simplify = FALSE)
    X_down_data <- sapply(idx_interest, function(x) {X_down_data[[x]]}, simplify = FALSE)
  } 

  ## 5.) Compute numerator: difference in conditional class probabilities. 
  if (inference) { # 5a.) If inference, extract weights and use these to predict.
    # 5aa.) Generate binary outcomes with honest sample. 
    honest_sample <- object$honest_data
    honest_outcomes <- list()
    counter <- 1
    for (m in y.classes) {
      honest_outcomes[[counter]] <- data.frame("y_m_honest" = ifelse(object$honest_data$Y <= m, 1, 0), "y_m_1_honest" = ifelse(object$honest_data$Y <= m -1, 1, 0))
      counter <- counter + 1
    }

    # 5ab.) Store forests.
    forests <- object$forests.info

    # 5ac.) Extract weights. List of lists: outer list concerns forests, inner lists concern shifted data.
    weights_up <- lapply(forests, function(x) {lapply(X_up_data, function(y) predict_forest_weights(x, honest_sample, test_sample = y))}) 
    weights_down <- lapply(forests, function(x) {lapply(X_down_data, function(y) predict_forest_weights(x, honest_sample, test_sample = y))}) 
    
    # 5ad.) Use weights for prediction. The j-th iteration uses data set with the j-th covariate shifted. Notice the normalization step.
    numerators <- list()
    
    for (j in seq_len(length(X_up_data))) {
      predictions_up <- mapply(function(x, y) {x[[j]] %*% (y$y_m_honest - y$y_m_1_honest)}, weights_up, honest_outcomes)
      predictions_down <- mapply(function(x, y) {x[[j]] %*% (y$y_m_honest - y$y_m_1_honest)}, weights_down, honest_outcomes)
      
      predictions_up <- matrix(predictions_up / rowSums(matrix(predictions_up, ncol = n.classes)), ncol = n.classes)
      predictions_down <- matrix(predictions_down / rowSums(matrix(predictions_down, ncol = n.classes)), ncol = n.classes)
      
      numerators[[j]] <- predictions_up - predictions_down
    }
  } else { # 5b.) If not inference, use standard prediction strategy. 
    # 5ba.) Predict by shifting the j-th covariate. Normalization step done within predict.ocf.
    predictions_up <- lapply(X_up_data, function(x) {predict(object, x)$probabilities})
    predictions_down <- lapply(X_down_data, function(x) {predict(object, x)$probabilities})
    numerators <- mapply(function(x, y) {x - y}, predictions_up, predictions_down, SIMPLIFY = FALSE)
  }
  
  ## 6.) Denominator: approximate the marginal change in the covariates, and enforce this to equal one for discrete covariates. Drop covariates not of interest.
  denominators <- 2 * bandwitdh * standard_deviations
  for (i in (union(X_categorical, X_dummy))) {
    denominators[[i]] <- 1
  }
  if (!is.null(which_covariates)) denominators <- denominators[idx_interest]
  
  ## 7.) Marginal effects. First, estimate them for each prediction point (we have only one if atmean/atmedian). Then, 
  ##     average for each class with colMeans (this does not affect results if atmean/atmedian).
  marginal_effects <- mapply(function(x, y) {x / y}, numerators, denominators, SIMPLIFY = FALSE)
  marginal_effects <- matrix(unlist(lapply(marginal_effects, function(x) {colMeans(x)}), use.names = FALSE), ncol = n.classes, byrow = TRUE)
  colnames(marginal_effects) <- paste0("P'(Y=", y.classes, ")")
  if (!is.null(which_covariates)) rownames(marginal_effects) <- colnames(X)[idx_interest] else rownames(marginal_effects) <- colnames(X)
  
  ## 8.) Variance of marginal effects, if necessary.
  if (inference) {
    # 8a.) Pre-allocating memory.
    variances <- matrix(NA, ncol = n.classes, nrow = length(denominators))
    colnames(variances) <- paste0("P'(Y=", y.classes, ")")
    if (!is.null(which_covariates)) rownames(variances) <- colnames(X)[idx_interest] else rownames(variances) <- colnames(X)
    
    # 8b.) Constants.
    denominators_squared <- 1 / denominators^2
    sample_correction <- n_honest / (n_honest - 1)
    
    # 8c.) Building the variance. Notice that colMeans picks the average weight of each honest unit if mean and does 
    #      nothing if atmean/atmedian.
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
  
  results$marginal.effects <- marginal_effects
  results$standard.errors <- if (inference) standard_errors else list()
  results$p.values <- if (inference) p_values else list()
  results$ci.upper <- if (inference) ci_upper else list()
  results$ci.lower <- if (inference) ci_lower else list()
  results$evaluation.type <- eval
  results$bandwitdh <- bandwitdh
  results$n.classes <- n.classes
  results$n.samples <- nrow(X)
  results$n.trees <- object$tuning.info$n.trees
  results$honesty <- object$tuning.info$honesty
  results$honesty.fraction <- object$tuning.info$honesty.fraction

  class(results) <- "ocf.marginal"
  
  ## Output.
  return(results)
}
