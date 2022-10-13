#' Prediction Method for Morf Objects
#'
#' Prediction method for class \code{\link{morf}}.
#'
#' @param object \code{morf} object.
#' @param data Data set of class \code{data.frame}. It must contain at least the same covariates used to train the forests. If \code{data} is \code{NULL}, then \code{object$full_sample} is used.
#' @param predict.all Return individual predictions for each tree instead of aggregated predictions for all trees (returns a matrix \code{n.samples} by \code{n.trees}). 
#' @param n.trees Number of trees used for prediction. The first \code{n.trees} in each forest are used. Default uses all trees in the forests.
#' @param type Type of prediction. One of \code{"response"} or \code{"terminalNodes"}. 
#' @param seed Random seed. Default is \code{NULL}, which generates the seed from \code{R}. Set to \code{0} to ignore the \code{R} seed. 
#' @param n.threads Number of threads. Default is number of CPUs available.
#' @param verbose Verbose output on or off.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return Object of class \code{morf.prediction} with elements:
#'   \item{\code{predictions}}{If \code{type = "response"}, predicted conditional class probabilities. If forests are honest, then these predictions are honest.
#'                             If \code{type = "terminalNodes"}, the IDs of the terminal node in each tree for each observation.}
#'   \item{\code{n.trees}}{Number of trees.} 
#'   \item{\code{n.covariates}}{Number of covariates.}
#'   \item{\code{n.samples}}{Number of samples.}
#'   
#' @details 
#' For \code{type = 'response'} (the default), the predicted conditional class probabilities are returned. If forests are 
#' honest, then these predictions are honest.\cr
#' 
#' For \code{type = 'terminalNodes'}, the IDs of the terminal node in each tree for each observation in the given 
#' dataset are returned.\cr
#' 
#' @references
#' \itemize{
#'   \item Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. J Stat Softw 77:1-17. \doi{10.18637/jss.v077.i01}.
#'   }
#'   
#' @seealso \code{\link{morf}}, \code{\link{marginal_effects}}
#' 
#' @importFrom stats predict
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
predict.morf <- function(object, data = NULL, type = "response",
                         predict.all = FALSE, n.trees = object$n.trees,
                         n.threads = NULL, verbose = TRUE, seed = NULL, ...) {
  ## Handling inputs and check. 
  if (is.null(data)) data <- object$full_data
  y.classes <- object$classes
  n.classes <- object$n.classes
  
  # Storing forests in a separate list. Forests are always the first n.classes elements of a morf object.
  forests <- list()
  for (m in seq_len(n.classes)) {
    forests[[m]] <- object[[m]]
  }
  
  ## Calling predict.morf.forest.
  prediction_output <- lapply(forests, function (x) {predict(x, data, type, predict.all, n.trees, object$inbag.counts, 
                                                             n.threads, verbose, seed)})
  
  ## Handling prediction output, according to prediction type.
  if (type == "response") {
    if (object$honesty) { 
      honest_outcomes <- list()
      counter <- 1
      for (m in y.classes) {
        honest_outcomes[[counter]] <- data.frame("y_m_honest" = ifelse(object$honest_data$y_honest <= m, 1, 0), "y_m_1_honest" = ifelse(object$honest_data$y_honest <= m -1, 1, 0))
        counter <- counter + 1
      }
      predictions <- mapply(function (x, y) {honest_predictions(x, object$honest_data, data, y$y_m_honest, y$y_m_1_honest)}, forests, honest_outcomes)
    } else {
      predictions <- matrix(unlist(lapply(prediction_output, function (x) {x$predictions}), use.names = FALSE), ncol = n.classes, byrow = FALSE)
    }
    
    # Normalization step.
    predictions <- matrix(predictions / rowSums(matrix(predictions, ncol = n.classes)), ncol = n.classes)
    colnames(predictions) <- paste("P(Y=", y.classes, ")", sep = "")
    
    # Output.
    out <- list("predictions" = predictions,
                "n.trees" = prediction_output[[1]]$n.trees,
                "n.covariates" = prediction_output[[1]]$n.covariates,
                "n.samples" = prediction_output[[1]]$num.samples)
  } else if (type == "terminalNodes") {
    node_ids <- lapply(prediction_output, function (x) {x$predictions})
    names(node_ids) <- paste("P(Y=", y.classes, ")", sep = "")
    out <- list("predictions" = node_ids,
                "n.trees" = prediction_output[[1]]$n.trees,
                "n.covariates" = prediction_output[[1]]$covariate.names,
                "n.samples" = prediction_output[[1]]$n.samples)
  }
  
  ## Output.
  class(out) <- "morf.prediction"
  return(out)
}


#' Prediction Method for Morf Forest Objects
#'
#' Prediction method for class \code{morf.forest}.
#'
#' @param object \code{morf} object.
#' @param data Data set of class \code{data.frame}. It must contain at least the same covariates used to train the forests.
#' @param predict.all Return individual predictions for each tree instead of aggregated predictions for all trees (returns a matrix \code{n.samples} by \code{n.trees}). 
#' @param n.trees Number of trees used for prediction. The first \code{n.trees} in the forest are used. Default uses all trees in the forest.
#' @param type Type of prediction. One of \code{"response"} or \code{"terminalNodes"}.
#' @param seed Random seed. Default is \code{NULL}, which generates the seed from \code{R}. Set to \code{0} to ignore the \code{R} seed. 
#' @param n.threads Number of threads. Default is number of CPUs available.
#' @param verbose Verbose output.
#' @param inbag.counts Number of times the observations are in-bag in the trees.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return Object of class \code{morf.prediction} with elements:
#'   \item{\code{predictions}}{If \code{type = "response"}, predicted conditional class probabilities. If forests are honest, then these predictions are honest.
#'                             If \code{type = "terminalNodes"}, the IDs of the terminal node in each tree for each observation.}
#'   \item{\code{n.trees}}{Number of trees.} 
#'   \item{\code{n.covariates}}{Number of covariates.}
#'   \item{\code{n.samples}}{Number of samples.}
#'   
#' @details 
#' For \code{type = 'response'} (the default), the predicted conditional class probabilities are returned. If forests are 
#' honest, then these predictions are honest.\cr
#' 
#' For \code{type = 'terminalNodes'}, the IDs of the terminal node in each tree for each observation in the given 
#' dataset are returned.\cr
#'   
#' @references
#' \itemize{
#'   \item Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. J Stat Softw 77:1-17. \doi{10.18637/jss.v077.i01}.
#'   }
#'   
#' @seealso \code{\link{morf}}, \code{\link{marginal_effects}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
predict.morf.forest <- function(object, data, type = "response",
                                predict.all = FALSE, n.trees = object$num.trees,
                                inbag.counts = NULL,
                                n.threads = NULL, verbose = TRUE, seed = NULL, ...) {
  ## Handling inputs and checks.
  # Generic checks.
  forest <- object
  if (!inherits(forest, "morf.forest")) stop("Invalid class of input object.", call. = FALSE) 
  if (is.null(forest$num.trees) || is.null(forest$child.nodeIDs) || is.null(forest$split.varIDs) ||
      is.null(forest$split.values) || is.null(forest$covariate.names) || is.null(forest$treetype)) stop("Invalid forest object.", call. = FALSE)
  
  if (type == "response") {
    prediction.type <- 1
  } else if (type == "terminalNodes") {
    prediction.type <- 2
  } else {
    stop("Invalid value for 'type'.", call. = FALSE)
  }
  
  if (is.null(n.threads)) {
    n.threads = 0
  } else if (!is.numeric(n.threads) || n.threads < 0) {
    stop("Invalid value for n.threads", call. = FALSE)
  }
  
  if (is.null(seed)) seed <- stats::runif(1 , 0, .Machine$integer.max)
  
  ## Data.
  # Handling and check.
  x <- data
  if (sum(!(forest$covariate.names %in% colnames(x))) > 0) stop("One or more covariates not found in 'data'.", call. = FALSE)
  
  # Subset to same columns as in training sample.
  if (length(colnames(x)) != length(forest$covariate.names) || 
      any(colnames(x) != forest$covariate.names)) x <- x[, forest$covariate.names, drop = FALSE]
  
  # Recode characters into factors.
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    char.columns <- sapply(x, is.character)
    if (length(char.columns) > 0) {
      x[char.columns] <- lapply(x[char.columns], factor)
    }
  }
  
  # Data type.
  if (is.list(x) && !is.data.frame(x)) x <- as.data.frame(x)
  if (!is.matrix(x) & !inherits(x, "Matrix")) x <- data.matrix(x)
  
  if (inherits(x, "dgCMatrix")) {
    sparse.x <- x
    x <- matrix(c(0, 0))
    use.sparse.data <- TRUE
  } else {
    sparse.x <- Matrix::Matrix(matrix(c(0, 0)))
    use.sparse.data <- FALSE
    x <- data.matrix(x)
  }
  
  # Missing values.
  if (any(is.na(x))) {
    offending_columns <- colnames(x)[colSums(is.na(x)) > 0]
    stop("Missing values in columns: ", paste0(offending_columns, collapse = ", "), ".", call. = FALSE)
  }
  
  ## Defaults for variables not needed.
  treetype <- 3; splitrule <- 1; mtry <- 0; max.depth <- 0; min.node.size <- 0; importance <- 0;
  prediction.mode <- TRUE; oob.error <- FALSE; y <- matrix(c(0, 0))
  split.select.weights <- list(c(0, 0)); use.split.select.weights <- FALSE
  always.split.variables <- c("0", "0"); use.always.split.variables <- FALSE
  write.forest <- FALSE; replace <- TRUE; sample.fraction <- 1
  probability <- FALSE
  unordered.factor.variables <- c("0", "0"); use.unordered.factor.variables <- FALSE;   order.snps <- FALSE
  save.memory <- FALSE
  alpha <- 0; minprop <- 0; num.random.splits <- 1
  case.weights <- c(0, 0); use.case.weights <- FALSE; class.weights <- c(0, 0)
  keep.inbag <- FALSE; holdout <- FALSE; inbag <- list(c(0,0)); use.inbag <- FALSE
  regularization.factor <- c(0, 0); use.regularization.factor <- FALSE; regularization.usedepth <- FALSE
  snp.data <- as.matrix(0); gwa.mode <- FALSE
  
  ## Calling Morf in prediction mode.
  result <- morfCpp(treetype, x, y, forest$covariate.names, mtry,
                    n.trees, verbose, seed, n.threads, write.forest, importance,
                    min.node.size, split.select.weights, use.split.select.weights,
                    always.split.variables, use.always.split.variables,
                    prediction.mode, forest, snp.data, replace, probability,
                    unordered.factor.variables, use.unordered.factor.variables, save.memory, splitrule,
                    case.weights, use.case.weights, class.weights, 
                    predict.all, keep.inbag, sample.fraction, alpha, minprop, holdout, 
                    prediction.type, num.random.splits, sparse.x, use.sparse.data,
                    order.snps, oob.error, max.depth, inbag, use.inbag, 
                    regularization.factor, use.regularization.factor, regularization.usedepth)
  if (length(result) == 0) stop("User interrupt or internal error.", call. = FALSE)
  
  ## Handling output.
  result$num.samples <- nrow(x)
  result$treetype <- forest$treetype
  
  if (predict.all) {
    if (is.list(result$predictions)) {
      result$predictions <- do.call(rbind, result$predictions)
    } else {
      result$predictions <- array(result$predictions, dim = c(1, length(result$predictions)))
    }
  } else {
    if (is.list(result$predictions)) {
      result$predictions <- do.call(rbind, result$predictions)
    } 
  }
  
  if (type == "terminalNodes") {
    if (is.vector(result$predictions)) {
      result$predictions <- matrix(result$predictions, nrow = 1)
    }
  }
  
  return(result)
}


#' Plot Method for Morf Objects
#' 
#' Plots a \code{morf} object.
#' 
#' @param x \code{morf} object.
#' @param multiple.panels Whether to plot each class in a separate panel.
#' @param ... Further arguments passed to or from other methods.
#'
#' @import ggplot2
#' @importFrom utils stack
#' 
#' @seealso \code{\link{morf}}, \code{\link{marginal_effects}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
plot.morf <- function(x, multiple.panels = FALSE, ...) {
  ## Handling inputs.
  probabilities <- stack(as.data.frame(x$predictions))
  
  if (multiple.panels) {
    ggplot(data = probabilities, aes(x = values, fill = ind)) +
      geom_density(alpha = 0.4) +
      facet_wrap(~ind) + 
      xlab("Predicted probability") + ylab("Density") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    ggplot(data = probabilities, aes(x = values, fill = ind)) +
      geom_density(alpha = 0.4) +
      xlab("Predicted probability") + ylab("Density") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  }
}


#' Summary Method for Morf Objects
#' 
#' Summarizes a \code{morf} object.
#' 
#' @param object \code{morf} object.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{morf}}, \code{\link{marginal_effects}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
summary.morf <- function(object, ...) {
  cat("Call: \n")
  cat(deparse(x$call), "\n\n")
  
  cat("Classes: \n")
  cat(object$classes, "\n\n")
  
  cat("Variable importance: \n")
  print(round(object$overall.importance, 3)); cat("\n")
  
  cat("Forests info: \n")
  cat("Sample size:       ", object$n.samples, "\n")
  cat("N.trees:           ", object$n.trees, "\n")
  cat("mtry:              ", object$mtry, "\n")
  cat("min.node.size      ", object$min.node.size, "\n")
  if (object$replace) cat("Subsampling scheme:     Bootstrap \n" ) else cat("Subsampling scheme: No replacement \n" )
  cat("Honesty:           ", object$honesty, "\n")
  if(object$honesty) cat("Honest fraction:   ", object$honesty.fraction)
  cat("\n\n")
  
  cat("In-sample accuracy: \n")
  cat("MSE: ", round(object$mean.squared.error, 3), "\n")
  cat("RPS: ", round(object$mean.ranked.score, 3))
}
 
 
#' Print Method for Morf Objects
#'
#' Prints a \code{morf} object.
#'
#' @param x \code{morf} object.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{morf}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
print.morf <- function(x, ...) {
  cat("Call: \n")
  cat(deparse(x$call), "\n\n")
  cat("Number of classes:               ", x$n.classes, "\n")
  cat("Number of trees:                 ", x$n.trees, "\n")
  cat("Sample size:                     ", x$n.samples, "\n")
  cat("Number of covariates:            ", x$n.covariates, "\n")
  cat("Mtry:                            ", x$mtry, "\n")
  cat("Minimum node size:               ", x$min.node.size, "\n")
  cat("Honesty:                         ", x$honesty, "\n")
  if (x$honesty) cat("Fraction honesty:                ", x$honesty.fraction, "\n")
  cat("MSE:                             ", round(x$mean.squared.error, 3), "\n")
  cat("RPS:                             ", round(x$mean.ranked.score, 3))
}


#' Print Method for Morf Marginal Effects
#'
#' Prints a \code{morf.marginal} object.
#'
#' @param x \code{morf.marginal} object.
#' @param latex If \code{TRUE}, prints a latex code for a table displaying the marginal effects.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{morf}} and \code{\link{marginal_effects}}.
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
print.morf.marginal <- function(x, latex = FALSE, ...) {
  if (!(latex %in% c(TRUE, FALSE))) stop("Invalid value of 'latex'.", call. = FALSE)
  
  cat("Morf marginal effects results \n\n")
  cat("Evaluation:                      ", x$evaluation, "\n")
  cat("Bandwidth:                       ", x$bandwitdh, "\n")
  cat("Number of classes:               ", x$n.classes, "\n")
  cat("Number of trees:                 ", x$n.trees, "\n")
  cat("Sample size:                     ", x$n.samples, "\n")
  cat("Honest forests:                  ", x$honesty, "\n")
  
  cat("\n\n")
  
  cat("Marginal Effects: \n\n")
  
  print(round(x$marginal.effects, 3))
}
