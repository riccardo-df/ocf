##' Prediction Method for Morf Objects
##'
##' Predicts conditional class probabilities with a \code{\link{morf}} object.
##'
##' @param object \code{morf} object.
##' @param data Data set of class \code{data.frame}. It must contain at least the same covariates used to train the forests. If \code{data} is \code{NULL}, predictions on the sample used to fit \code{object} are provided.
##' @param predict.all Return individual predictions for each tree instead of aggregated predictions for all trees (returns a matrix \code{n.samples} by \code{n.trees}). 
##' @param n.trees Number of trees used for prediction. The first \code{n.trees} in each forest are used. Default uses all trees in the forests.
##' @param type Type of prediction. One of \code{"response"} or \code{"terminalNodes"}. 
##' @param seed Random seed. Default is \code{NULL}, which generates the seed from \code{R}. Set to \code{0} to ignore the \code{R} seed. 
##' @param n.threads Number of threads. Default is number of CPUs available.
##' @param verbose Verbose output on or off.
##' @param ... Further arguments passed to or from other methods.
##' 
##' @return Object of class \code{morf.prediction} with elements:
##'   \item{\code{predictions}}{If \code{type = "response"}, predicted conditional class probabilities. 
##'                             If \code{type = "terminalNodes"}, the IDs of the terminal node in each tree for each observation.}
##'   \item{\code{n.trees}}{Number of trees.} 
##'   \item{\code{n.covariates}}{Number of covariates.}
##'   \item{\code{n.samples}}{Number of samples.}
##'   
##' @details 
##' For \code{type = 'response'} (the default), the predicted conditional class probabilities are returned.\cr
##' 
##' For \code{type = 'terminalNodes'}, the IDs of the terminal node in each tree for each observation in the given 
##' dataset are returned.\cr
##' 
##' If the forests are honest, then \code{predict.morf} yields honest predictions.
##' 
##' @references
##' \itemize{
##'   \item Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. J Stat Softw 77:1-17. \doi{10.18637/jss.v077.i01}.
##'   }
##'   
##' @seealso \code{\link{morf}}, , \code{\link{marginal_effects}}
##' 
##' @importFrom stats predict
##' 
##' @author Riccardo Di Francesco
##' 
##' @export
predict.morf <- function(object, data = NULL, type = "response",
                         predict.all = FALSE, n.trees = object$n.trees,
                         n.threads = NULL,
                         verbose = TRUE, seed = NULL, ...) {
  ## Handling inputs and check. It is enough to check the first forest.
  n.classes <- object$n.classes
  
  if (is.null(object$forest.1)) stop("No saved forest in morf object. Please set write.forest to TRUE when calling morf.", call. = FALSE)
  if (is.null(data)) data <- object$full_data
  
  ## Predicting for each class.
  if (type == "response") {
    predictions <- matrix(NA, nrow = dim(data)[1], ncol = n.classes) 
    
    # If the forest is honest, use honest prediction method from rcpp. This is required as we replace leaf estimates each time.
    if (object$honesty) { 
      for (m in seq_len(n.classes)) { # The morf.forest objects are always the first M elements of a morf object.
        temp <- predict(object[[m]], data, type, predict.all, n.trees, object$inbag.counts, n.threads, verbose, seed, ...) # Needed for output consistency.
        predictions[, m] <- honest_predictions(object[[m]], object$honest_data, data, m) # morf objects always store the honest sample, so we can use it for predictions.
      }
    } else { # Use default prediction method.
      for (m in seq_len(n.classes)) { 
        temp <- predict(object[[m]], data, type, predict.all, n.trees, object$inbag.counts, n.threads, verbose, seed, ...)
        predictions[, m] <- temp$predictions
      }
    }
    
    ## Normalization step.
    predictions <- matrix(apply(predictions, 1, function(x) (x)/(sum(x))), ncol = n.classes, byrow = T)
    colnames(predictions) <- sapply(seq_len(n.classes), function(x) paste("class", x, sep = "."))
    
    ## Handling output.
    out <- list("predictions" = predictions,
                "n.trees" = temp$n.trees,
                "n.covariates" = temp$n.covariates,
                "n.samples" = temp$num.samples)
  } else if (type == "terminalNodes") {
    node_ids <- list()
    for (m in seq_len(n.classes)) {
      temp <- predict(object[[m]], data, type, predict.all, n.trees, object$inbag.counts, n.threads, verbose, seed, ...)
      node_ids[[m]] <- temp$predictions
      names(node_ids)[m] <- paste("class", m, sep = ".")
      
      ## Handling output.
      out <- list("predictions" = node_ids,
                  "n.trees" = temp$n.trees,
                  "n.covariates" = temp$covariate.names,
                  "n.samples" = temp$n.samples)
    }
  }
  
  ## Output.
  class(out) <- "morf.prediction"
  return(out)
}


##' Prediction Method for Morf Forest Objects
##'
##' Predicts conditional class probabilities with a \code{morf.forest} object.
##'
##' @param object \code{morf} object.
##' @param data Data set of class \code{data.frame}. It must contain at least the same covariates used to train the forests.
##' @param predict.all Return individual predictions for each tree instead of aggregated predictions for all trees (returns a matrix \code{n.samples} by \code{n.trees}). 
##' @param n.trees Number of trees used for prediction. The first \code{n.trees} in the forest are used. Default uses all trees in the forest.
##' @param type Type of prediction. One of \code{"response"} or \code{"terminalNodes"}.
##' @param seed Random seed. Default is \code{NULL}, which generates the seed from \code{R}. Set to \code{0} to ignore the \code{R} seed. 
##' @param n.threads Number of threads. Default is number of CPUs available.
##' @param verbose Verbose output.
##' @param inbag.counts Number of times the observations are in-bag in the trees.
##' @param ... Further arguments passed to or from other methods.
##' 
##' @return Object of class \code{morf.prediction} with elements:
##'   \item{\code{predictions}}{If \code{type = "response"}, predicted conditional class probabilities. 
##'                             If \code{type = "terminalNodes"}, the IDs of the terminal node in each tree for each observation.}
##'   \item{\code{n.trees}}{Number of trees.}
##'   \item{\code{n.covariates}}{Number of covariates.}
##'   \item{\code{n.samples}}{Number of samples.}
##'   
##' @details 
##' For \code{type = 'response'} (the default), the predicted conditional class probabilities are returned. \cr
##' 
##' For \code{type = 'terminalNodes'}, the IDs of the terminal node in each tree for each observation in the given 
##' dataset are returned. 
##'   
##' @references
##' \itemize{
##'   \item Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. J Stat Softw 77:1-17. \doi{10.18637/jss.v077.i01}.
##'   }
##'   
##' @seealso \code{\link{morf}}, \code{\link{marginal_effects}}
##' 
##' @author Riccardo Di Francesco
##' 
##' @export
predict.morf.forest <- function(object, data, type = "response",
                                predict.all = FALSE, n.trees = object$num.trees,
                                inbag.counts = NULL,
                                n.threads = NULL, verbose = TRUE, seed = NULL, ...) {
  ## Handling inputs and checks.
  if (!inherits(object, "morf.forest")) stop("Invalid class of input object.", call. = FALSE) 
  
  forest <- object
  
  if (is.null(forest$num.trees) || is.null(forest$child.nodeIDs) || is.null(forest$split.varIDs) ||
        is.null(forest$split.values) || is.null(forest$covariate.names) || is.null(forest$treetype)) stop("Invalid forest object.", call. = FALSE)
  
  # Prediction type.
  if (type == "response") {
    prediction.type <- 1
  } else if (type == "terminalNodes") {
    prediction.type <- 2
  } else {
    stop("Invalid value for 'type'. Use 'response' or 'terminalNodes'.", call. = FALSE)
  }
  
  # Covariates.
  x <- data
  if (sum(!(forest$covariate.names %in% colnames(x))) > 0) stop("One or more covariates not found in 'data'.", call. = FALSE)

  # Subsetting to same columns as in training sample if necessary.
  if (length(colnames(x)) != length(forest$covariate.names) || 
      any(colnames(x) != forest$covariate.names)) x <- x[, forest$covariate.names, drop = FALSE]

  # Recoding characters into factors.
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    char.columns <- sapply(x, is.character)
    if (length(char.columns) > 0) {
      x[char.columns] <- lapply(x[char.columns], factor)
    }
  }
  
  # Handling data set.
  if (is.list(x) && !is.data.frame(x)) x <- as.data.frame(x)

  # Converting to data matrix.
  if (!is.matrix(x) & !inherits(x, "Matrix")) x <- data.matrix(x)
  
  # Missing values.
  if (any(is.na(x))) {
    offending_columns <- colnames(x)[colSums(is.na(x)) > 0]
    stop("Missing values in columns: ", paste0(offending_columns, collapse = ", "), ".", call. = FALSE)
  }

  # Number of threads.
  if (is.null(n.threads)) {
    n.threads = 0
  } else if (!is.numeric(n.threads) || n.threads < 0) {
    stop("Invalid value for n.threads", call. = FALSE)
  }

  # Seed.
  if (is.null(seed)) seed <- stats::runif(1 , 0, .Machine$integer.max)

  # Defaults for variables not needed.
  treetype <- 3
  mtry <- 0
  importance <- 0
  min.node.size <- 0
  split.select.weights <- list(c(0, 0))
  use.split.select.weights <- FALSE
  always.split.variables <- c("0", "0")
  use.always.split.variables <- FALSE
  prediction.mode <- TRUE
  write.forest <- FALSE
  replace <- TRUE
  probability <- FALSE
  unordered.factor.variables <- c("0", "0")
  use.unordered.factor.variables <- FALSE
  save.memory <- FALSE
  splitrule <- 1
  alpha <- 0
  minprop <- 0
  case.weights <- c(0, 0)
  use.case.weights <- FALSE
  class.weights <- c(0, 0)
  keep.inbag <- FALSE
  sample.fraction <- 1
  holdout <- FALSE
  num.random.splits <- 1
  order.snps <- FALSE
  oob.error <- FALSE
  max.depth <- 0
  inbag <- list(c(0,0))
  use.inbag <- FALSE
  y <- matrix(c(0, 0))
  regularization.factor <- c(0, 0)
  use.regularization.factor <- FALSE
  regularization.usedepth <- FALSE
  snp.data <- as.matrix(0)
  gwa.mode <- FALSE
  
  # Use sparse matrix.
  if (inherits(x, "dgCMatrix")) {
    sparse.x <- x
    x <- matrix(c(0, 0))
    use.sparse.data <- TRUE
  } else {
    sparse.x <- Matrix::Matrix(matrix(c(0, 0)))
    use.sparse.data <- FALSE
    x <- data.matrix(x)
  }
  
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
