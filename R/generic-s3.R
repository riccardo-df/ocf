#' Prediction Method for morf Objects
#'
#' Prediction method for class \code{\link{morf}}.
#'
#' @param object An \code{\link{morf}} object.
#' @param data Data set of class \code{data.frame}. It must contain at least the same covariates used to train the forests. If \code{data} is \code{NULL}, then \code{object$full_data} is used.
#' @param type Type of prediction. Either \code{"response"} or \code{"terminalNodes"}. 
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Desired predictions.
#' 
#' @examples 
#' \donttest{
#' ## Load data from orf package.
#' set.seed(1986)
#' 
#' library(orf)
#' data(odata)
#' 
#' y <- as.numeric(odata[, 1])
#' X <- as.matrix(odata[, -1])
#' 
#' ## Training-test split.
#' train_idx <- sample(seq_len(length(y)), length(y)/2)
#' 
#' y_tr <- y[train_idx]
#' X_tr <- X[train_idx, ]
#' 
#' y_test <- y[-train_idx]
#' X_test <- X[-train_idx, ]
#' 
#' ## Fit morf on training sample.
#' forests <- morf(y_tr, X_tr)
#' 
#' ## Predict on test sample.
#' predictions <- predict(forests, X_test)
#' head(predictions$probabilities)
#' predictions$classification
#' 
#' ## Get terminal nodes.
#' predictions <- predict(forests, X_test, type = "terminalNodes")
#' predictions$forest.1[1:10, 1:20] # Rows are observations, columns are forests.}
#'
#' @details 
#' If \code{type == "response"}, the routine returns the predicted conditional class probabilities and the predicted class 
#' labels. If forests are honest, the predicted probabilities are honest.\cr
#' 
#' If \code{type == "terminalNodes"}, the IDs of the terminal node in each tree for each observation in \code{data} are returned.\cr
#'   
#' @seealso \code{\link{morf}}, \code{\link{marginal_effects}}
#' 
#' @importFrom stats predict
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
predict.morf <- function(object, data = NULL, type = "response", ...) {
  ## 0.) Default for variables not needed.
  predict.all <- FALSE; n.trees <- object$tuning.info$n.trees; n.threads <- NULL; verbose <- TRUE; seed <- NULL
  
  ## 1.) Handling inputs and check. 
  if (is.null(data)) data <- object$full_data
  y.classes <- sort(unique(object$full_data[, 1]))
  n.classes <- length(y.classes)
  forests.info <- object$forests.info
  
  ## 2.) Calling predict.morf.forest.
  prediction_output <- lapply(forests.info, function (x) {predict(x, data, type)})
  
  ## 3.) Handling prediction output, according to prediction type.
  # 3a.) Handle response type.
  if (type == "response") {
    # 3aa.) If honest forests, generate data sets and call honest function. Else, call standard function.
    if (object$tuning.info$honesty) { 
      honest_outcomes <- list()
      counter <- 1
      for (m in y.classes) {
        honest_outcomes[[counter]] <- data.frame("y_m_honest" = ifelse(object$honest_data$y <= m, 1, 0), "y_m_1_honest" = ifelse(object$honest_data$y <= m -1, 1, 0))
        counter <- counter + 1
      }
      
      predictions <- mapply(function (x, y) {honest_predictions(x, object$honest_data, data, y$y_m_honest, y$y_m_1_honest)}, forests.info, honest_outcomes)
    } else {
      predictions <- matrix(unlist(lapply(prediction_output, function (x) {x$predictions}), use.names = FALSE), ncol = n.classes, byrow = FALSE)
    }
    
    # Normalization step.
    predictions <- matrix(predictions / rowSums(matrix(predictions, ncol = n.classes)), ncol = n.classes)
    colnames(predictions) <- paste("P(Y=", y.classes, ")", sep = "")
    
    # Classification
    classification <- apply(predictions, 1, which.max)
    
    # Output.
    out <- list("probabilities" = predictions,
                "classification" = classification,
                "honesty" = object$tuning.info$honesty)
  } else if (type == "terminalNodes") {
    node_ids <- lapply(prediction_output, function (x) {x$predictions})
    names(node_ids) <- paste("forest.", y.classes, sep = "")
    out <- node_ids
  }
  
  ## Output.
  return(out)
}


#' Prediction Method for morf.forest Objects (Internal Use)
#'
#' Prediction method for class \code{morf.forest}.
#'
#' @param object An \code{morf.forest} object.
#' @param data Data set of class \code{data.frame}. It must contain at least the same covariates used to train the forests.
#' @param type Type of prediction. Either \code{"response"} or \code{"terminalNodes"}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Object of class \code{morf.prediction} with elements:
#'   \item{\code{predictions}}{If \code{type = "response"}, predicted conditional class probabilities.
#'                             If \code{type = "terminalNodes"}, the IDs of the terminal node in each tree for each observation.}
#'   
#' @details 
#' If \code{type === "response"} (the default), the predicted conditional class probabilities are returned. If forests are 
#' honest, these predictions are honest.\cr
#' 
#' If \code{type == "terminalNodes"}, the IDs of the terminal node in each tree for each observation in \code{data} are returned.
#'   
#' @seealso \code{\link{morf}}, \code{\link{marginal_effects}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
predict.morf.forest <- function(object, data, type = "response", ...) {
  ## 0.) Default for variables not needed.
  predict.all <- FALSE; n.trees <- object$num.trees; n.threads <- NULL; verbose <- TRUE; seed <- NULL; inbag.counts <- NULL
  
  ## 1.) Handling inputs and checks.
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
  
  ## 2.) Data.
  # 2a.) Handling and check.
  x <- data
  if (sum(!(forest$covariate.names %in% colnames(x))) > 0) stop("One or more covariates not found in 'data'.", call. = FALSE)
  
  # 2b.) Subset to same columns as in training sample.
  if (length(colnames(x)) != length(forest$covariate.names) || 
      any(colnames(x) != forest$covariate.names)) x <- x[, forest$covariate.names, drop = FALSE]
  
  # 2c.) Recode characters into factors.
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    char.columns <- sapply(x, is.character)
    if (length(char.columns) > 0) {
      x[char.columns] <- lapply(x[char.columns], factor)
    }
  }
  
  # 2d.) Data type.
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
  
  # 2e.) Missing values.
  if (any(is.na(x))) {
    offending_columns <- colnames(x)[colSums(is.na(x)) > 0]
    stop("Missing values in columns: ", paste0(offending_columns, collapse = ", "), ".", call. = FALSE)
  }
  
  ## 3.) Defaults for variables not needed.
  treetype <- 3; splitrule <- 1; mtry <- 0; max.depth <- 0; min.node.size <- 0; importance <- 0;
  prediction.mode <- TRUE; oob.error <- FALSE; y <- matrix(c(0, 0)); alpha_balance <- 0.2
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
  
  ## 4.) Calling morf in prediction mode.
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
                    regularization.factor, use.regularization.factor, regularization.usedepth,
                    alpha_balance)
  if (length(result) == 0) stop("User interrupt or internal error.", call. = FALSE)
  
  ## 5.) Handling output.
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


#' Summary Method for morf Objects
#' 
#' Summarizes an \code{\link{morf}} object.
#' 
#' @param object An \code{\link{morf}} object.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{morf}}, \code{\link{marginal_effects}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
summary.morf <- function(object, ...) {
  cat("Call: \n")
  cat(deparse(object$tuning.info$call), "\n\n")
  
  cat("Data info: \n")
  cat("Full sample size:  ", dim(object$full_data)[1], "\n")
  cat("N. covariates:     ", dim(object$full_data)[2], "\n")
  cat("Classes:           ", sort(unique(object$full_data[, 1])), "\n\n")
  
  cat("Relative variable importance: \n")
  print(round(object$importance, 3)); cat("\n")
  
  cat("Tuning parameters: \n")
  cat("N. trees:          ", object$tuning.info$n.trees, "\n")
  cat("mtry:              ", object$tuning.info$mtry, "\n")
  cat("min.node.size      ", object$tuning.info$min.node.size, "\n")
  if (object$tuning.info$replace) cat("Subsampling scheme:     Bootstrap \n" ) else cat("Subsampling scheme: No replacement \n" )
  cat("Honesty:           ", object$tuning.info$honesty, "\n")
  if(object$tuning.info$honesty) cat("Honest fraction:   ", object$tuning.info$honesty.fraction)
}
 
 
#' Print Method for morf Objects
#'
#' Prints an \code{\link{morf}} object.
#'
#' @param x An \code{\link{morf}} object.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{morf}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
print.morf <- function(x, ...) {
  summary.morf(x, ...)
}


#' Print Method for morf.marginal Objects
#'
#' Prints an \code{morf.marginal} object.
#'
#' @param x An \code{morf.marginal} object.
#' @param latex If \code{TRUE}, prints LATEX code.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @details 
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}. If
#' standard errors have been estimated, they are printed in parenthesis below each point estimate.
#' 
#' 
#' @seealso \code{\link{morf}}, \code{\link{marginal_effects}}.
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
print.morf.marginal <- function(x, latex = FALSE, ...) {
  if (!(latex %in% c(TRUE, FALSE))) stop("Invalid value of 'latex'.", call. = FALSE)
  
  if (latex) {
    table_names <- rename_latex(rownames(x$marginal.effects))
    
    cat("\\begingroup
    \\setlength{\\tabcolsep}{8pt}
    \\renewcommand{\\arraystretch}{1.1}
    \\begin{table}[H]
        \\centering
        \\begin{adjustbox}{width = 0.75\\textwidth}
        \\begin{tabular}{@{\\extracolsep{5pt}}l", rep(" c", x$n.classes), "}
        \\\\[-1.8ex]\\hline
        \\hline \\\\[-1.8ex]
        &", paste0(" Class ", seq_len(x$n.classes)[-x$n.classes], " &"), paste0(" Class ", x$n.classes), " \\\\
        \\addlinespace[2pt]
        \\hline \\\\[-1.8ex] \n\n", sep = "")
    
    if (inherits(x$standard.errors, "list")) {
      for (i in seq_len(nrow(x$marginal.effects))) {
        cat("        \\texttt{", table_names[i], "} & ", paste0(round(x$marginal.effects[i, 1:(ncol(x$marginal.effects)-1)], 3), " & "), round(x$marginal.effects[i, ncol(x$marginal.effects)], 3), " \\\\ \n", sep = "")
      }
    } else {
      for (i in seq_len(nrow(x$marginal.effects))) {
        cat("        \\texttt{", table_names[i], "} & ", paste0(round(x$marginal.effects[i, 1:(ncol(x$marginal.effects)-1)], 3), " & "), round(x$marginal.effects[i, ncol(x$marginal.effects)], 3), " \\\\ \n", sep = "")
        cat("                     & ", paste0("(", round(x$standard.errors[i, 1:(ncol(x$standard.errors)-1)], 3), ")", " & "), paste0("(", round(x$standard.errors[i, ncol(x$standard.errors)], 3), ")"), " \\\\ \n", sep = "")
      }
    }
    
    cat("\n        \\addlinespace[3pt]
        \\\\[-1.8ex]\\hline
        \\hline \\\\[-1.8ex]
        \\end{tabular}
        \\end{adjustbox}
        \\caption{https://soundcloud.com/theplasticchairband}
        \\label{table:descriptive.stats}
    \\end{table}
\\endgroup")
  } else {
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
}
