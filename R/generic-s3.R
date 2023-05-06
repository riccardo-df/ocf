#' Prediction Method for ocf Objects
#'
#' Prediction method for class \code{\link{ocf}}.
#'
#' @param object An \code{\link{ocf}} object.
#' @param data Data set of class \code{data.frame}. It must contain at least the same covariates used to train the forests. If \code{data} is \code{NULL}, then \code{object$full_data} is used.
#' @param type Type of prediction. Either \code{"response"} or \code{"terminalNodes"}. 
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Desired predictions.
#' 
#' @examples 
#' ## Load data from orf package.
#' set.seed(1986)
#' 
#' library(orf)
#' data(odata)
#' odata <- odata[1:200, ] # Subset to reduce elapsed time.
#' 
#' y <- as.numeric(odata[, 1])
#' X <- as.matrix(odata[, -1])
#' 
#' ## Training-test split.
#' train_idx <- sample(seq_len(length(y)), floor(length(y) * 0.5))
#' 
#' y_tr <- y[train_idx]
#' X_tr <- X[train_idx, ]
#' 
#' y_test <- y[-train_idx]
#' X_test <- X[-train_idx, ]
#' 
#' ## Fit ocf on training sample.
#' forests <- ocf(y_tr, X_tr)
#' 
#' ## Predict on test sample.
#' predictions <- predict(forests, X_test)
#' head(predictions$probabilities)
#' predictions$classification
#' 
#' ## Get terminal nodes.
#' predictions <- predict(forests, X_test, type = "terminalNodes")
#' predictions$forest.1[1:10, 1:20] # Rows are observations, columns are forests.
#'
#' @details 
#' If \code{type == "response"}, the routine returns the predicted conditional class probabilities and the predicted class 
#' labels. If forests are honest, the predicted probabilities are honest.\cr
#' 
#' If \code{type == "terminalNodes"}, the IDs of the terminal node in each tree for each observation in \code{data} are returned.\cr
#'   
#' @seealso \code{\link{ocf}}, \code{\link{marginal_effects}}
#' 
#' @importFrom stats predict
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
predict.ocf <- function(object, data = NULL, type = "response", ...) {
  ## 0.) Default for variables not needed.
  predict.all <- FALSE; n.trees <- object$tuning.info$n.trees; n.threads <- NULL; verbose <- TRUE; seed <- NULL
  
  ## 1.) Handling inputs and check. 
  if (is.null(data)) data <- object$full_data
  y.classes <- sort(unique(object$full_data[, 1]))
  n.classes <- length(y.classes)
  forests.info <- object$forests.info
  
  ## 2.) Handling prediction output, according to prediction type.
  if (type == "response") { # 2a.) Under "response", generate binary outcomes and call honest function if honest forests, else call standard function. In any case, proceeds with normalization step, predict labels, and construct output.
    if (object$tuning.info$honesty) { 
      honest_outcomes <- list()
      counter <- 1
      for (m in y.classes) {
        honest_outcomes[[counter]] <- data.frame("y_m_honest" = ifelse(object$honest_data$y <= m, 1, 0), "y_m_1_honest" = ifelse(object$honest_data$y <= m -1, 1, 0))
        counter <- counter + 1
      }
      
      predictions <- mapply(function (x, y) {honest_predictions(x, object$honest_data, data, y$y_m_honest, y$y_m_1_honest)}, forests.info, honest_outcomes)
    } else {
      prediction_output <- lapply(forests.info, function (x) {predict(x, data, type)})
      predictions <- matrix(unlist(lapply(prediction_output, function (x) {x$predictions}), use.names = FALSE), ncol = n.classes, byrow = FALSE)
    }
    
    predictions <- matrix(predictions / rowSums(matrix(predictions, ncol = n.classes)), ncol = n.classes)
    colnames(predictions) <- paste("P(Y=", y.classes, ")", sep = "")
    classification <- apply(predictions, 1, which.max)
    
    out <- list("probabilities" = predictions,
                "classification" = classification,
                "honesty" = object$tuning.info$honesty)
  } else if (type == "terminalNodes") { # 2b.) Under "terminalNodes", extract leaf ids.
    prediction_output <- lapply(forests.info, function (x) {predict(x, data, type)})
    node_ids <- lapply(prediction_output, function (x) {x$predictions})
    names(node_ids) <- paste("forest.", y.classes, sep = "")
    out <- node_ids
  }
  
  ## Output.
  return(out)
}


#' Prediction Method for ocf.forest Objects
#'
#' Prediction method for class \code{ocf.forest}.
#'
#' @param object An \code{ocf.forest} object.
#' @param data Data set of class \code{data.frame}. It must contain at least the same covariates used to train the forests.
#' @param type Type of prediction. Either \code{"response"} or \code{"terminalNodes"}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Prediction results.
#' 
#' @keywords internal
#' 
#' @details 
#' If \code{type === "response"} (the default), the predicted conditional class probabilities are returned. If forests are 
#' honest, these predictions are honest.\cr
#' 
#' If \code{type == "terminalNodes"}, the IDs of the terminal node in each tree for each observation in \code{data} are returned.
#'   
#' @seealso \code{\link{ocf}}, \code{\link{marginal_effects}}
#' 
#' @author Riccardo Di Francesco
predict.ocf.forest <- function(object, data, type = "response", ...) {
  ## 0.) Default for variables not needed.
  predict.all <- FALSE; n.trees <- object$num.trees; n.threads <- 0; verbose <- TRUE; inbag.counts <- NULL; treetype <- 3; mtry <- 0;
  seed <- stats::runif(1 , 0, .Machine$integer.max); splitrule <- 1; max.depth <- 0; min.node.size <- 0; importance <- 0; 
  prediction.mode <- TRUE; oob.error <- FALSE; y <- matrix(c(0, 0)); alpha_balance <- 0; split.select.weights <- list(c(0, 0)); 
  use.split.select.weights <- FALSE; always.split.variables <- c("0", "0"); use.always.split.variables <- FALSE 
  write.forest <- FALSE; replace <- TRUE; sample.fraction <- 1; probability <- FALSE; unordered.factor.variables <- c("0", "0"); 
  use.unordered.factor.variables <- FALSE;   order.snps <- FALSE; save.memory <- FALSE; alpha <- 0; minprop <- 0;
  num.random.splits <- 1; case.weights <- c(0, 0); use.case.weights <- FALSE; class.weights <- c(0, 0); keep.inbag <- FALSE; 
  holdout <- FALSE; inbag <- list(c(0,0)); use.inbag <- FALSE; regularization.factor <- c(0, 0); use.regularization.factor <- FALSE; 
  regularization.usedepth <- FALSE; snp.data <- as.matrix(0); gwa.mode <- FALSE
  
  ## 1.) Handling inputs and checks.
  forest <- object
  if (!inherits(forest, "ocf.forest")) stop("Invalid 'object'.", call. = FALSE) 
  if (is.null(forest$num.trees) || is.null(forest$child.nodeIDs) || is.null(forest$split.varIDs) ||
      is.null(forest$split.values) || is.null(forest$covariate.names) || is.null(forest$treetype)) stop("Invalid 'object'.", call. = FALSE)
  
  if (type == "response") {
    prediction.type <- 1
  } else if (type == "terminalNodes") {
    prediction.type <- 2
  } else {
    stop("Invalid 'type'.", call. = FALSE)
  }
  
  ## 2.) Data.
  # 2a.) Check.
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

  ## 3.) Calling ocf in prediction mode.
  result <- ocfCpp(treetype, x, y, forest$covariate.names, mtry,
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
  
  ## 4.) Handling output.
  if (is.list(result$predictions)) result$predictions <- do.call(rbind, result$predictions)
  
  if (type == "terminalNodes") {
    if (is.vector(result$predictions)) {
      result$predictions <- matrix(result$predictions, nrow = 1)
    }
  }
  
  ##5.) Output.
  return(result)
}


#' Summary Method for ocf Objects
#' 
#' Summarizes an \code{\link{ocf}} object.
#' 
#' @param object An \code{\link{ocf}} object.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Summarizes an \code{\link{ocf}} object.
#' 
#' @examples 
#' ## Load data from orf package.
#' set.seed(1986)
#' 
#' library(orf)
#' data(odata)
#' odata <- odata[1:200, ] # Subset to reduce elapsed time.
#' 
#' y <- as.numeric(odata[, 1])
#' X <- as.matrix(odata[, -1])
#' 
#' ## Fit ocf.
#' forests <- ocf(y, X)
#' 
#' ## Summary.
#' summary(forests)
#' 
#' @seealso \code{\link{ocf}}, \code{\link{marginal_effects}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
summary.ocf <- function(object, ...) {
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
  cat("Honest fraction:   ", object$tuning.info$honesty.fraction)
}
 
 
#' Print Method for ocf Objects
#'
#' Prints an \code{\link{ocf}} object.
#'
#' @param x An \code{\link{ocf}} object.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Prints an \code{\link{ocf}} object.
#' 
#' @examples 
#' \donttest{
#' ## Load data from orf package.
#' set.seed(1986)
#' 
#' library(orf)
#' data(odata)
#' odata <- odata[1:200, ] # Subset to reduce elapsed time.
#' 
#' y <- as.numeric(odata[, 1])
#' X <- as.matrix(odata[, -1])
#' 
#' ## Fit ocf.
#' forests <- ocf(y, X)
#' 
#' ## Print.
#' print(forests)}
#' 
#' @seealso \code{\link{ocf}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
print.ocf <- function(x, ...) {
  summary.ocf(x, ...)
}


#' Summary Method for ocf.marginal Objects
#'
#' Summarizes an \code{ocf.marginal} object.
#'
#' @param object An \code{ocf.marginal} object.
#' @param latex If \code{TRUE}, prints LATEX code.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Summarizes an \code{ocf.marginal} object.
#' 
#' @examples 
#' ## Load data from orf package.
#' set.seed(1986)
#' 
#' library(orf)
#' data(odata)
#' odata <- odata[1:200, ] # Subset to reduce elapsed time.
#' 
#' y <- as.numeric(odata[, 1])
#' X <- as.matrix(odata[, -1])
#' 
#' ## Fit ocf. Use large number of trees.
#' forests <- ocf(y, X, n.trees = 4000)
#' 
#' ## Marginal effects at the mean.
#' me <- marginal_effects(forests, eval = "atmean")
#' summary(me)
#' summary(me, latex = TRUE)
#' 
#' \donttest{
#' ## Add standard errors.
#' honest_forests <- ocf(y, X, n.trees = 4000, honesty = TRUE)
#' honest_me <- marginal_effects(honest_forests, eval = "atmean", inference = TRUE)
#' summary(honest_me, latex = TRUE)}
#' 
#' @details 
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}. If
#' standard errors have been estimated, they are printed in parenthesis below each point estimate.
#' 
#' @seealso \code{\link{ocf}}, \code{\link{marginal_effects}}.
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
summary.ocf.marginal <- function(object, latex = FALSE, ...) {
  if (!(latex %in% c(TRUE, FALSE))) stop("Invalid value of 'latex'.", call. = FALSE)
  
  if (latex) {
    table_names <- rename_latex(rownames(object$marginal.effects))
    
    cat("\\begingroup
    \\setlength{\\tabcolsep}{8pt}
    \\renewcommand{\\arraystretch}{1.1}
    \\begin{table}[H]
        \\centering
        \\begin{adjustbox}{width = 0.75\\textwidth}
        \\begin{tabular}{@{\\extracolsep{5pt}}l", rep(" c", object$n.classes), "}
        \\\\[-1.8ex]\\hline
        \\hline \\\\[-1.8ex]
        &", paste0(" Class ", seq_len(object$n.classes)[-object$n.classes], " &"), paste0(" Class ", object$n.classes), " \\\\
        \\addlinespace[2pt]
        \\hline \\\\[-1.8ex] \n\n", sep = "")
    
    if (inherits(object$standard.errors, "list")) {
      for (i in seq_len(nrow(object$marginal.effects))) {
        cat("        \\texttt{", table_names[i], "} & ", paste0(round(object$marginal.effects[i, 1:(ncol(object$marginal.effects)-1)], 3), " & "), round(object$marginal.effects[i, ncol(object$marginal.effects)], 3), " \\\\ \n", sep = "")
      }
    } else {
      for (i in seq_len(nrow(object$marginal.effects))) {
        cat("        \\texttt{", table_names[i], "} & ", paste0(round(object$marginal.effects[i, 1:(ncol(object$marginal.effects)-1)], 3), " & "), round(object$marginal.effects[i, ncol(object$marginal.effects)], 3), " \\\\ \n", sep = "")
        cat("                     & ", paste0("(", round(object$standard.errors[i, 1:(ncol(object$standard.errors)-1)], 3), ")", " & "), paste0("(", round(object$standard.errors[i, ncol(object$standard.errors)], 3), ")"), " \\\\ \n", sep = "")
      }
    }
    
    cat("\n        \\addlinespace[3pt]
        \\\\[-1.8ex]\\hline
        \\hline \\\\[-1.8ex]
        \\end{tabular}
        \\end{adjustbox}
        \\caption{Marginal effects.}
        \\label{table:ocf.marginal.effects}
    \\end{table}
\\endgroup")
  } else {
    cat("ocf marginal effects results \n\n")
    
    cat("Data info: \n")
    cat("Number of classes:   ", object$n.classes, "\n")
    cat("Sample size:         ", object$n.samples, "\n\n")
    
    cat("Tuning parameters: \n")
    cat("Evaluation:          ", object$evaluation, "\n")
    cat("Bandwidth:           ", object$bandwitdh, "\n")
    cat("Number of trees:     ", object$n.trees, "\n")
    cat("Honest forests:      ", object$honesty, "\n")
    cat("Honesty fraction:    ",  object$honesty.fraction, "\n\n")
    
    cat("Marginal Effects: \n")
    
    print(round(object$marginal.effects, 3))
  }
}


#' Print Method for ocf.marginal Objects
#'
#' Prints an \code{ocf.marginal} object.
#'
#' @param x An \code{ocf.marginal} object.
#' @param latex If \code{TRUE}, prints LATEX code.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Prints an \code{ocf.marginal} object.
#' 
#' @examples 
#' ## Load data from orf package.
#' set.seed(1986)
#' 
#' library(orf)
#' data(odata)
#' odata <- odata[1:200, ] # Subset to reduce elapsed time.
#' 
#' y <- as.numeric(odata[, 1])
#' X <- as.matrix(odata[, -1])
#' 
#' ## Fit ocf. Use large number of trees.
#' forests <- ocf(y, X, n.trees = 4000)
#' 
#' ## Marginal effects at the mean.
#' me <- marginal_effects(forests, eval = "atmean")
#' print(me)
#' print(me, latex = TRUE)
#' 
#' \donttest{
#' ## Add standard errors.
#' honest_forests <- ocf(y, X, n.trees = 4000, honesty = TRUE)
#' honest_me <- marginal_effects(honest_forests, eval = "atmean", inference = TRUE)
#' print(honest_me, latex = TRUE)}
#' 
#' @details 
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}. If
#' standard errors have been estimated, they are printed in parenthesis below each point estimate.
#' 
#' @seealso \code{\link{ocf}}, \code{\link{marginal_effects}}.
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
print.ocf.marginal <- function(x, latex = FALSE, ...) {
  summary.ocf.marginal(x, latex, ...)
}


#' Prediction Method for mml Objects
#'
#' Prediction method for class \code{mml}.
#'
#' @param object An \code{mml} object.
#' @param data Data set of class \code{data.frame}. It must contain the same covariates used to train the base learners. If \code{data} is \code{NULL}, then \code{object$X} is used.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Matrix of predictions.
#' 
#' @examples 
#' ## Load data from orf package.
#' set.seed(1986)
#' 
#' library(orf)
#' data(odata)
#' odata <- odata[1:200, ] # Subset to reduce elapsed time.
#' 
#' y <- as.numeric(odata[, 1])
#' X <- as.matrix(odata[, -1])
#' 
#' ## Training-test split.
#' train_idx <- sample(seq_len(length(y)), floor(length(y) * 0.5))
#' 
#' y_tr <- y[train_idx]
#' X_tr <- X[train_idx, ]
#' 
#' y_test <- y[-train_idx]
#' X_test <- X[-train_idx, ]
#' 
#' ## Fit multinomial machine learning on training sample using two different learners.
#' multinomial_forest <- multinomial_ml(y_tr, X_tr, learner = "forest")
#' multinomial_l1 <- multinomial_ml(y_tr, X_tr, learner = "l1")
#' 
#' ## Predict out of sample.
#' predictions_forest <- predict(multinomial_forest, X_test)
#' predictions_l1 <- predict(multinomial_l1, X_test)
#' 
#' ## Compare predictions.
#' cbind(head(predictions_forest), head(predictions_l1))
#' 
#' @details
#' If \code{object$learner == "l1"}, then \code{data} is scaled to have zero mean and unit variance.
#' 
#' @seealso \code{\link{multinomial_ml}}, \code{\link{ordered_ml}}
#' 
#' @importFrom stats predict
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
predict.mml <- function(object, data = NULL, ...) {
  ## 0.) Handling inputs and checks.
  if (is.null(data)) data <- object$X
  learner <- object$learner
  estimators <- object$estimators
  n_categories <- length(unique(object$y))
  
  ## 1.) Get predictions.
  if (learner == "forest") {
    predictions <- lapply(estimators, function(x) {predict(x, data)$predictions}) 
  } else if (learner == "l1") {
    data_design <- stats::model.matrix(y ~ ., data = data.frame("y" = 1, data))[, -1]
    data_scaled <- as.matrix(scale(data_design))
    predictions <- lapply(estimators, function(model) {as.numeric(predict(model, data_scaled, type = "response"))}) 
  }
  
  ## 2.) Normalize and put into matrix.
  predictions_final <- sapply(predictions, function(x) as.matrix(x))
  predictions_final <- matrix(apply(predictions_final, 1, function(x) (x) / (sum(x))), ncol = n_categories, byrow = T)
  colnames(predictions_final) <- paste0("P(Y=", seq_len(n_categories), ")")
  
  ## 3.) Output.
  return(predictions_final)
}


#' Prediction Method for oml Objects
#'
#' Prediction method for class \code{oml}.
#'
#' @param object An \code{oml} object.
#' @param data Data set of class \code{data.frame}. It must contain the same covariates used to train the base learners. If \code{data} is \code{NULL}, then \code{object$X} is used.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Matrix of predictions.
#' 
#' @examples 
#' ## Load data from orf package.
#' set.seed(1986)
#' 
#' library(orf)
#' data(odata)
#' odata <- odata[1:200, ] # Subset to reduce elapsed time.
#' 
#' y <- as.numeric(odata[, 1])
#' X <- as.matrix(odata[, -1])
#' 
#' ## Training-test split.
#' train_idx <- sample(seq_len(length(y)), floor(length(y) * 0.5))
#' 
#' y_tr <- y[train_idx]
#' X_tr <- X[train_idx, ]
#' 
#' y_test <- y[-train_idx]
#' X_test <- X[-train_idx, ]
#' 
#' ## Fit ordered machine learning on training sample using two different learners.
#' ordered_forest <- ordered_ml(y_tr, X_tr, learner = "forest")
#' ordered_l1 <- ordered_ml(y_tr, X_tr, learner = "l1")
#' 
#' ## Predict out of sample.
#' predictions_forest <- predict(ordered_forest, X_test)
#' predictions_l1 <- predict(ordered_l1, X_test)
#' 
#' ## Compare predictions.
#' cbind(head(predictions_forest), head(predictions_l1))
#' 
#' @details
#' If \code{object$learner == "l1"}, then \code{data} is scaled to have zero mean and unit variance.
#' 
#' @seealso \code{\link{multinomial_ml}}, \code{\link{ordered_ml}}
#' 
#' @importFrom stats predict
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
predict.oml <- function(object, data = NULL, ...) {
  ## 0.) Handling inputs and checks.
  if (is.null(data)) data <- object$X
  learner <- object$learner
  estimators <- object$estimators
  n_categories <- length(unique(object$y))
  n <- dim(data)[1]
  
  ## 1.) Get predictions.
  if (learner == "forest") {
    predictions <- lapply(estimators, function(x) {predict(x, data)$predictions}) 
  } else if (learner == "l1") {
    data_design <- stats::model.matrix(y ~ ., data = data.frame("y" = 1, data))[, -1]
    data_scaled <- as.matrix(scale(data_design))
    predictions <- lapply(estimators, function(model) {as.numeric(predict(model, data_scaled, type = "response"))}) 
  }
  
  ## 2.) Pick differences.
  predictions1 <- append(predictions, list(rep(1, n))) 
  predictions0 <- append(list(rep(0, n)), predictions) 
  differences <- as.list(mapply(function(x, y) x - y, predictions1, predictions0, SIMPLIFY = FALSE))
  
  ## 3.) Truncate, normalize and put into matrix.
  predictions_final <- lapply(differences, function(x) ifelse((x < 0), 0, x))
  predictions_final <- sapply(predictions_final, function(x) as.matrix(x))
  predictions_final <- matrix(apply(predictions_final, 1, function(x) (x) / (sum(x))), ncol = n_categories, byrow = T)
  colnames(predictions_final) <- paste0("P(Y=", seq_len(n_categories), ")")
  
  ## 4.) Output.
  return(predictions_final)
}
