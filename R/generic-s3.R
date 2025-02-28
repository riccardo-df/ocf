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
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#' 
#' data <- generate_ordered_data(100)
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
#' @importFrom stats predict
#' 
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2025). Ordered Correlation Forest. Econometric Reviews, 1–17. \doi{10.1080/07474938.2024.2429596}.
#' }
#'
#' @seealso \code{\link{ocf}}, \code{\link{marginal_effects}}
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
        honest_outcomes[[counter]] <- data.frame("y_m_honest" = ifelse(object$honest_data$Y <= m, 1, 0), "y_m_1_honest" = ifelse(object$honest_data$Y <= m -1, 1, 0))
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
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2025). Ordered Correlation Forest. Econometric Reviews, 1–17. \doi{10.1080/07474938.2024.2429596}.
#' }
#' 
#' @seealso \code{\link{ocf}}, \code{\link{marginal_effects}}
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
#' ## Summary.
#' summary(forests)}
#' 
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2025). Ordered Correlation Forest. Econometric Reviews, 1–17. \doi{10.1080/07474938.2024.2429596}.
#' }
#'
#' @seealso \code{\link{ocf}}, \code{\link{marginal_effects}}
#' 
#' @export
summary.ocf <- function(object, ...) {
  cat("Call: \n")
  cat(deparse(object$tuning.info$call), "\n\n")
  
  cat("Data info: \n")
  cat("Full sample size:  ", dim(object$full_data)[1], "\n")
  cat("N. covariates:     ", dim(object$full_data)[2] - 1, "\n")
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
#' ## Print.
#' print(forests)}
#' 
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2025). Ordered Correlation Forest. Econometric Reviews, 1–17. \doi{10.1080/07474938.2024.2429596}.
#' }
#'
#' @seealso \code{\link{ocf}}
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
#' summary(me)
#' summary(me, latex = TRUE)
#' 
#' ## Add standard errors.
#' honest_forests <- ocf(Y, X, honesty = TRUE)
#' honest_me <- marginal_effects(honest_forests, eval = "atmean", inference = TRUE)
#' summary(honest_me, latex = TRUE)}
#' 
#' @details 
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}. If
#' standard errors have been estimated, they are printed in parenthesis below each point estimate.
#' 
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2025). Ordered Correlation Forest. Econometric Reviews, 1–17. \doi{10.1080/07474938.2024.2429596}.
#' }
#'
#' @seealso \code{\link{ocf}}, \code{\link{marginal_effects}}.
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
#' print(me)
#' print(me, latex = TRUE)
#' 
#' ## Add standard errors.
#' honest_forests <- ocf(Y, X, honesty = TRUE)
#' honest_me <- marginal_effects(honest_forests, eval = "atmean", inference = TRUE)
#' print(honest_me, latex = TRUE)}
#' 
#' @details 
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}. If
#' standard errors have been estimated, they are printed in parenthesis below each point estimate.
#' 
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2025). Ordered Correlation Forest. Econometric Reviews, 1–17. \doi{10.1080/07474938.2024.2429596}.
#' }
#'
#' @seealso \code{\link{ocf}}, \code{\link{marginal_effects}}.
#' 
#' @export
print.ocf.marginal <- function(x, latex = FALSE, ...) {
  summary.ocf.marginal(x, latex, ...)
}


#' Plot Method for ocf.marginal Objects
#'
#' Plots an \code{ocf.marginal} object.
#'
#' @param x An \code{ocf.marginal} object.
#' @param class_names Character vector of length equal to \code{x$n.classes} to set class names in the plot.
#' @param point_size Controls the points' size.
#' @param facet_text_size Controls the facets' labels' size.
#' @param legend_text_size Controls the legends' size.
#' 
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Plots an \code{ocf.marginal} object.
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
#' plot(me)
#' 
#' ## Add standard errors.
#' honest_forests <- ocf(Y, X, honesty = TRUE)
#' honest_me <- marginal_effects(honest_forests, eval = "atmean", inference = TRUE)
#' plot(honest_me)}
#' 
#' @details 
#' If standard errors have been estimated, 95\% confidence intervals are shown.
#' 
#' @import tidyr ggplot2 ggthemes
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom magrittr %>%
#' 
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2025). Ordered Correlation Forest. Econometric Reviews, 1–17. \doi{10.1080/07474938.2024.2429596}.
#' }
#'
#' @seealso \code{\link{ocf}}, \code{\link{marginal_effects}}.
#' 
#' @export
plot.ocf.marginal <- function(x, class_names = NULL, point_size = 2, facet_text_size = 12, legend_text_size = 10, ...) {
  ## Handling inputs and checks.
  CI_lower <- NULL
  CI_upper <- NULL 
  covariate <- NULL 
  marginal_effect <- NULL 
  standard_error <- NULL
  
  n_classes <- x$n.classes
  
  if (is.null(class_names)) {
    class_names <- paste0("Class ", seq_len(n_classes))
  } else {
    if (!is.character(class_names) || length(class_names) != n_classes) stop(paste("class_names must be a character vector of length", n_classes))
  }
  
  ## Pivot longer for marginal effects and standard errors (latter only if honesty is TRUE).
  long_me <- x$marginal.effects %>%
    as.data.frame() %>%
    dplyr::mutate(covariate = rownames(x$marginal.effects)) %>%
    tidyr::pivot_longer(cols = starts_with("P"), names_to = "class", values_to = "marginal_effect")
  
  if (x$honesty) {
    long_se <- x$standard.errors %>%
      as.data.frame() %>%
      dplyr::mutate(covariate = rownames(x$standard.errors)) %>%
      tidyr::pivot_longer(cols = starts_with("P"), names_to = "class", values_to = "standard_error")
  }

  ## Arrange plotting data. If honesty is FALSE, set standard errors to zero. Construct 95% CIs.
  if (x$honesty) {
    plot_dta <- long_me %>%
      dplyr::left_join(long_se, by = c("covariate", "class"))
  } else {
    plot_dta <- long_me %>%
      dplyr::mutate(standard_error = 0)
  }
  
  plot_dta <- plot_dta %>%
    dplyr::mutate(CI_upper = marginal_effect + 1.96 * standard_error,
                  CI_lower = marginal_effect - 1.96 * standard_error)
  
  ## Rename classes.
  for (m in seq_len(n_classes)) {
    plot_dta$class[grepl(m, plot_dta$class)] <- class_names[m]
  }
  
  ## Generate plot.
  plot_dta %>%
    dplyr::mutate(class = factor(class, levels = class_names),
                  class_reversed = factor(class, levels = rev(levels(class)))) %>%
    ggplot2::ggplot(ggplot2::aes(x = marginal_effect, y = interaction(class_reversed, covariate), color = class)) +
    ggplot2::geom_point(size = point_size, shape = 4, position = ggplot2::position_dodge(width = 0.7)) +
    ggplot2::geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, position = ggplot2::position_dodge(width = 0.7)) + 
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::facet_grid(covariate ~ ., switch = "y", scales = "free_y", space = "free_y") +
    ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::ggtitle("") +
    ggthemes::theme_tufte() + 
    ggplot2::theme(legend.position = "right", 
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = legend_text_size),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   strip.text = element_text(face = "italic", size = facet_text_size),
                   strip.background = element_rect(fill = "gray90", color = "black"))
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
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#' 
#' data <- generate_ordered_data(100)
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
#' ## Fit multinomial machine learning on training sample using two different learners.
#' multinomial_forest <- multinomial_ml(Y_tr, X_tr, learner = "forest")
#' multinomial_l1 <- multinomial_ml(Y_tr, X_tr, learner = "l1")
#' 
#' ## Predict out of sample.
#' predictions_forest <- predict(multinomial_forest, X_test)
#' predictions_l1 <- predict(multinomial_l1, X_test)
#' 
#' ## Compare predictions.
#' cbind(head(predictions_forest), head(predictions_l1))}
#' 
#' @details
#' If \code{object$learner == "l1"}, then \code{\link[stats]{model.matrix}} is used to handle non-numeric covariates. If we also
#' have \code{object$scaling == TRUE}, then \code{data} is scaled to have zero mean and unit variance.
#' 
#' @importFrom stats predict
#' 
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2025). Ordered Correlation Forest. Econometric Reviews, 1–17. \doi{10.1080/07474938.2024.2429596}.
#' }
#'
#' @seealso \code{\link{multinomial_ml}}, \code{\link{ordered_ml}}
#' 
#' @export
predict.mml <- function(object, data = NULL, ...) {
  ## 0.) Handling inputs and checks.
  if (is.null(data)) data <- object$X
  learner <- object$learner
  estimators <- object$estimators
  scale <- object$scaling
  n_categories <- length(unique(object$Y))
  
  ## 1.) Get predictions.
  if (learner == "forest") {
    predictions <- lapply(estimators, function(x) {predict(x, data)$predictions}) 
  } else if (learner == "l1") {
    data_design <- stats::model.matrix(y ~ ., data = data.frame("y" = 1, data))[, -1]
    if (scale) data_design <- as.matrix(scale(data_design))
    predictions <- lapply(estimators, function(x) {as.numeric(predict(x, data_design, s = "lambda.min", type = "response"))}) 
  }
  
  ## 2.) Put into matrix and normalize.
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
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#' 
#' data <- generate_ordered_data(100)
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
#' ## Fit ordered machine learning on training sample using two different learners.
#' ordered_forest <- ordered_ml(Y_tr, X_tr, learner = "forest")
#' ordered_l1 <- ordered_ml(Y_tr, X_tr, learner = "l1")
#' 
#' ## Predict out of sample.
#' predictions_forest <- predict(ordered_forest, X_test)
#' predictions_l1 <- predict(ordered_l1, X_test)
#' 
#' ## Compare predictions.
#' cbind(head(predictions_forest), head(predictions_l1))}
#' 
#' @details
#' If \code{object$learner == "l1"}, then \code{\link[stats]{model.matrix}} is used to handle non-numeric covariates. If we also
#' have \code{object$scaling == TRUE}, then \code{data} is scaled to have zero mean and unit variance.
#' 
#' @importFrom stats predict
#' 
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2025). Ordered Correlation Forest. Econometric Reviews, 1–17. \doi{10.1080/07474938.2024.2429596}.
#' }
#'
#' @seealso \code{\link{multinomial_ml}}, \code{\link{ordered_ml}}
#' 
#' @export
predict.oml <- function(object, data = NULL, ...) {
  ## 0.) Handling inputs and checks.
  if (is.null(data)) data <- object$X
  learner <- object$learner
  estimators <- object$estimators
  scale <- object$scaling
  n_categories <- length(unique(object$Y))
  n <- dim(data)[1]
  
  ## 1.) Get predictions.
  if (learner == "forest") {
    predictions <- lapply(estimators, function(x) {predict(x, data)$predictions}) 
  } else if (learner == "l1") {
    data_design <- stats::model.matrix(y ~ ., data = data.frame("y" = 1, data))[, -1]
    if (scale) data_design <- as.matrix(scale(data_design))
    predictions <- lapply(estimators, function(x) {as.numeric(predict(x, data_design, s = "lambda.min", type = "response"))}) 
  }
  
  ## 2.) Pick differences.
  predictions1 <- append(predictions, list(rep(1, n))) 
  predictions0 <- append(list(rep(0, n)), predictions) 
  differences <- as.list(mapply(function(x, y) x - y, predictions1, predictions0, SIMPLIFY = FALSE))
  
  ## 3.) Truncate, put into matrix, and normalize.
  predictions_final <- lapply(differences, function(x) ifelse((x < 0), 0, x))
  predictions_final <- sapply(predictions_final, function(x) as.matrix(x))
  predictions_final <- matrix(apply(predictions_final, 1, function(x) (x) / (sum(x))), ncol = n_categories, byrow = T)
  colnames(predictions_final) <- paste0("P(Y=", seq_len(n_categories), ")")
  
  ## 4.) Output.
  return(predictions_final)
}
