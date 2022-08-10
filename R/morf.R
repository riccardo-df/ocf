##' Modified Ordered Random Forest for Estimating the Ordered Choice Model
##' 
##' Non-parametric estimation of the ordered choice model using random forests.
##'
##' @param x Covariate matrix (no intercept).
##' @param y Outcome vector.
##' @param n.trees Number of trees.
##' @param mtry Number of covariates to possibly split at in each node. Default is the (rounded down) square root of the number of covariates. Alternatively, one can pass a single-argument function returning an integer, where the argument is the number of covariates.
##' @param write.forest Save \code{morf.forest} object, required for prediction. Set to \code{FALSE} to reduce memory usage if no prediction intended.
##' @param min.node.size Minimal node size.
##' @param max.depth Maximal tree depth. A value of 0 (the default) corresponds to unlimited depth, 1 to "stumps" (one split per tree).
##' @param replace If \code{TRUE}, grow trees on bootstrap subsamples. Otherwise, trees are grown on random subsamples drawn without replacement. 
##' @param sample.fraction Fraction of observations to sample. Default is 1 for bootstrap sampling and 0.632 for sampling without replacement. 
##' @param case.weights Weights for sampling of training observations. Observations with larger weights will be drawn with higher probability in the subsamples for the trees.
##' @param split.select.weights Numeric vector with weights between 0 and 1, used to calculate the probability to select variables for splitting. Alternatively, a list of size \code{n.trees}, containing split select weight vectors for each tree can be used.  
##' @param always.split.variables Character vector with variable names to be always selected in addition to the \code{mtry} variables tried for splitting.
##' @param keep.inbag Save how often observations are in-bag in each tree. 
##' @param inbag Manually set observations per tree. List of size \code{n.trees}, containing in-bag counts for each observation. Can be used for stratified sampling.
##' @param holdout Hold-out mode. Hold-out all samples with zero \code{case.weights} and use these for variable importance and prediction error.
##' @param oob.error Compute out-of-bag prediction error. Set to \code{FALSE} to save computation time.
##' @param n.threads Number of threads. Default is number of CPUs available.
##' @param save.memory Use memory saving (but slower) splitting mode. Warning: This option slows down the tree growing, use only if you encounter memory problems.
##' @param verbose Show computation status and estimated runtime.
##' @param seed Random seed. Default is \code{NULL}, which generates the seed from \code{R}. Set to \code{0} to ignore the \code{R} seed. 
##' 
##' @details 
##' \subsection{Splitting Criterion}{
##' \code{morf} fits \code{M} random forests, where \code{M} is the number of classes of \code{y}. 
##' Each forest computes the conditional probabilities of the \code{m}-th class:
##' 
##' \deqn{p_m (x) = P ( Y_i = m \, | \, X_i = x), \,\,\,  m = 1, ..., M}
##' 
##' To estimate this quantity, \code{morf} exploits the following:
##' 
##' \deqn{p_m (x) = E [\, 1 (Y_i \leq m) \, | \, X_i = x] - E[\, 1 ( Y_i \leq m - 1 ) \, | \, X_i = x]  
##'                 = \mu_m (x) - \mu_{m-1} (x)}
##' 
##' with \code{1(.)} an indicator of the truth of its argument. It would be tempting to estimate the 
##' conditional probabilities as:
##' 
##' \deqn{\hat{p}_m (x) = \hat{\mu}_m (x) - \hat{\mu}_{m-1} (x)}
##' 
##' where the expectations can be estimated via separate regression forests. However, this strategy 
##' ignores potential correlation in the estimation errors of the two surfaces, as shown by the following expansion:
##' 
##' \deqn{MSE[\hat{p}_m (x)] = MSE[\hat{\mu}_m (x)] + MSE [\hat{\mu}_{m-1} (x)] - 2 MCE[\hat{\mu}_m (x), \hat{\mu}_{m-1} (x)]}
##' 
##' where the last term is the mean correlation error of the estimation. The splitting criterion used by \code{morf} 
##' seeks to greedily minimize the above expression, thereby accounting for the correlation term. Therefore, each split
##' is chosen to minimize:
##' 
##' \deqn{Var( 1 (Y_i \leq m)) + Var( 1 (Y_i \leq m - 1)) - 2 * Cor(1 (Y_i \leq m), 1 (Y_i \leq m - 1))}
##' }
##' 
##' \subsection{Predictions}{
##' Predictions in the \code{l}-th leaf are computed as:
##' 
##' \deqn{\frac{1}{\{i: x_i \in L_l\}} \sum_{\{i: x_i \in L_l\}} 1 (Y_i \leq m) - \frac{1}{\{i: x_i \in L_l\}} \sum_{\{i: x_i \in L_l\}} 1 (Y_i \leq m - 1)}
##' 
##' Notice that a normalization step may be needed to ensure that the estimated probabilities sum up to one across classes.
##' }
##' 
##' \subsection{Variable Importance}{
##' For each covariate, an overall variable importance measure is provided. The \code{m}-th forest computes the 
##' importance of the j-th covariate for the \code{m}-th class by recording the improvement in the splitting criterion 
##' at each split placed on such covariate. Summing over all such splits of the trees in the \code{m}-th 
##' forest gives the j-th covariate's importance for the \code{m}-th class. The overall variable importance measure of this 
##' covariate is then defined as the mean of its importances in each class.
##' }
##' 
##' @return 
##' Object of class \code{morf} with elements:
##'   \item{\code{forest.1}}{\code{morf.forest} object of the first class.}
##'   \item{\code{forest.2}}{\code{morf.forest} object of the second class.} 
##'   \item{\code{...}}{}
##'   \item{\code{forest.M}}{\code{morf.forest} object of the last class.}
##'   \item{\code{predictions}}{Matrix of predicted conditional class probabilities, based on out-of-bag samples. Requires \code{oob.error = TRUE} (the default).}
##'   \item{\code{mean.squared.error}}{Mean squared error of the model, based on out-of-bag predictions. Requires \code{oob.error = TRUE} (the default).}
##'   \item{\code{mean.ranked.score}}{Mean squared error of the model, based on out-of-bag predictions. Requires \code{oob.error = TRUE} (the default).}
##'   \item{\code{overall.importance}}{Measure of overall variable importance, computed as the mean of the importance of each variable across classes. Relative importance is provided.}
##'   \item{\code{n.classes}}{Number of classes.}
##'   \item{\code{n.samples}}{Number of observations.}
##'   \item{\code{n.covariates}}{Number of covariates.}
##'   \item{\code{n.trees}}{Number of trees of each forest.}
##'   \item{\code{mtry}}{Number of covariates considered for splitting at each step.}
##'   \item{\code{min.node.size}}{Minimum node size.}
##'   \item{\code{replace}}{Whether the subsamples to grow trees are drawn with replacement (bootstrap).}
##'   \item{\code{data}}{Training sample.}
##'   \item{\code{call}}{The system call.}
##' 
##' @importFrom Rcpp evalCpp
##' @import utils
##' @useDynLib morf
##' 
##' @seealso \code{\link{predict.morf}}, \code{\link{marginal_effects}}, \code{\link{mean_squared_error}},
##' \code{\link{mean_ranked_score}}.
##' 
##' @author Riccardo Di Francesco
##' 
##' @references
##' \itemize{
##'   \item Lechner, M., & Okasa, G. (2019). Random forest estimation of the ordered choice model. arXiv preprint arXiv:1907.02436. \doi{10.48550/arXiv.1907.02436}.
##'   \item Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. \doi{10.18637/jss.v077.i01}.
##' }
##' @export
morf <- function(x = NULL, y = NULL, 
                 n.trees = 1000, mtry = NULL,
                 min.node.size = 5, max.depth = 0, 
                 replace = FALSE, sample.fraction = ifelse(replace, 1, 0.632), case.weights = NULL,
                 split.select.weights = NULL, always.split.variables = NULL,
                 keep.inbag = FALSE, inbag = NULL, holdout = FALSE, oob.error = TRUE,
                 n.threads = NULL, save.memory = FALSE,
                 write.forest = TRUE, verbose = TRUE, seed = NULL) {
  ## Handling inputs and checks.
  # Outcomes and covariates.
  if (is.null(x)) stop("Input x is required.", call. = FALSE)
  if (is.null(y)) stop("Input y is required.", call. = FALSE)

  # Missing values.
  if (any(is.na(x))) {
    offending.columns <- colnames(x)[colSums(is.na(x)) > 0]
    stop("'x' has missing values in columns: ", paste0(offending.columns, collapse = ", "), ".", call. = FALSE)
  } else if (any(is.na(y))) {
    stop("'y' has missing values.", call. = FALSE)
  }

  # Sparse matrix.
  if (inherits(x, "Matrix")) {
    if (!inherits(x, "dgCMatrix")) stop("Error: Currently only sparse data of class 'dgCMatrix' supported.", call. = FALSE)
  }

  # Recoding characters as factors.
  covariate.levels <- NULL
  if (!is.matrix(x) && !inherits(x, "Matrix") && ncol(x) > 0) {
    character.idx <- sapply(x, is.character)
    x[character.idx] <- lapply(x[character.idx], factor)
    
    if (any(sapply(x, is.factor))) {
      covariate.levels <- lapply(x, levels)
    }     
  }

  # Error if no covariates.
  independent.variable.names <- colnames(x)
  all.independent.variable.names <- independent.variable.names
  if (length(all.independent.variable.names) < 1) stop("No covariates found.", call. = FALSE)

  # Number of trees.
  if (!is.numeric(n.trees) || n.trees < 1) stop("Invalid value for 'n.trees'.", call. = FALSE)

  # mtry.
  if (is.function(mtry)) {
    nv <- length(independent.variable.names)

    if (length(formals(mtry)) > 1) stop("'mtry' function requires single argument (the number of covariates in the model).", call. = FALSE)
  
    mtry <- try(mtry(nv), silent = TRUE)

    if (inherits(mtry, "try-error")) {
      message("The 'mtry' function produced the error: ", mtry)
      stop("'mtry' function evaluation resulted in an error.", call. = FALSE)
    }

    if (!is.numeric(mtry) || length(mtry) != 1) {
      stop("'mtry' function should return a single integer or numeric.")
    } else {
      mtry <- as.integer(mtry)
    }

    if (mtry < 1 || mtry > nv) {
      stop("'mtry' function should evaluate to a value not less than 1 and not greater than the number of covariates ( = ", nv, " )", call. = FALSE)
    }
  }

  if (is.null(mtry)) {
    mtry <- 0
  } else if (!is.numeric(mtry) || mtry < 0) {
    stop("Invalid value for 'mtry'")
  }

  # Seed.
  if (is.null(seed)) seed <- stats::runif(1 , 0, .Machine$integer.max)
  
  # Keep in-bag.
  if (!is.logical(keep.inbag)) stop("Invalid value for 'keep.inbag'", call. = FALSE)
  
  # Number of threads.
  if (is.null(n.threads)) {
    n.threads <- 0
  } else if (!is.numeric(n.threads) || n.threads < 0) {
    stop("Invalid value for 'n.threads'", call. = FALSE)
  }

  # Minimum node size.
  if (!is.numeric(min.node.size) || min.node.size <= 0) stop("Invalid value for 'min.node.size'", call. = FALSE)
  
  # Tree depth.
  if (!is.numeric(max.depth) || max.depth < 0) stop("Invalid value for 'max.depth'. Please give a positive integer.", call. = FALSE)

  # Sample fraction.
  if (!is.numeric(sample.fraction)) stop("Invalid value for 'sample.fraction'. Please give a value in (0,1].", call. = FALSE)
  if (sample.fraction <= 0 || sample.fraction > 1) stop("Invalid value for 'sample.fraction' Please give a value in (0,1].", call. = FALSE)

  # Case weights: NULL for no weights or all weights equal.
  if (is.null(case.weights) || length(unique(case.weights)) == 1) {
    case.weights <- c(0, 0)
    use.case.weights <- FALSE
    
    if (holdout) stop("Case weights required to use holdout mode.", call. = FALSE)
  } else {
    use.case.weights <- TRUE

    if (holdout) sample.fraction <- sample.fraction * mean(case.weights > 0)
    if (!replace && sum(case.weights > 0) < sample.fraction * nrow(x)) stop("Fewer non-zero case weights than observations to sample.", call. = FALSE)
  }

  # Manual in-bag selection.
  if (is.null(inbag)) {
    inbag <- list(c(0 ,0))
    use.inbag <- FALSE
  } else if (is.list(inbag)) {
    use.inbag <- TRUE
    
    if (use.case.weights) stop("Combination of 'case.weights' and 'inbag' not supported.", call. = FALSE)
    if (length(sample.fraction) > 1) stop("Combination of class-wise sampling and inbag not supported.", call. = FALSE)
    if (length(inbag) != n.trees) stop("Size of inbag list not equal to the number of trees.", call. = FALSE)
  } else {
    stop("Invalid inbag, expects list of vectors of size n.trees.", call. = FALSE)
  }

  # Split select weights: NULL for no weights.
  if (is.null(split.select.weights)) {
    split.select.weights <- list(c(0, 0))
    use.split.select.weights <- FALSE
  } else if (is.numeric(split.select.weights)) {
    if (length(split.select.weights) != length(independent.variable.names)) stop("Number of split select weights not equal to number of covariates.", call. = FALSE)
    split.select.weights <- list(split.select.weights)
    use.split.select.weights <- TRUE
  } else if (is.list(split.select.weights)) {
    if (length(split.select.weights) != n.trees) stop("Size of split select weights list not equal to number of trees.", call. = FALSE)
    use.split.select.weights <- TRUE
  } else {
    stop("Invalid split select weights.", call. = FALSE)
  }

  # Always split variables: NULL for no variables.
  if (is.null(always.split.variables)) {
    always.split.variables <- c("0", "0")
    use.always.split.variables <- FALSE
  } else {
    use.always.split.variables <- TRUE
  }

  # Use sparse matrix.
  if (inherits(x, "dgCMatrix")) {
    sparse.x <- x
    x <- matrix(c(0, 0))
    use.sparse.data <- TRUE
  } else {
    sparse.x <- Matrix::Matrix(matrix(c(0, 0)))
    use.sparse.data <- FALSE
    if (is.data.frame(x)) {
      x <- data.matrix(x)
    }
  }

  # Defaults for variables not needed.
  splitrule.num <- 1
  treetype <- 3
  probability <- FALSE
  
  importance.mode <- 1
  scale.permutation.importance <- FALSE
  local.importance <- FALSE
  
  respect.unordered.factors <- "ignore" 
  unordered.factor.variables <- c("0", "0")
  use.unordered.factor.variables <- FALSE
  order.snps <- FALSE
  
  prediction.mode <- FALSE
  predict.all <- FALSE
  prediction.type <- 1
  
  loaded.forest <- list() 
  
  snp.data <- as.matrix(0)
  
  class.weights <- rep(1, nlevels(y))
  alpha <- 0.5
  minprop <- 0.1
  num.random.splits <- 1
  
  use.regularization.factor <- FALSE
  regularization.factor <- c(0, 0)
  regularization.usedepth <- FALSE

  # Preparing for loop. 
  y.classes <- sort(unique(y))
  n.classes <- length(y.classes)
  
  class.probabilities <- matrix(NA, nrow = length(y), ncol = n.classes)
  overall.importance <- 0
  
  results <- list()
  
  ## Growing an ordered forest for each class.
  for (m in seq_len(n.classes)) {
    if (verbose) cat(crayon::silver("Growing forest for class "), crayon::silver(m), 
                     crayon::silver("/"), crayon::silver(n.classes), crayon::silver("\n"), sep = "")
    
    y_m_1 <- ifelse(y <= m - 1, 1, 0)
    y_m <- ifelse(y <= m, 1, 0)
    
    y.mat <- matrix(c(as.numeric(y_m_1), as.numeric(y_m)), ncol = 2) 
    
    class_result <- morfCpp(treetype, x, y.mat, independent.variable.names, mtry,
                     n.trees, verbose, seed, n.threads, write.forest, importance.mode,
                     min.node.size, split.select.weights, use.split.select.weights,
                     always.split.variables, use.always.split.variables,
                     prediction.mode, loaded.forest, snp.data,
                     replace, probability, unordered.factor.variables, use.unordered.factor.variables, 
                     save.memory, splitrule.num, case.weights, use.case.weights, class.weights, 
                     predict.all, keep.inbag, sample.fraction, alpha, minprop, holdout, prediction.type, 
                     num.random.splits, sparse.x, use.sparse.data, order.snps, oob.error, max.depth, 
                     inbag, use.inbag, 
                     regularization.factor, use.regularization.factor, regularization.usedepth)
    
    if (length(class_result) == 0) stop("User interrupt or internal error.", call. = FALSE)
    
    ## Handling class output.
    # Names for variable importance.
    names(class_result$variable.importance) <- all.independent.variable.names
    overall.importance <- overall.importance + class_result$variable.importance
  
    # Writing forest object.
    if (write.forest) {
      if (is.factor(y)) {
        class_result$forest$levels <- levels(y)
      }
      
      class_result$forest$covariate.names <- independent.variable.names
      class_result$forest$treetype <- "ordered"
      
      if (!is.null(covariate.levels)) {
        class_result$forest$covariate.levels <- covariate.levels
      }
    }
    
    # Saving and removing unnecessary element.
    if (oob.error) class.probabilities[, m] <- class_result$predictions
    n.trees <- class_result$num.trees
    num.covariates <- class_result$num.covariates
    mtry <- class_result$mtry
    min.node.size <- class_result$min.node.size
    
    class_result <- class_result[names(class_result) %in% c("forest")]
    
    # Object class.
    class(class_result$forest) <- "morf.forest"
    
    # Storing for output.
    results[[m]] <- class_result$forest
    names(results)[m] <- paste("forest", m, sep = ".")
  }
  
  ## Normalization step.
  class.probabilities <- matrix(apply(class.probabilities, 1, function(x) (x)/(sum(x))), ncol = n.classes, byrow = T)
  colnames(class.probabilities) <- sapply(y.classes, function(x) paste("class", x, sep = "."))
  
  ## Overall variable importance.
  overall.importance <- overall.importance / n.classes
  overall.importance <- overall.importance / sum(overall.importance)
  
  ## Handling output.
  # Predictions.
  if (oob.error) {
    results$predictions <- class.probabilities
  } else{
    results$predictions <- list()
  }
  
  # Prediction error.
  results$mean.squared.error <- mean_squared_error(y, class.probabilities)
  results$mean.ranked.score <- mean_ranked_score(y, class.probabilities)
  
  # Overall variable importance.
  results$overall.importance <- overall.importance
  
  # Number of classes.
  results$n.classes <- n.classes
  
  # Number of units.
  if (use.sparse.data) {
    results$n.samples <- nrow(sparse.x)
  } else {
    results$n.samples <- nrow(x)
  }
  
  # Number of covariates.
  results$n.covariates <- num.covariates
  
  # Number of trees.
  results$n.trees <- n.trees
  
  # Mtry.
  results$mtry <- mtry
  
  # Minimum node size.
  results$min.node.size <- min.node.size
  
  # Resampling mode.
  results$replace <- replace
  
  # Training data.
  results$data <-data.frame(x)
  
  # Call.
  results$call <- sys.call()

  class(results) <- "morf"
  
  ## Output.
  return(results)
}
