#' Modified Ordered Random Forest
#' 
#' Nonparametric estimator of the ordered choice model using random forests. The estimator modifies a standard random forest
#' splitting criterion to build a collection of forests, each estimating the conditional probability of a single class.
#'
#' @param y Outcome vector.
#' @param X Covariate matrix (no intercept).
#' @param n.trees Number of trees.
#' @param alpha Controls the balance of each split. Each split leaves at least a fraction \code{alpha} of observations in the parent node on each side of the split.
#' @param honesty Whether to grow honest forests.
#' @param honesty.fraction Fraction of honest sample. Ignored if \code{honesty = FALSE}.
#' @param inference Whether to extract weights and compute standard errors. The weights extraction considerably slows down the program. \code{honesty = TRUE} is required for valid inference.
#' @param mtry Number of covariates to possibly split at in each node. Default is the square root of the number of covariates.
#' @param min.node.size Minimal node size.
#' @param max.depth Maximal tree depth. A value of 0 corresponds to unlimited depth, 1 to "stumps" (one split per tree).
#' @param replace If \code{TRUE}, grow trees on bootstrap subsamples. Otherwise, trees are grown on random subsamples drawn without replacement. 
#' @param sample.fraction Fraction of observations to sample. 
#' 
#' @return 
#' Object of class \code{morf}.
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
#' ## We have compatibility with generic S3-methods.
#' print(forests)
#' summary(forests)
#' predictions <- predict(forests, X_test)
#' head(predictions$probabilities)
#' table(y_test, predictions$classification)
#' 
#' ## Compute standard errors. This requires honest forests.
#' honest_forests <- morf(y_tr, X_tr, honesty = TRUE, inference = TRUE)
#' honest_forests$predictions$standard.errors}
#' 
#' @import utils stats orf
#' @importFrom Rcpp evalCpp
#' @useDynLib morf
#' 
#' @seealso \code{\link{marginal_effects}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
morf <- function(y = NULL, X = NULL,
                 honesty = FALSE, honesty.fraction = 0.5, inference = FALSE, alpha = 0,
                 n.trees = 2000, mtry = ceiling(sqrt(ncol(X))), min.node.size = 5, max.depth = 0, 
                 replace = FALSE, sample.fraction = ifelse(replace, 1, 0.5)) {
  ## 0.) Defaults for variables not needed.
  splitrule.num <- 1; treetype <- 3; probability <- FALSE
  importance.mode <- 1; scale.permutation.importance <- FALSE; local.importance <- FALSE
  respect.unordered.factors <- "ignore"; unordered.factor.variables <- c("0", "0"); use.unordered.factor.variables <- FALSE
  order.snps <- FALSE; prediction.mode <- FALSE; predict.all <- FALSE; prediction.type <- 1; oob.error <- FALSE
  loaded.forest <- list(); write.forest <- TRUE; snp.data <- as.matrix(0); class.weights <- rep(1, nlevels(y))
  alpha_balance <- alpha; alpha <- 0.5; minprop <- 0.1; num.random.splits <- 1; use.regularization.factor <- FALSE; 
  regularization.factor <- c(0, 0); regularization.usedepth <- FALSE; case.weights <- c(0, 0); use.case.weights <- FALSE; 
  split.select.weights <- list(c(0, 0)); use.split.select.weights <- FALSE; always.split.variables <- c("0", "0"); 
  use.always.split.variables <- FALSE; inbag <- list(c(0 ,0)); use.inbag <- FALSE; keep.inbag <- FALSE; holdout <- FALSE; 
  n.threads <- 0; save.memory <- FALSE; verbose <- TRUE; seed <- stats::runif(1 , 0, .Machine$integer.max)
  
  ## 1.) Handling inputs and checks.
  # 1a.) Generic checks.
  check_x_y(X, y)
  check_honesty_inference(honesty, honesty.fraction, inference)
  check_ntrees(n.trees)
  check_alpha(alpha)
  mtry <- check_mtry(mtry, colnames(X))
  check_minnodesize(min.node.size)
  check_maxdepth(max.depth)
  check_samplefraction(sample.fraction)
  
  # 1b.) Handle factor outcomes.
  if (is.factor(y)) y <- as.numeric(y)
  
  # 1c.) Store useful variables.
  y.classes <- sort(unique(y))
  n.classes <- length(y.classes)
  
  # 1d.) Recode characters as factors.
  covariate.levels <- NULL
  if (!is.matrix(X) && !inherits(X, "Matrix") && ncol(X) > 0) {
    character.idx <- sapply(X, is.character)
    X[character.idx] <- lapply(X[character.idx], factor)
    
    if (any(sapply(X, is.factor))) {
      covariate.levels <- lapply(X, levels)
    }     
  }
  
  # 1e.) Error if no covariates (X must have colnames).
  independent.variable.names <- colnames(X)
  all.independent.variable.names <- independent.variable.names
  if (length(all.independent.variable.names) < 1) stop("No covariates found. Maybe 'X' is missing? colnames?", call. = FALSE)
  
  ## 2.) Handling training data.
  # 2a.) Honest split if necessary. 
  if (honesty) {
    honest_split <- class_honest_split(data.frame(y, X), honesty.fraction)
    
    train_sample <- honest_split$train_sample
    honest_sample <- honest_split$honest_sample
    colnames(train_sample) <- c("y", independent.variable.names)
    colnames(honest_sample) <- c("y", independent.variable.names)
    
    y_train <- train_sample[, 1]
    x_train <- as.data.frame(train_sample[, -1])
    colnames(x_train) <- independent.variable.names
    
    y_honest <- honest_sample[, 1]
    x_honest <- as.data.frame(honest_sample[, -1])
    colnames(x_honest) <- independent.variable.names
  } else { 
    train_sample <- data.frame(y, X)
    honest_sample <- list()
    colnames(train_sample) <- c("y", independent.variable.names)
    
    y_train <- train_sample[, 1]
    x_train <- as.data.frame(train_sample[, -1])
    colnames(x_train) <- independent.variable.names
    
    y_honest <- NULL
    x_honest <- NULL
  }
  
  # 2b.) Transform training data into sparse matrix.
  if (inherits(x_train, "dgCMatrix")) {
    sparse.x <- x_train
    x_train <- matrix(c(0, 0))
    use.sparse.data <- TRUE
  } else {
    sparse.x <- Matrix::Matrix(matrix(c(0, 0)))
    use.sparse.data <- FALSE
    
    if (is.data.frame(x_train)) {
      x_train <- data.matrix(x_train)
    }
  }
  
  ## 3.) Generating binary outcomes for each class. The m-th element stores the indicator variables relative to the m-th class.
  train_outcomes <- list()
  honest_outcomes <- list()
  counter <- 1
  for (m in y.classes) {
    train_outcomes[[counter]] <- data.frame("y_m_train" = ifelse(y_train <= m, 1, 0), "y_m_1_train" = ifelse(y_train <= m -1, 1, 0))
    honest_outcomes[[counter]] <- data.frame("y_m_honest" = ifelse(y_honest <= m, 1, 0), "y_m_1_honest" = ifelse(y_honest <= m -1, 1, 0))
    counter <- counter + 1
  }
  
  ## 4.) Fitting a modified ordered forest to each pair of binary outcomes.
  forest_output <- lapply(train_outcomes, function(x) {
    morfCpp(treetype, x_train, matrix(c(as.numeric(x$y_m_1_train), as.numeric(x$y_m_train)), ncol = 2), 
            independent.variable.names, mtry, n.trees, verbose, seed, n.threads, write.forest, importance.mode, min.node.size, 
            split.select.weights, use.split.select.weights, always.split.variables, use.always.split.variables, 
            prediction.mode, loaded.forest, snp.data, replace, probability, unordered.factor.variables, 
            use.unordered.factor.variables, save.memory, splitrule.num, case.weights, use.case.weights, class.weights, 
            predict.all, keep.inbag, sample.fraction, alpha, minprop, holdout, prediction.type, num.random.splits, sparse.x, 
            use.sparse.data, order.snps, oob.error, max.depth, inbag, use.inbag, regularization.factor,
            use.regularization.factor, regularization.usedepth, alpha_balance)})
  if (any(sapply(forest_output, function(x) {length(x) == 0}))) stop("User interrupt or internal error.", call. = FALSE)
  
  ## 5.) Handling forests' output.
  # 5a.) Names for variable importance.
  forest_output <- lapply(forest_output, function(x) {
    names(x$variable.importance) <- all.independent.variable.names
    x})
  
  # 5b.) Write forest objects.
  if (is.factor(y_train)) {
    forest_output <- lapply(forest_output, function(x) {
      x$forest$levels <- levels(y_train)
      x})
  }
  forest_output <- lapply(forest_output, function(x) {
    x$forest$covariate.names <- independent.variable.names
    x$forest$treetype <- "modified.ordered"
    x})
  if (!is.null(covariate.levels)) {
    forest_output <- lapply(forest_output, function(x) {
      x$forest$covariate.levels <- covariate.levels
      x})
  }
  forests <- lapply(forest_output, function(x) {x$forest})
  forests <- lapply(forests, function(x) {
    class(x) <- "morf.forest"
    x})
  names(forests) <- paste("forest", y.classes, sep = ".")
  
  # 5c.) Save and remove unnecessary element.
  n.trees <- lapply(forest_output, function(x) {x$num.trees})[[1]]
  n.covariates <- lapply(forest_output, function(x) {x$num.covariates})[[1]]
  mtry <- lapply(forest_output, function(x) {x$mtry})[[1]]
  min.node.size <- lapply(forest_output, function(x) {x$min.node.size})[[1]]
  overall.importance <- rowMeans(as.matrix(sapply(forest_output, function(x) {x$variable.importance})))
  relative.importance <- round(overall.importance / sum(overall.importance), 3)
  
  ## 6) Predictions. 
  # 6a.) If inference, extract weights and use these to predict. Otherwise, use standard, faster prediction method.
  if (inference) { 
    weights <- lapply(forests, function(x) {forest_weights_fitted(x, honest_sample, train_sample)})
    class.probabilities <- mapply(function(x, y) {x %*% (y$y_m_honest - y$y_m_1_honest)}, weights, honest_outcomes)
    
    n_honest <- dim(honest_sample)[1]
    sample_correction <- n_honest / (n_honest - 1)
    products <- mapply(function(x, y) {t(apply(x, 1, function(z) {z * (y$y_m_honest - y$y_m_1_honest)}))}, weights, honest_outcomes, SIMPLIFY = FALSE)
    sums_squares <- lapply(products, function(x) {rowSums((x - rowMeans(x))^2)})
    variances <- matrix(unlist(lapply(sums_squares, function(x) {sample_correction * x}), use.names = FALSE), ncol = n.classes)
    colnames(variances) <- paste("Y=", y.classes, sep = "")
  } else if (honesty) {
    rownames(x_train) <- rownames(honest_split$train_sample)
    rownames(x_honest) <- rownames(honest_split$honest_sample)
    class.probabilities <- mapply(function(x, y) {honest_fitted(x, x_train, x_honest, y$y_m_honest, y$y_m_1_honest)}, forests, honest_outcomes)
    variances <- list()
  } else {
    class.probabilities <- matrix(unlist(lapply(forests, function(x) {predict(x, data = x_train)$predictions}), use.names = FALSE), ncol = n.classes)
    variances <- list()
  } 
  
  # 6b.) Normalization step.
  class.probabilities <- matrix(apply(class.probabilities, 1, function(x) (x)/(sum(x))), ncol = n.classes, byrow = TRUE)
  colnames(class.probabilities) <- paste("P(Y=", y.classes, ")", sep = "")
  
  # 6c.) Classification.
  classification <- apply(class.probabilities, 1, which.max)
  
  ## 7.) Construct morf object.
  output <- list()
  output$forests.info <- forests
  output$predictions <- list("probabilities" = class.probabilities, "standard.errors" = if (inference) sqrt(variances) else list(),
                             "classification" = classification)
  output$importance <- relative.importance
  output$tuning.info <- list("n.trees" = n.trees, "mtry" = mtry, "min.node.size" = min.node.size, "replace" = replace,
                             "honesty" = honesty, "honesty.fraction" = if (honesty) honesty.fraction else 0, "call" = sys.call())
  output$full_data <- data.frame(y, X)
  colnames(output$full_data) <- colnames(train_sample)
  output$honest_data <- if (honesty) data.frame(y_honest, x_honest) else list()
  if (honesty) colnames(output$honest_data) <- colnames(train_sample)

  class(output) <- "morf"
  
  ## 8.) Output.
  return(output)
}
