#' Modified Ordered Random Forest
#' 
#' Nonparametric estimation of the ordered choice model using random forests.
#'
#' @param x Covariate matrix (no intercept).
#' @param y Outcome vector.
#' @param n.trees Number of trees.
#' @param alpha Controls the balance of each split. Each split leaves at least a fraction \code{alpha} of observations in the parent node on each side of the split.
#' @param honesty Whether to grow honest forests.
#' @param honesty.fraction Fraction of honest sample. Ignored if \code{honesty = FALSE}.
#' @param inference Whether to conduct weight-based inference. The weights' extraction considerably slows down the program. \code{honesty = TRUE} is required for valid inference.
#' @param mtry Number of covariates to possibly split at in each node. Default is the square root of the number of covariates.
#' @param min.node.size Minimal node size.
#' @param max.depth Maximal tree depth. A value of 0 corresponds to unlimited depth, 1 to "stumps" (one split per tree).
#' @param replace If \code{TRUE}, grow trees on bootstrap subsamples. Otherwise, trees are grown on random subsamples drawn without replacement. 
#' @param sample.fraction Fraction of observations to sample. 
#' @param case.weights Weights for sampling training observations. Observations with larger weights will be drawn with higher probability in the subsamples for the trees.
#' @param split.select.weights Numeric vector with weights between 0 and 1, used to calculate the probability to select variables for splitting. Alternatively, one can use a list of size \code{n.trees} containing \code{split.select.weights} vectors, one for each tree.  
#' @param always.split.variables Character vector with variable names to be always selected in addition to the \code{mtry} variables tried for splitting.
#' @param keep.inbag Save how often observations are in-bag in each tree. 
#' @param inbag Manually set observations per tree. List of size \code{n.trees}, containing in-bag counts for each observation. Can be used for stratified sampling.
#' @param holdout Hold-out mode. Hold-out all samples with zero \code{case.weights} and use these for variable importance and prediction error.
#' @param n.threads Number of threads. Default is number of CPUs available.
#' @param save.memory Use memory saving splitting mode. It slows down the tree growing, use only if you encounter memory problems.
#' @param verbose Show computation status and estimated runtime.
#' @param seed Random seed. Default is \code{NULL}, which generates the seed from \code{R}. Set to \code{0} to ignore the \code{R} seed. 
#' 
#' @details 
#' \subsection{Splitting Criterion}{
#' \code{morf} fits \code{M} separated random forests, where \code{M} is the number of classes of \code{y}. 
#' Each forest computes the \code{m}-th class conditional choice probabilities:
#' 
#' \deqn{p_m (x) := P ( Y_i = m \, | \, X_i = x)}
#' 
#' To estimate this quantity, \code{morf} exploits the following:
#' 
#' \deqn{p_m (x) = E [\, 1 (Y_i \leq m) \, | \, X_i = x] - E[\, 1 ( Y_i \leq m - 1 ) \, | \, X_i = x]  
#'                 = \mu_m (x) - \mu_{m-1} (x)}
#' 
#' with \code{1(.)} an indicator of the truth of its argument. A straightforward estimator consists of
#' estimating the two conditional expectations separately and taking the difference:
#' 
#' \deqn{\hat{p}_m (x) = \hat{\mu}_m (x) - \hat{\mu}_{m-1} (x)}
#' 
#' However, this strategy ignores potential correlation in the estimation errors of the two surfaces. An alternative 
#' approach is to tackle the minimization of the mean squared error (MSE) of the particular estimation problem:
#' 
#' \deqn{MSE[\hat{p}_m (x)] = MSE[\hat{\mu}_m (x)] + MSE [\hat{\mu}_{m-1} (x)] - 2 MCE[\hat{\mu}_m (x), \hat{\mu}_{m-1} (x)]}
#' 
#' where the last term is the mean correlation error of the estimation. \code{morf} grows forests that tie the estimation 
#' of the two conditional expectations together to make the correlation term positive. For this purpose, trees are built 
#' by greedily minimizing the following expression:
#' 
#' \deqn{Var( 1 (Y_i \leq m)) + Var( 1 (Y_i \leq m - 1)) - 2 * Cor(1 (Y_i \leq m), 1 (Y_i \leq m - 1))}
#' }
#' 
#' \subsection{Predictions}{
#' Predictions in the \code{l}-th leaf are computed as:
#' 
#' \deqn{\frac{1}{\{i: x_i \in L_l\}} \sum_{\{i: x_i \in L_l\}} 1 (Y_i \leq m) - \frac{1}{\{i: x_i \in L_l\}} \sum_{\{i: x_i \in L_l\}} 1 (Y_i \leq m - 1)}
#' 
#' Notice that a normalization step may be needed to ensure that the estimated probabilities sum up to one across classes.
#' }
#' 
#' \subsection{Variable Importance}{
#' For each covariate, an overall variable importance measure is provided. The \code{m}-th forest computes the 
#' importance of the j-th covariate for the \code{m}-th class by recording the improvement in the splitting criterion 
#' at each split placed on such covariate. Summing over all such splits of the trees in the \code{m}-th 
#' forest gives the j-th covariate's importance for the \code{m}-th class. The overall variable importance measure of this 
#' covariate is then defined as the mean of its importances in each class.
#' }
#' 
#' \subsection{Honest Forests}{
#' Growing honest forests is a necessary requirements to conduct valid inference. \code{morf} implements honest estimation 
#' as follows. The data set is split into a training sample and a honest sample. Forests are grown using only the training 
#' sample. Then, for each prediction point, the weights relative to units from the honest sample are extracted, and 
#' predictions are based on the honest outcomes (relying on the prediction method outlined above). This way, weights and
#' outcomes are independent of each other, thereby allowing for valid weight-based inference.
#' }
#' 
#' @return 
#' Object of class \code{morf} with elements:
#'   \item{\code{forest.1}}{\code{morf.forest} object of the first class.}
#'   \item{\code{forest.2}}{\code{morf.forest} object of the second class.} 
#'   \item{\code{...}}{}
#'   \item{\code{forest.M}}{\code{morf.forest} object of the last class.}
#'   \item{\code{predictions}}{Matrix of predicted conditional class probabilities. If \code{honesty = TRUE}, these are honest predictions.}
#'   \item{\code{standard.errors}}{Standard error of the predictions. Requires \code{inference = TRUE}.}
#'   \item{\code{mean.squared.error}}{Mean squared error of the model, based on \code{predictions}.}
#'   \item{\code{mean.ranked.score}}{Mean ranked probability score of the model, based on \code{predictions}.}
#'   \item{\code{overall.importance}}{Measure of overall variable importance, computed as the mean of the importance of each variable across classes. Relative importance is provided.}
#'   \item{\code{classes}}{Possible outcome values.}
#'   \item{\code{n.classes}}{Number of classes.}
#'   \item{\code{n.samples}}{Number of observations.}
#'   \item{\code{n.covariates}}{Number of covariates.}
#'   \item{\code{n.trees}}{Number of trees of each forest.}
#'   \item{\code{mtry}}{Number of covariates considered for splitting at each step.}
#'   \item{\code{min.node.size}}{Minimum node size.}
#'   \item{\code{replace}}{Whether the subsamples to grow trees are drawn with replacement.}
#'   \item{\code{honesty}}{Whether forests are honest.}
#'   \item{\code{honesty.fraction}}{Fraction of units allocated to honest sample.}
#'   \item{\code{full_data}}{Whole sample.}
#'   \item{\code{honest_data}}{Honest sample.}
#'   \item{\code{call}}{System call.}
#' 
#' @import utils stats
#' @importFrom Rcpp evalCpp
#' @useDynLib morf
#' 
#' @seealso \code{\link{predict.morf}}, \code{\link{marginal_effects}}
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
morf <- function(x = NULL, y = NULL,
                 n.trees = 2000, mtry = ceiling(sqrt(ncol(x))), min.node.size = 5, max.depth = 0, alpha = 0.2,
                 replace = FALSE, sample.fraction = ifelse(replace, 1, 0.5), case.weights = NULL,
                 honesty = TRUE, honesty.fraction = 0.5, inference = FALSE,
                 split.select.weights = NULL, always.split.variables = NULL,
                 keep.inbag = FALSE, inbag = NULL, holdout = FALSE,
                 n.threads = NULL, save.memory = FALSE,
                 verbose = TRUE, seed = NULL) {
  ## Handling inputs and checks.
  # Generic checks.
  check_x_y(x, y)
  check_honesty_inference(honesty, honesty.fraction, inference)
  check_ntrees(n.trees)
  check_alpha(alpha)
  alpha_balance <- alpha
  mtry <- check_mtry(mtry, colnames(x))
  seed <- check_seed(seed)
  check_keepinbag(keep.inbag)
  n.threads <- check_nthreads(n.threads)
  check_minnodesize(min.node.size)
  check_maxdepth(max.depth)
  check_samplefraction(sample.fraction)
  
  # Store useful variables.
  y.classes <- sort(unique(y))
  n.classes <- length(y.classes)
  
  # Recode characters as factors.
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
  
  ## Handling training data.
  # Honest split. 
  if (honesty) {
    honest_split <- class_honest_split(data.frame(y, x), honesty.fraction)
    
    train_sample <- honest_split$train_sample
    honest_sample <- honest_split$honest_sample
    
    y_train <- train_sample[, 1]
    x_train <- as.data.frame(train_sample[, -1])
    colnames(x_train) <- independent.variable.names
    
    y_honest <- honest_sample[, 1]
    x_honest <- as.data.frame(honest_sample[, -1])
    colnames(x_honest) <- independent.variable.names
  } else { 
    train_sample <- data.frame(y, x)
    honest_sample <- list()
    
    y_train <- train_sample[, 1]
    x_train <- as.data.frame(train_sample[, -1])
    colnames(x_train) <- independent.variable.names
    
    y_honest <- NULL
    x_honest <- NULL
  }
  
  # Use sparse matrix.
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
  
  ## Defaults for variables not needed.
  splitrule.num <- 1; treetype <- 3; probability <- FALSE
  importance.mode <- 1; scale.permutation.importance <- FALSE; local.importance <- FALSE
  respect.unordered.factors <- "ignore"; unordered.factor.variables <- c("0", "0"); use.unordered.factor.variables <- FALSE
  order.snps <- FALSE; prediction.mode <- FALSE; predict.all <- FALSE; prediction.type <- 1; oob.error <- FALSE
  loaded.forest <- list(); write.forest <- TRUE
  snp.data <- as.matrix(0)
  class.weights <- rep(1, nlevels(y))
  alpha <- 0.5; minprop <- 0.1; num.random.splits <- 1
  use.regularization.factor <- FALSE; regularization.factor <- c(0, 0); regularization.usedepth <- FALSE
  
  ## Generating outcomes for each class. The m-th element stores the indicator variables relative to the m-th class.
  train_outcomes <- list()
  honest_outcomes <- list()
  counter <- 1
  for (m in y.classes) {
    train_outcomes[[counter]] <- data.frame("y_m_train" = ifelse(y_train <= m, 1, 0), "y_m_1_train" = ifelse(y_train <= m -1, 1, 0))
    honest_outcomes[[counter]] <- data.frame("y_m_honest" = ifelse(y_honest <= m, 1, 0), "y_m_1_honest" = ifelse(y_honest <= m -1, 1, 0))
    counter <- counter + 1
  }
  
  ## Fitting a modified ordered forest to each data frame.
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
  
  ## Handling forests' output.
  # Names for variable importance.
  forest_output <- lapply(forest_output, function(x) {
    names(x$variable.importance) <- all.independent.variable.names
    x})
  
  # Write forest objects.
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
  
  # Save and remove unnecessary element.
  n.trees <- lapply(forest_output, function(x) {x$num.trees})[[1]]
  n.covariates <- lapply(forest_output, function(x) {x$num.covariates})[[1]]
  mtry <- lapply(forest_output, function(x) {x$mtry})[[1]]
  min.node.size <- lapply(forest_output, function(x) {x$min.node.size})[[1]]
  overall.importance <- rowMeans(as.matrix(sapply(forest_output, function(x) {x$variable.importance})))
  
  ## Predictions. If inference, extract weights and use these to predict. Otherwise, use standard, faster prediction method.
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
    class.probabilities <- mapply(function(x, y) {honest_fitted(x, train_sample, honest_sample, y$y_m_honest, y$y_m_1_honest)}, forests, honest_outcomes)
    variances <- list()
  } else {
    class.probabilities <- matrix(unlist(lapply(forests, function(x) {predict(x, data = x_train)$predictions}), use.names = FALSE), ncol = n.classes)
    variances <- list()
  } 
  
  # Normalization step.
  class.probabilities <- matrix(apply(class.probabilities, 1, function(x) (x)/(sum(x))), ncol = n.classes, byrow = TRUE)
  colnames(class.probabilities) <- paste("P(Y=", y.classes, ")", sep = "")
  
  ## Constructing morf object.
  output <- forests
  output$predictions <- class.probabilities
  output$standard.errors <- if (inference) sqrt(variances) else list()
  output$mean.squared.error <- mean_squared_error(y, class.probabilities)
  output$mean.ranked.score <- mean_ranked_score(y, class.probabilities)
  output$overall.importance <- overall.importance
  output$classes <- y.classes
  output$n.classes <- n.classes
  output$n.samples <- length(y)
  output$n.covariates <- n.covariates
  output$n.trees <- n.trees
  output$mtry <- mtry
  output$min.node.size <- min.node.size
  output$replace <- replace
  output$honesty <- honesty
  output$honesty.fraction <- if (honesty) honesty.fraction else 0
  output$full_data <- data.frame(y, x)
  output$honest_data <- if (honesty) data.frame(y_honest, x_honest) else list()
  output$call <- sys.call()
  
  class(output) <- "morf"
  
  ## Output.
  return(output)
}
