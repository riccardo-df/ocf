##' Check Arguments x and y (Internal Use)
##' 
##' @param x Covariate matrix (no intercept).
##' @param y Outcome vector.
check_x_y <- function(x, y) {
  # Null.
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
    if (!inherits(x, "dgCMatrix")) stop("Currently only sparse data of class 'dgCMatrix' supported.", call. = FALSE)
  }
}


##' Check Arguments honesty, honesty.fraction and inference (Internal Use)
##' 
##' @param honesty Whether to grow honest forests.
##' @param honesty.fraction Fraction of honest sample.
##' @param inference Whether to conduct weight-based inference.
check_honesty_inference <- function(honesty, honesty.fraction, inference) {
  if (!(honesty %in% c(TRUE, FALSE))) stop("Invalid value of 'honesty'.", call. = FALSE)
  if (!(inference %in% c(TRUE, FALSE))) stop("Invalid value of 'inference'.", call. = FALSE)
  
  if (honesty.fraction < 0 | honesty.fraction >= 1) stop("Invalid value of 'honesty.fraction'. This must be a value in [0, 1)", call. = FALSE)
  if(!honesty & inference) stop("Valid inference requires honest forests.", call. = FALSE)
}


##' Check Argument n.trees (Internal Use)
##' 
##' @param n.trees Number of trees.
check_ntrees <- function(n.trees) {
  if (!is.numeric(n.trees) || n.trees < 1) stop("Invalid value for 'n.trees'.", call. = FALSE)
}


##' Check Argument mtry (Internal Use)
##' 
##' @param mtry Number of covariates to possibly split at in each node. Default is the (rounded down) square root of the number of covariates. Alternatively, one can pass a single-argument function returning an integer, where the argument is the number of covariates.
##' @param nv Number of covariates.
##' @return 
##' Appropriate value of \code{mtry}.
check_mtry <- function(mtry, nv) {
  if (is.function(mtry)) {
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
  
  return(mtry)
}


##' Check Argument seed (Internal Use)
##' 
##' @param seed Random seed. 
##' 
##' @return 
##' Appropriate value of \code{seed}.
check_seed <- function(seed) {
  if (is.null(seed)) seed <- stats::runif(1 , 0, .Machine$integer.max)
  
  return(seed)
}


##' Check Argument keep.inbag (Internal Use)
##' 
##' @param keep.inbag Save how often observations are in-bag in each tree. 
check_keepinbag <- function(keep.inbag) {
  if (!is.logical(keep.inbag)) stop("Invalid value for 'keep.inbag'", call. = FALSE)
}


##' Check Argument n.threads (Internal Use)
##' 
##' @param n.threads Number of threads. Default is number of CPUs available.
##' 
##' @return 
##' Appropriate value of \code{n.threads}.
check_nthreads <- function(n.threads) {
  if (is.null(n.threads)) {
    n.threads <- 0
  } else if (!is.numeric(n.threads) || n.threads < 0) {
    stop("Invalid value for 'n.threads'", call. = FALSE)
  }
  
  return(n.threads)
}


##' Check Argument min.node.size (Internal Use)
##' 
##' @param min.node.size Minimal node size.
check_minnodesize <- function(min.node.size) {
  if (!is.numeric(min.node.size) || min.node.size <= 0) stop("Invalid value for 'min.node.size'", call. = FALSE)
}


##' Check Argument max.depth (Internal Use)
##' 
##' @param max.depth Maximal tree depth. A value of 0 corresponds to unlimited depth, 1 to "stumps" (one split per tree).
check_maxdepth <- function(max.depth) {
  if (!is.numeric(max.depth) || max.depth < 0) stop("Invalid value for 'max.depth'. Please give a positive integer.", call. = FALSE)
}


##' Check Argument sample.fraction (Internal Use)
##' 
##' @param sample.fraction Fraction of observations to sample. 
check_samplefraction <- function(sample.fraction) {
  if (!is.numeric(sample.fraction)) stop("Invalid value for 'sample.fraction'. Please give a value in (0,1].", call. = FALSE)
  if (sample.fraction <= 0 || sample.fraction > 1) stop("Invalid value for 'sample.fraction' Please give a value in (0,1].", call. = FALSE)
}


##' Check Argument alpha (Internal Use)
##' 
##' @param alpha Fraction of observations that must lie on each side of each split. 
check_alpha <- function(alpha) {
  if (!is.numeric(alpha)) stop("Invalid value for 'alpha'. Please give a value in [0, 0.5].", call. = FALSE)
  if (alpha < 0 | alpha > 0.5) stop("Invalid value for 'alpha'. Please give a value in [0, 0.5].", call. = FALSE)
}
