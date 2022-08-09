##' Print Method for Morf Objects
##'
##' Prints a \code{morf} object.
##'
##' @param x \code{morf} object.
##' @param ... Further arguments passed to or from other methods.
##' 
##' @seealso \code{\link{morf}}
##' 
##' @author Riccardo Di Francesco
##' 
##' @export
print.morf <- function(x, ...) {
  cat("Morf results \n\n")
  cat("Call:\n", deparse(x$call), "\n\n")
  cat("Number of classes:               ", x$n.classes, "\n")
  cat("Number of trees:                 ", x$num.trees, "\n")
  cat("Training sample size:            ", x$num.samples, "\n")
  cat("Number of covariates:            ", x$num.covariates, "\n")
  cat("Mtry:                            ", x$mtry, "\n")
  cat("Minimum node size:               ", x$min.node.size, "\n")
  cat("MSE:                             ", x$mean.squared.error, "\n")
}
