##' Print Method for Morf xs
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
  cat("Number of trees:                 ", x$n.trees, "\n")
  cat("Training sample size:            ", x$n.samples, "\n")
  cat("Number of covariates:            ", x$n.covariates, "\n")
  cat("Mtry:                            ", x$mtry, "\n")
  cat("Minimum node size:               ", x$min.node.size, "\n")
  cat("Honesty:                         ", x$honesty, "\n")
  cat("Fraction honesty:                ", x$honesty.fraction, "\n")
  cat("MSE:                             ", x$mean.squared.error, "\n")
  cat("RPS:                             ", x$mean.ranked.score)
}


##' Print Method for Morf Marginal Effects
##'
##' Prints a \code{morf.marginal} object.
##'
##' @param x \code{morf.marginal} object.
##' @param latex If \code{TRUE}, prints a latex code for a table displaying the marginal effects.
##' @param ... Further arguments passed to or from other methods.
##' 
##' @seealso \code{\link{morf}} and \code{\link{marginal_effects}}.
##' 
##' @author Riccardo Di Francesco
##' 
##' @export
print.morf.marginal <- function(x, latex = FALSE, ...) {
  if (!(latex %in% c(TRUE, FALSE))) stop("Invalid value of 'latex'.", call. = FALSE)
  
  cat("Morf marginal effects results \n\n")
  cat("Evaluation:                      ", x$evaluation, "\n")
  cat("Bandwidth:                       ", x$bandwitdh, "\n")
  cat("Number of classes:               ", x$n.classes, "\n")
  cat("Number of trees:                 ", x$n.trees, "\n")
  cat("Sample size:                     ", x$n.samples, "\n")
  cat("Honest forest:                   ", x$honesty, "\n")
  
  cat("\n\n")
  
  cat("Marginal Effects: \n\n")
  
  print(x$marginal.effects)
}
