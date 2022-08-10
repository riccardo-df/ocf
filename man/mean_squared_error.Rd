% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation-metrics.R
\name{mean_squared_error}
\alias{mean_squared_error}
\title{Mean Squared Error}
\usage{
mean_squared_error(y, predictions, use.true = FALSE)
}
\arguments{
\item{y}{Either the observed outcome vector or a matrix of true probabilities (\code{n.samples} x \code{M}). See the parameter \code{use.true} and the details section.}

\item{predictions}{Matrix of predictions (\code{n.samples} x \code{M}).}

\item{use.true}{If \code{FALSE} (the default), then the program assumes that \code{y} stores the observed outcome vector, otherwise it treats \code{y} as a matrix of true probabilities.}
}
\value{
The mean squared error of the method.
}
\description{
Computes the mean squared error for evaluating the accuracy of ordered probability predictions.
}
\details{
If \code{use.true = FALSE}, the mean squared error is computed as follows:

\deqn{\frac{1}{n} \sum_{i = 1}^n \sum_{m = 1}^M (1 (Y_i = m) - \hat{p}_m (x))^2}

otherwise:

\deqn{\frac{1}{n} \sum_{i = 1}^n \sum_{m = 1}^M (p_m (x) - \hat{p}_m (x))^2}

where:

\deqn{p_m (x) = P(Y_i = m \, | \, X_i = x)}

The second formula is useful for simulation studies.
}