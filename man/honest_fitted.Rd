% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/honesty.R
\name{honest_fitted}
\alias{honest_fitted}
\title{Honest In-Sample Predictions}
\usage{
honest_fitted(forest, train_sample, honest_sample, y_m_honest, y_m_1_honest)
}
\arguments{
\item{forest}{An \code{ocf.forest} object.}

\item{train_sample}{Training sample.}

\item{honest_sample}{Honest sample.}

\item{y_m_honest}{Indicator variable, whether the outcome is smaller than or equal to the m-th class.}

\item{y_m_1_honest}{Indicator variable, whether the outcome is smaller than or equal to the (m-1)-th class.}
}
\value{
In-sample honest predictions.
}
\description{
Computes honest in-sample predictions for an \code{ocf.forest} object.
}
\details{
\code{forest} must have been grown using only the training sample. \code{honest_fitted} replaces the leaf estimates 
using the outcome from the honest sample (using the prediction method of \code{\link{ocf}}).
}
\keyword{internal}
