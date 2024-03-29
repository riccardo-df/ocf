% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weights.R
\name{predict_forest_weights}
\alias{predict_forest_weights}
\title{Forest Out-of-Sample Weights}
\usage{
predict_forest_weights(forest, honest_sample, test_sample)
}
\arguments{
\item{forest}{An \code{ocf.forest} object.}

\item{honest_sample}{Honest sample.}

\item{test_sample}{Test sample.}
}
\value{
Matrix of out-of-sample honest weights.
}
\description{
Computes forest out-of-sample honest weights for an \code{ocf.forest} object.
}
\details{
\code{forest} must have been grown using only the training sample.
}
\keyword{internal}
