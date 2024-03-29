% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{forest_weights_fitted_cpp}
\alias{forest_weights_fitted_cpp}
\title{Forest In-Sample Honest Weights}
\usage{
forest_weights_fitted_cpp(
  leaf_IDs_train_list,
  leaf_IDs_honest_list,
  leaf_size_honest_list
)
}
\arguments{
\item{leaf_IDs_train_list}{List of size \code{n.trees}, storing the leaf of each tree where training units fall into.}

\item{leaf_IDs_honest_list}{List of size \code{n.trees}, storing the leaf of each tree where honest units fall into.}

\item{leaf_size_honest_list}{List of size \code{n.trees}, storing the size of the leaves of each tree computed with honest units.}
}
\description{
Computes forest in-sample honest weights for a \code{ocf.forest} object relative to the m-th class.
}
\keyword{internal}
