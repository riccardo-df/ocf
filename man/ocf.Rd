% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ocf.R
\name{ocf}
\alias{ocf}
\title{Ordered Correlation Forest}
\usage{
ocf(
  y = NULL,
  X = NULL,
  honesty = FALSE,
  honesty.fraction = 0.5,
  inference = FALSE,
  alpha = 0,
  n.trees = 2000,
  mtry = ceiling(sqrt(ncol(X))),
  min.node.size = 5,
  max.depth = 0,
  replace = FALSE,
  sample.fraction = ifelse(replace, 1, 0.5),
  n.threads = 1
)
}
\arguments{
\item{y}{Outcome vector.}

\item{X}{Covariate matrix (no intercept).}

\item{honesty}{Whether to grow honest forests.}

\item{honesty.fraction}{Fraction of honest sample. Ignored if \code{honesty = FALSE}.}

\item{inference}{Whether to extract weights and compute standard errors. The weights extraction considerably slows down the routine. \code{honesty = TRUE} is required for valid inference.}

\item{alpha}{Controls the balance of each split. Each split leaves at least a fraction \code{alpha} of observations in the parent node on each side of the split.}

\item{n.trees}{Number of trees.}

\item{mtry}{Number of covariates to possibly split at in each node. Default is the square root of the number of covariates.}

\item{min.node.size}{Minimal node size.}

\item{max.depth}{Maximal tree depth. A value of 0 corresponds to unlimited depth, 1 to "stumps" (one split per tree).}

\item{replace}{If \code{TRUE}, grow trees on bootstrap subsamples. Otherwise, trees are grown on random subsamples drawn without replacement.}

\item{sample.fraction}{Fraction of observations to sample.}

\item{n.threads}{Number of threads. Zero corresponds to the number of CPUs available.}
}
\value{
Object of class \code{ocf}.
}
\description{
Nonparametric estimator for ordered non-numeric outcomes. The estimator modifies a standard random forest
splitting criterion to build a collection of forests, each estimating the conditional probability of a single class.
}
\examples{
## Load data from orf package.
set.seed(1986)

library(orf)
data(odata)
odata <- odata[1:200, ] # Subset to reduce elapsed time.

y <- as.numeric(odata[, 1])
X <- as.matrix(odata[, -1])

## Training-test split.
train_idx <- sample(seq_len(length(y)), floor(length(y) * 0.5))

y_tr <- y[train_idx]
X_tr <- X[train_idx, ]

y_test <- y[-train_idx]
X_test <- X[-train_idx, ]

## Fit ocf on training sample.
forests <- ocf(y_tr, X_tr)

## We have compatibility with generic S3-methods.
print(forests)
summary(forests)
predictions <- predict(forests, X_test)
head(predictions$probabilities)
table(y_test, predictions$classification)
\donttest{
## Compute standard errors. This requires honest forests.
honest_forests <- ocf(y_tr, X_tr, honesty = TRUE, inference = TRUE)
head(honest_forests$predictions$standard.errors)}

}
\seealso{
\code{\link{marginal_effects}}
}
\author{
Riccardo Di Francesco
}