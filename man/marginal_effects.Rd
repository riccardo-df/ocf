% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal-effects.R
\name{marginal_effects}
\alias{marginal_effects}
\title{Marginal Effects for Ordered Correlation Forest}
\usage{
marginal_effects(
  object,
  data = NULL,
  these_covariates = NULL,
  eval = "atmean",
  bandwitdh = 0.1,
  inference = FALSE
)
}
\arguments{
\item{object}{An \code{\link{ocf}} object.}

\item{data}{Data set of class \code{data.frame} to estimate marginal effects. It must contain at least the same covariates used to train the forests. If \code{NULL}, marginal effects are estimated on \code{object$full_data}.}

\item{these_covariates}{Named list with covariates' names as keys and strings denoting covariates' types as entries. Strings must be either \code{"continuous"} or \code{"discrete"}. The names of the list indicate the covariates for which marginal effect estimation is desired. If \code{NULL} (the default), marginal effects are estimated for all covariates and covariates' types are inferred by the routine.}

\item{eval}{Evaluation point for marginal effects. Either \code{"mean"}, \code{"atmean"} or \code{"atmedian"}.}

\item{bandwitdh}{How many standard deviations \code{x_up} and \code{x_down} differ from \code{x}.}

\item{inference}{Whether to extract weights and compute standard errors. The weights extraction considerably slows down the program.}
}
\value{
Object of class \code{ocf.marginal}.
}
\description{
Nonparametric estimation of marginal effects using an \code{\link{ocf}} object.
}
\details{
\code{\link{marginal_effects}} can estimate mean marginal effects, marginal effects at the mean, or marginal effects at the
median, according to the \code{eval} argument.\cr 

If \code{these_covariates} is \code{NULL} (the default), the routine assumes that covariates with with at most ten unique values are categorical and treats the remaining covariates as continuous.\cr
}
\examples{
\donttest{## Generate synthetic data.
set.seed(1986)

data <- generate_ordered_data(100)
sample <- data$sample
Y <- sample$Y
X <- sample[, -1]

## Fit ocf.
forests <- ocf(Y, X)

## Marginal effects at the mean.
me <- marginal_effects(forests, eval = "atmean")

print(me)
print(me, latex = TRUE)
plot(me)

## Compute standard errors. This requires honest forests.
honest_forests <- ocf(Y, X, honesty = TRUE)
honest_me <- marginal_effects(honest_forests, eval = "atmean", inference = TRUE)

print(honest_me, latex = TRUE)
plot(honest_me)

## Subset covariates and select covariates' types.
my_covariates <- list("x1" = "continuous", "x2" = "discrete", "x4" = "discrete")
honest_me <- marginal_effects(honest_forests, eval = "atmean", inference = TRUE,
                              these_covariates = my_covariates)
print(honest_me)
plot(honest_me)}

}
\references{
\itemize{
  \item Di Francesco, R. (2025). Ordered Correlation Forest. Econometric Reviews, 1–17. \doi{10.1080/07474938.2024.2429596}.
}
}
\seealso{
\code{\link{ocf}}
}
\author{
Riccardo Di Francesco
}
