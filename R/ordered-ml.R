#' Ordered Machine Learning
#' 
#' Estimation strategy to estimate conditional choice probabilities for ordered non-numeric outcomes.
#'
#' @param Y Outcome vector.
#' @param X Covariate matrix (no intercept).
#' @param learner String, either \code{"forest"} or \code{"l1"}. Selects the base learner to estimate each expectation.
#' @param scale Logical, whether to scale the covariates. Ignored if \code{learner} is not \code{"l1"}.
#' 
#' @return 
#' Object of class \code{oml}.
#' 
#' @examples 
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#' 
#' data <- generate_ordered_data(100)
#' sample <- data$sample
#' Y <- sample$Y
#' X <- sample[, -1]
#' 
#' ## Training-test split.
#' train_idx <- sample(seq_len(length(Y)), floor(length(Y) * 0.5))
#' 
#' Y_tr <- Y[train_idx]
#' X_tr <- X[train_idx, ]
#' 
#' Y_test <- Y[-train_idx]
#' X_test <- X[-train_idx, ]
#' 
#' ## Fit ordered machine learning on training sample using two different learners.
#' ordered_forest <- ordered_ml(Y_tr, X_tr, learner = "forest")
#' ordered_l1 <- ordered_ml(Y_tr, X_tr, learner = "l1")
#' 
#' ## Predict out of sample.
#' predictions_forest <- predict(ordered_forest, X_test)
#' predictions_l1 <- predict(ordered_l1, X_test)
#' 
#' ## Compare predictions.
#' cbind(head(predictions_forest), head(predictions_l1))}
#' 
#' @details
#' Ordered machine learning expresses conditional choice probabilities as the difference between the cumulative probabilities 
#' of two adjacent classes, which in turn can be expressed as conditional expectations of binary variables:
#' 
#' \deqn{p_m \left( X_i \right) = \mathbb{E} \left[ 1 \left( Y_i \leq m \right) | X_i \right] - \mathbb{E} \left[ 1 \left( Y_i \leq m - 1 \right) | X_i \right]}
#' 
#' Then we can separately estimate each expectation using any regression algorithm and pick the difference between the m-th and the
#' (m-1)-th estimated surfaces to estimate conditional probabilities.\cr
#' 
#' \code{\link{ordered_ml}} combines this strategy with either regression forests or penalized logistic regressions with an L1 penalty,
#' according to the user-specified parameter \code{learner}.\cr
#' 
#' If \code{learner == "forest"}, then the \code{\link[orf]{orf}}
#' function is called from an external package, as this estimator has already been proposed by Lechner and Okasa (2019).\cr
#' 
#' If \code{learner == "l1"}, 
#' the penalty parameters are chosen via 10-fold cross-validation and \code{\link[stats]{model.matrix}} is used to handle non-numeric covariates. 
#' Additionally, if \code{scale == TRUE}, the covariates are scaled to have zero mean and unit variance.
#' 
#' @import ranger glmnet orf
#' 
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2023). Ordered Correlation Forest. arXiv preprint \href{https://arxiv.org/abs/2309.08755}{arXiv:2309.08755}.
#' }
#'
#' @seealso \code{\link{multinomial_ml}}, \code{\link{ocf}}
#' 
#' @export
ordered_ml <- function(Y = NULL, X = NULL,
                           learner = "forest", scale = TRUE) {
  ## 0.) Handling inputs and checks.
  y <- Y
  check_x_y(X, y)
  if (!(learner %in% c("forest", "l1"))) stop("Invalid 'learner'. This must be either 'forest' or 'l1'.", call. = FALSE)
  y <- as.numeric(y)
  n <- length(y)
  n_categories <- length(unique(y))
  y_classes <- sort(unique(y))
  
  ## 1.) Fit the estimator and get predictions.
  if (learner == "forest") {
    estimators <- orf::orf(X, y, num.trees = 2000, honesty = FALSE)$forests
    predictions <- lapply(estimators, function(x) {predict(x, X)$predictions}) 
  } else if (learner == "l1") {
    # Generate binary outcomes for each class. The m-th element stores the indicator variables relative to the m-th class. We do not need the M-th element.
    train_outcomes <- list()
    counter <- 1
    for (m in y_classes[-n_categories]) {
      train_outcomes[[counter]] <- ifelse(y <= m, 1, 0)
      counter <- counter + 1
    }  
    
    # Construct design matrix and scale if necessary.
    X_design <- stats::model.matrix(y ~ ., data = data.frame(y, X))[, -1]
    if (scale) X_design <- scale(as.matrix(X_design))
    
    # 10-fold CV to find best penalization parameters.
    estimators <- lapply(train_outcomes, function(outcome) {glmnet::cv.glmnet(x = X_design, y = outcome, alpha = 1, family = "binomial")})
    predictions <- lapply(estimators, function(x) {as.numeric(predict(x, X_design, s = "lambda.min", type = "response"))}) 
  }
  
  ## 2.) Pick differences.
  predictions1 <- append(predictions, list(rep(1, n))) 
  predictions0 <- append(list(rep(0, n)), predictions) 
  differences <- as.list(mapply(function(x, y) x - y, predictions1, predictions0, SIMPLIFY = FALSE))
  
  ## 3.) Truncate, put into matrix, and normalize. 
  predictions_final <- lapply(differences, function(x) ifelse(x < 0, 0, x))
  predictions_final <- sapply(predictions_final, function(x) as.matrix(x))
  predictions_final <- matrix(apply(predictions_final, 1, function(x) (x) / (sum(x))), ncol = n_categories, byrow = T)
  colnames(predictions_final) <- paste0("P(Y=", seq_len(n_categories), ")")
  
  ## 4.) Output.
  output <- list("estimators" = estimators,
                 "predictions" = predictions_final,
                 "scaling" = scale,
                 "learner" = learner,
                 "X" = X,
                 "Y" = y)
  class(output) <- "oml"
  
  ## 8.) Output.
  return(output)
}
