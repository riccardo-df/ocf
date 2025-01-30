#' Multinomial Machine Learning
#' 
#' Estimation strategy to estimate conditional choice probabilities for ordered non-numeric outcomes.
#'
#' @param Y Outcome vector.
#' @param X Covariate matrix (no intercept).
#' @param learner String, either \code{"forest"} or \code{"l1"}. Selects the base learner to estimate each expectation.
#' @param scale Logical, whether to scale the covariates. Ignored if \code{learner} is not \code{"l1"}.
#' 
#' @return 
#' Object of class \code{mml}.
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
#' ## Fit multinomial machine learning on training sample using two different learners.
#' multinomial_forest <- multinomial_ml(Y_tr, X_tr, learner = "forest")
#' multinomial_l1 <- multinomial_ml(Y_tr, X_tr, learner = "l1")
#' 
#' ## Predict out of sample.
#' predictions_forest <- predict(multinomial_forest, X_test)
#' predictions_l1 <- predict(multinomial_l1, X_test)
#' 
#' ## Compare predictions.
#' cbind(head(predictions_forest), head(predictions_l1))}
#' 
#' @details
#' Multinomial machine learning expresses conditional choice probabilities as expectations of binary variables:
#' 
#' \deqn{p_m \left( X_i \right) = \mathbb{E} \left[ 1 \left( Y_i = m \right) | X_i \right]}
#' 
#' This allows us to estimate each expectation separately using any regression algorithm to get an estimate of conditional probabilities.\cr
#' 
#' \code{\link{multinomial_ml}} combines this strategy with either regression forests or penalized logistic regressions with an L1 penalty,
#' according to the user-specified parameter \code{learner}.\cr
#' 
#' If \code{learner == "l1"}, the penalty parameters are chosen via 10-fold cross-validation 
#' and \code{\link[stats]{model.matrix}} is used to handle non-numeric covariates. Additionally, if \code{scale == TRUE}, the covariates are scaled to 
#' have zero mean and unit variance.
#' 
#' @import ranger glmnet
#' 
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2025). Ordered Correlation Forest. Econometric Reviews, 1–17. \href{https://doi.org/10.1080/07474938.2024.2429596}{https://doi.org/10.1080/07474938.2024.2429596}.
#' }
#' 
#' @seealso \code{\link{ordered_ml}}, \code{\link{ocf}}
#' 
#' @export
multinomial_ml <- function(Y = NULL, X = NULL,
                           learner = "forest", scale = TRUE) {
  ## 0.) Handling inputs and checks.
  y <- Y
  check_x_y(X, y)
  if (!(learner %in% c("forest", "l1"))) stop("Invalid 'learner'. This must be either 'forest' or 'l1'.", call. = FALSE)
  if (!is.logical(scale)) stop("Invalid 'scale'. This must be logical.", call. = FALSE)
  
  n_categories <- length(unique(y))
  y_classes <- sort(unique(y))
  
  ## 1.) Generate binary outcomes for each class. The m-th element stores the indicator variable relative to the m-th class. 
  train_outcomes <- list()
  counter <- 1
  for (m in y_classes) {
    train_outcomes[[counter]] <- ifelse(y == m, 1, 0)
    counter <- counter + 1
  }  

  ## 2.) Fit the estimator on each dummy and get predictions.
  if (learner == "forest") {
    estimators <- lapply(train_outcomes, function(outcome) {ranger::ranger(y = outcome, x = X, num.trees = 2000)})
    predictions <- lapply(estimators, function(x) {predict(x, X)$predictions}) 
  } else if (learner == "l1") {
    # Construct design matrix and scale if necessary.
    X_design <- stats::model.matrix(y ~ ., data = data.frame(y, X))[, -1]
    if (scale) X_design <- scale(as.matrix(X_design))
    
    # 10-fold CV to find best penalization parameters.
    estimators <- lapply(train_outcomes, function(outcome) {glmnet::cv.glmnet(x = X_design, y = outcome, alpha = 1, family = "binomial")})
    predictions <- lapply(estimators, function(x) {as.numeric(predict(x, X_design, s = "lambda.min", type = "response"))}) 
  }
  
  ## 3.) Put into matrix and normalize.
  predictions_final <- sapply(predictions, function(x) as.matrix(x))
  predictions_final <- matrix(apply(predictions_final, 1, function(x) (x) / (sum(x))), ncol = n_categories, byrow = T)
  colnames(predictions_final) <- paste0("P(Y=", seq_len(n_categories), ")")
  
  ## 4.) Output.
  output <- list("estimators" = estimators,
                 "predictions" = predictions_final,
                 "learner" = learner,
                 "scaling" = scale,
                 "X" = X,
                 "Y" = y)
  class(output) <- "mml"
  
  ## 8.) Output.
  return(output)
}
