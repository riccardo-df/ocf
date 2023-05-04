#' Multinomial Machine Learning
#' 
#' Estimation strategy to estimate conditional choice probabilities for discrete outcomes.
#'
#' @param y Outcome vector.
#' @param X Covariate matrix (no intercept).
#' @param learner String, either \code{"forest"} or \code{"l1"}. Selects the base learner to estimate each expectation.
#' 
#' @return 
#' Object of class \code{mml}.
#' 
#' @examples 
#' ## Load data from orf package.
#' set.seed(1986)
#' 
#' library(orf)
#' data(odata)
#' odata <- odata[1:200, ] # Subset to reduce elapsed time.
#' 
#' y <- as.numeric(odata[, 1])
#' X <- as.matrix(odata[, -1])
#' 
#' ## Training-test split.
#' train_idx <- sample(seq_len(length(y)), floor(length(y) * 0.5))
#' 
#' y_tr <- y[train_idx]
#' X_tr <- X[train_idx, ]
#' 
#' y_test <- y[-train_idx]
#' X_test <- X[-train_idx, ]
#' 
#' ## Fit multinomial machine learning on training sample using two different learners.
#' multinomial_forest <- multinomial_ml(y_tr, X_tr, learner = "forest")
#' multinomial_l1 <- multinomial_ml(y_tr, X_tr, learner = "l1")
#' 
#' ## Predict out of sample.
#' predictions_forest <- predict(multinomial_forest, X_test)
#' predictions_l1 <- predict(multinomial_l1, X_test)
#' 
#' ## Compare predictions.
#' cbind(head(predictions_forest), head(predictions_l1))
#' 
#' @details
#' Multinomial machine learning expresses conditional choice probabilities as expectations of binary variables:
#' 
#' \deqn{p_m \left( X_i \right) = \mathbb{E} \left[ 1 \left( Y_i = m \right) | X_i \right]}
#' 
#' This allows us to estimate each expectation separately using any regression algorithm to get an estimate of conditional 
#' probabilities.\cr
#' 
#' \code{\link{multinomial_ml}} combines this strategy with either regression forests or penalized logistic regression with an L1 penalty,
#' according to the user-specified parameter \code{learner}. If \code{learner == "l1"}, the covariates are scaled to have zero mean
#' and unit variance, and the penalty parameters are chosen via 10-fold cross-validation.\cr
#' 
#' @import ranger glmnet orf
#'  
#' @seealso \code{\link{ordered_ml}}, \code{\link{ocf}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
multinomial_ml <- function(y = NULL, X = NULL,
                           learner = "forest") {
  ## 0.) Handling inputs and checks.
  check_x_y(X, y)
  if (!(learner %in% c("forest", "l1"))) stop("Invalid 'learner'. This must be either 'forest' or 'l1'.", call. = FALSE)
  n_categories <- length(unique(y))
  
  ## 1.) Generate binary outcomes for each class. The m-th element stores the indicator variables relative to the m-th class. 
  train_outcomes <- list()
  counter <- 1
  for (m in sort(unique(y))) {
    train_outcomes[[counter]] <- ifelse(y == m, 1, 0)
    counter <- counter + 1
  }  

  ## 2.) Fit the estimator on each dummy and get predictions.
  if (learner == "forest") {
    estimates <- lapply(train_outcomes, function(x) {ranger::ranger(y = x, x = X, num.trees = 2000)})
    predictions <- lapply(estimates, function(x) {predict(x, X)$predictions}) 
  } else if (learner == "l1") {
    X_scaled <- scale(as.matrix(X))
    cv_lassos <- lapply(train_outcomes, function(outcome) {glmnet::cv.glmnet(x = X_scaled, y = outcome, alpha = 1, family = "binomial")})
    best_lambdas <- lapply(cv_lassos, function(x) {x$lambda.min})
    
    estimates <- list()
    for (m in sort(unique(y))) {
      estimates[[m]] <- glmnet::glmnet(x = X_scaled, y = train_outcomes[[m]], alpha = 1, family = "binomial", lambda = best_lambdas[[m]])
    }
    predictions <- lapply(estimates, function(model) {as.numeric(predict(model, X_scaled, type = "response"))}) 
  }
  
  ## 3.) Normalize and put into matrix.
  predictions_final <- sapply(predictions, function(x) as.matrix(x))
  predictions_final <- matrix(apply(predictions_final, 1, function(x) (x) / (sum(x))), ncol = n_categories, byrow = T)
  colnames(predictions_final) <- paste0("P(Y=", seq_len(n_categories), ")")
  
  ## 4.) Output.
  output <- list("estimators" = estimates,
                 "predictions" = predictions_final,
                 "learner" = learner,
                 "X" = X,
                 "y" = y)
  class(output) <- "mml"
  
  ## 8.) Output.
  return(output)
}