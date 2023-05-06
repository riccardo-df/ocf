#' Ordered Machine Learning
#' 
#' Estimation strategy to estimate conditional choice probabilities for ordered non-numeric outcomes.
#'
#' @param y Outcome vector.
#' @param X Covariate matrix (no intercept).
#' @param learner String, either \code{"forest"} or \code{"l1"}. Selects the base learner to estimate each expectation.
#' 
#' @return 
#' Object of class \code{oml}.
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
#' ## Fit ordered machine learning on training sample using two different learners.
#' ordered_forest <- ordered_ml(y_tr, X_tr, learner = "forest")
#' ordered_l1 <- ordered_ml(y_tr, X_tr, learner = "l1")
#' 
#' ## Predict out of sample.
#' predictions_forest <- predict(ordered_forest, X_test)
#' predictions_l1 <- predict(ordered_l1, X_test)
#' 
#' ## Compare predictions.
#' cbind(head(predictions_forest), head(predictions_l1))
#' 
#' @details
#' Ordered machine learning expresses conditional choice probabilities as the difference between the cumulative probabilities 
#' of two adjacent classes, which in turn can be expressed as conditional expectations of binary variables:
#' 
#' \deqn{p_m \left( X_i \right) = \mathbb{E} \left[ 1 \left( Y_i \leq m \right) | X_i \right] - \mathbb{E} \left[ 1 \left( Y_i \leq m - 1 \right) | X_i \right]}
#' 
#' Then we can separately estimate each expectation using any regression algorithm and pick the difference between the m-th and the
#' (m-1)-th estimated surfaces to estimate conditional probabilities.
#' 
#' \code{\link{ordered_ml}} combines this strategy with either regression forests or penalized logistic regression with an L1 penalty,
#' according to the user-specified parameter \code{learner}. If \code{learner == "forest"}, then the \code{\link[orf]{orf}}
#' function is called from an external package, as this estimator has already been proposed by Lechner and Okasa (2019). If 
#' \code{learner == "l1"}, the covariates are scaled to have zero mean and unit variance, and the penalty parameters are chosen 
#' via 10-fold cross-validation. Also, \code{\link[stats]{model.matrix}} is used to handle non-numeric covariates.\cr
#' 
#' 
#' @import ranger glmnet orf
#'  
#' @seealso \code{\link{multinomial_ml}}, \code{\link{ocf}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
ordered_ml <- function(y = NULL, X = NULL,
                           learner = "forest") {
  ## 0.) Handling inputs and checks.
  check_x_y(X, y)
  if (!(learner %in% c("forest", "l1"))) stop("Invalid 'learner'. This must be either 'forest' or 'l1'.", call. = FALSE)
  n <- length(y)
  n_categories <- length(unique(y))
  y_classes <- sort(unique(y))
  
  ## 1.) Fit the estimator and get predictions.
  if (learner == "forest") {
    estimates <- orf::orf(X, y, num.trees = 2000, honesty = FALSE)$forests
    predictions <- lapply(estimates, function(x) {predict(x, X)$predictions}) 
  } else if (learner == "l1") {
    # Generate binary outcomes for each class. The m-th element stores the indicator variables relative to the m-th class. We do not need the M-th element.
    train_outcomes <- list()
    counter <- 1
    for (m in y_classes[-n_categories]) {
      train_outcomes[[counter]] <- ifelse(y == m, 1, 0)
      counter <- counter + 1
    }  
    
    # Fit the estimator on each dummy.
    X_design <- stats::model.matrix(y ~ ., data = data.frame(y, X))[, -1]
    X_scaled <- scale(as.matrix(X_design))    
    cv_lassos <- lapply(train_outcomes, function(outcome) {glmnet::cv.glmnet(x = X_scaled, y = outcome, alpha = 1, family = "binomial")})
    best_lambdas <- lapply(cv_lassos, function(x) {x$lambda.min})
    
    estimates <- list()
    for (m in 1:(n_categories-1)) {
      estimates[[m]] <- glmnet::glmnet(x = X_scaled, y = train_outcomes[[m]], alpha = 1, family = "binomial", lambda = best_lambdas[[m]])
    }
    
    predictions <- lapply(estimates, function(model) {as.numeric(predict(model, X_scaled, type = "response"))}) 
  }
  
  ## 2.) Pick differences.
  predictions1 <- append(predictions, list(rep(1, n))) 
  predictions0 <- append(list(rep(0, n)), predictions) 
  differences <- as.list(mapply(function(x, y) x - y, predictions1, predictions0, SIMPLIFY = FALSE))
  
  ## 3.) Truncate, normalize and put into matrix.
  predictions_final <- lapply(differences, function(x) ifelse((x < 0), 0, x))
  predictions_final <- sapply(predictions_final, function(x) as.matrix(x))
  predictions_final <- matrix(apply(predictions_final, 1, function(x) (x) / (sum(x))), ncol = n_categories, byrow = T)
  colnames(predictions_final) <- paste0("P(Y=", seq_len(n_categories), ")")
  
  ## 4.) Output.
  output <- list("estimators" = estimates,
                 "predictions" = predictions_final,
                 "learner" = learner,
                 "X" = X,
                 "y" = y)
  class(output) <- "oml"
  
  ## 8.) Output.
  return(output)
}
