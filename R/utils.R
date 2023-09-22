#' Tree Information in Readable Format
#' 
#' Extracts tree information from a \code{ocf.forest} object. 
#' 
#' @param object \code{ocf.forest} object.
#' @param tree Number of the tree of interest.
#' 
#' @return A \code{data.frame} with the following columns:
#'    \item{\code{nodeID}}{Node IDs.} 
#'    \item{\code{leftChild}}{IDs of the left child node.} 
#'    \item{\code{rightChild}}{IDs of the right child node.} 
#'    \item{\code{splitvarID}}{IDs of the splitting variable.}
#'    \item{\code{splitvarName}}{Name of the splitting variable.}
#'    \item{\code{splitval}}{Splitting value.} 
#'    \item{\code{terminal}}{Logical, TRUE for terminal nodes.} 
#'    \item{\code{prediction}}{One column with the predicted conditional class probabilities.} 
#'    
#' @examples 
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#' 
#' data <- generate_ordered_data(1000)
#' sample <- data$sample
#' Y <- sample$Y
#' X <- sample[, -1]
#' 
#' ## Fit ocf.
#' forests <- ocf(Y, X)
#' 
#' ## Extract information from tenth tree of first forest.
#' info <- tree_info(forests$forests.info$forest.1, tree = 10)
#' head(info)}
#'   
#' @details 
#' Nodes and variables IDs are 0-indexed, i.e., node 0 is the root node. \cr
#' 
#' All values smaller than or equal to \code{splitval} go to the left and all values larger go to the right. 
#' 
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2023). Ordered Correlation Forest. arXiv preprint \href{https://arxiv.org/abs/2309.08755}{arXiv:2309.08755}.
#' }
#'
#' @seealso \code{\link{ocf}}
#' 
#' @export
tree_info <- function(object, tree = 1) {
  ## Handling inputs and checks.
  if (!inherits(object, "ocf.forest")) stop("Invalid class of input object.", call. = FALSE)
  
  forest <- object
  
  if (is.null(forest$num.trees) ||
      is.null(forest$child.nodeIDs) || is.null(forest$split.varIDs) ||
      is.null(forest$split.values) || is.null(forest$covariate.names) ||
      is.null(forest$treetype)) stop("Invalid forest object.", call. = FALSE)
  
  if (tree > forest$num.trees) stop("Requesting tree ", tree, ", but forest has only ", forest$num.trees," trees.", call. = FALSE)
  
  ## Handling output.
  result <- data.frame(nodeID = 0:(length(forest$split.values[[tree]]) - 1),
                       leftChild = forest$child.nodeIDs[[tree]][[1]], 
                       rightChild = forest$child.nodeIDs[[tree]][[2]], 
                       splitvarID = forest$split.varIDs[[tree]], 
                       splitvarName = "X",
                       splitval = forest$split.values[[tree]], 
                       terminal = FALSE)
  
  result$leftChild[result$leftChild == 0] <- NA
  result$rightChild[result$rightChild == 0] <- NA
  result$terminal[is.na(result$leftChild)] <- TRUE
  result$splitvarID[result$terminal] <- NA
  result$splitvarName[result$terminal] <- NA
  result$splitval[result$terminal] <- NA
  result$splitvarName <- forest$independent.variable.names[result$splitvarID + 1]
  
  # Unordered splitting.
  idx.unordered <- !result$terminal & !forest$is.ordered[result$splitvarID + 1]
  if (any(idx.unordered)) {
    if (any(result$splitval[idx.unordered] > (2^31 - 1))) {
      warning("Unordered splitting levels can only be shown for up to 31 levels.")
      result$splitval[idx.unordered] <- NA
    } else {
      result$splitval[idx.unordered] <- sapply(result$splitval[idx.unordered], function(x) {
        paste(which(as.logical(intToBits(x))), collapse = ",")
      })
    }
  }
  
  # Prediction.
  result$prediction <- forest$split.values[[tree]]
  result$prediction[!result$terminal] <- NA
  
  ## Output.
  return(result)
}


#' Generate Ordered Data
#' 
#' Generate a synthetic data set with an ordered non-numeric outcome, together with conditional probabilities and covariates' marginal effects.
#' 
#' @param n Sample size.
#' 
#' @return A list storing a data frame with the observed data, a matrix of true conditional probabilities, 
#' and a matrix of true marginal effects at the mean of the covariates.
#'    
#' @examples 
#' \donttest{## Generate synthetic data.
#' set.seed(1986)
#' 
#' data <- generate_ordered_data(1000)
#' 
#' head(data$true_probs)
#' data$me_at_mean
#' 
#' sample <- data$sample
#' Y <- sample$Y
#' X <- sample[, -1]
#' 
#' ## Fit ocf.
#' forests <- ocf(Y, X)}
#'   
#' @details 
#' First, a latent outcome is generated as follows:
#' 
#' \deqn{Y_i^* = g ( X_i ) + \epsilon_i}
#' 
#' with:
#' 
#' \deqn{g ( X_i ) = X_i^T \beta}
#' 
#' \deqn{X_i := (X_{i, 1}, X_{i, 2}, X_{i, 3}, X_{i, 4}, X_{i, 5}, X_{i, 6})}
#' 
#' \deqn{X_{i, 1}, X_{i, 3}, X_{i, 5} \sim \mathcal{N} \left( 0, 1 \right)}
#' 
#' \deqn{X_{i, 2}, X_{i, 4}, X_{i, 6} \sim \textit{Bernoulli} \left( 0, 1 \right)}
#' 
#' \deqn{\beta = \left( 1, 1, 1/2, 1/2, 0, 0 \right)}
#' 
#' \deqn{\epsilon_i \sim logistic (0, 1)}
#'
#' Second, the observed outcomes are obtained by discretizing the latent outcome into three classes using uniformly spaced threshold parameters.
#' 
#' Third, the conditional probabilities and the covariates' marginal effects at the mean are generated using standard textbook formulas. Marginal
#' effects are approximated using a sample of 1,000,000 observations.
#' 
#' @author Riccardo Di Francesco
#' 
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2023). Ordered Correlation Forest. arXiv preprint \href{https://arxiv.org/abs/2309.08755}{arXiv:2309.08755}.
#' }
#' 
#' @seealso \code{\link{ocf}}
#' 
#' @export
generate_ordered_data <- function(n) {
  ## 0) Handling inputs and checks.
  if (n <= 0 | n %% 1 != 0) stop("Invalid 'n'. This must be a positive integer.", call. = FALSE)
  
  n_pop <- 1000000
  
  ## 1.) Generate covariates. By DGP, we know that x1, x3, and x5 have zero mean, and that x2, x4, and X6 have mean equal to 0.4.
  # 1a.) Signal covariates.
  X_signal_tr <- data.frame("x1" = rnorm(n), "x2" = rbinom(n, 1, 0.4), "x3" = rnorm(n), "x4" = rbinom(n, 1, 0.4))
  X_signal_pop <- data.frame("x1" = rnorm(n_pop), "x2" = rbinom(n_pop, 1, 0.4), "x3" = rnorm(n_pop), "x4" = rbinom(n_pop, 1, 0.4))
  
  # 1b.) Noise covariates.
  X_noise_tr <- data.frame("x5" = rnorm(n), "x6" = rbinom(n, 1, 0.4))
  X_noise_pop <- data.frame("x5" = rnorm(n_pop), "x6" = rbinom(n_pop, 1, 0.4))
  
  # 1c.) Bind columns.
  X_tr <- as.matrix(data.frame(X_signal_tr, X_noise_tr))
  X_pop <- as.matrix(data.frame(X_signal_pop, X_noise_pop))
  
  ## 2.) Generate latent outcome.
  # 2a.) Error term.
  epsilon_tr <- rlogis(n)
  epsilon_pop <- rlogis(n_pop)
  
  # 2b.) Regression function.
  betas <- c(1, 1, 1/2, 1/2, 0, 0)
  
  g_function_tr <- X_tr[, 1] * betas[1] + X_tr[, 2] * betas[2] + X_tr[, 3] * betas[3] + X_tr[, 4] * betas[4] + X_tr[, 5] * betas[5] + X_tr[, 6] * betas[6]
  g_function_pop <- X_pop[, 1] * betas[1] + X_pop[, 2] * betas[2] + X_pop[, 3] * betas[3] + X_pop[, 4] * betas[4] + X_pop[, 5] * betas[5] + X_pop[, 6] * betas[6]
  g_function_pop_at_mean <- 0 * betas[1] + 0.4 * betas[2] + 0 * betas[3] + 0.4 * betas[4] + 0 * betas[5] + 0.4 * betas[6]

  # 2c.) Latent outcome.
  Y_star_tr <- g_function_tr + epsilon_tr
  Y_star_pop <- g_function_pop + epsilon_pop
  
  ## 3.) Threshold parameters. First, define the cumulative probability of each class. Then, find the corresponding quantiles.
  probabilities <- seq(1, 2 ,1) / 3
  zetas <- quantile(Y_star_pop, probabilities)
  
  ## 4.) True choice probabilities.
  real_probs_tr <- cbind(plogis(zetas[1] - g_function_tr),
                         plogis(zetas[2] - g_function_tr) - plogis(zetas[1] - g_function_tr),
                         1 - plogis(zetas[2] - g_function_tr))
  colnames(real_probs_tr) <- paste0("P(Y=", 1:3, ")")
  
  # 5a.) For continuous covariates, compute the derivative of regression function wrt to x_j and evaluate it at the mean of the covariates.
  nabla_x1_g_at_mean <- betas[1]
  nabla_x3_g_at_mean <- betas[3]
  nabla_x5_g_at_mean <- betas[5]
  
  # 5b.) For discrete covariates, evaluate regression function at ceiling{\bar{w_j}} and the floor{\bar{w_j}}. 
  g_function_pop_at_mean_ceil_x2 <- 0 * betas[1] + 1 * betas[2] + 0 * betas[3] + 0.4 * betas[4] + 0 * betas[5] + 0.4 * betas[6] 
  g_function_pop_at_mean_ceil_x4 <- 0 * betas[1] + 0.4 * betas[2] + 0 * betas[3] + 1 * betas[4] + 0 * betas[5] + 0.4 * betas[6] 
  g_function_pop_at_mean_ceil_x6 <- 0 * betas[1] + 0.4 * betas[2] + 0 * betas[3] + 0.4 * betas[4] + 0 * betas[5] + 1 * betas[6] 
  g_function_pop_at_mean_floor_x2 <- 0 * betas[1] + 0 * betas[2] + 0 * betas[3] + 0.4 * betas[4] + 0 * betas[5] + 0.4 * betas[6] 
  g_function_pop_at_mean_floor_x4 <- 0 * betas[1] + 0.4 * betas[2] + 0 * betas[3] + 0 * betas[4] + 0 * betas[5] + 0.4 * betas[6] 
  g_function_pop_at_mean_floor_x6 <- 0 * betas[1] + 0.4 * betas[2] + 0 * betas[3] + 0.4 * betas[4] + 0 * betas[5] + 0 * betas[6] 
  
  g_function_pop_at_median_ceil_x2 <- 0 * betas[1] + 1 * betas[2] + 0 * betas[3] + 0 * betas[4] + 0 * betas[5] + 0 * betas[6] 
  g_function_pop_at_median_ceil_x4 <- 0 * betas[1] + 0 * betas[2] + 0 * betas[3] + 1 * betas[4] + 0 * betas[5] + 0 * betas[6] 
  g_function_pop_at_median_ceil_x6 <- 0 * betas[1] + 0 * betas[2] + 0 * betas[3] + 0 * betas[4] + 0 * betas[5] + 1 * betas[6] 
  g_function_pop_at_median_floor_x2 <- 0 * betas[1] + 0 * betas[2] + 0 * betas[3] + 0 * betas[4] + 0 * betas[5] + 0 * betas[6] 
  g_function_pop_at_median_floor_x4 <- 0 * betas[1] + 0 * betas[2] + 0 * betas[3] + 0 * betas[4] + 0 * betas[5] + 0 * betas[6] 
  g_function_pop_at_median_floor_x6 <- 0 * betas[1] + 0 * betas[2] + 0 * betas[3] + 0 * betas[4] + 0 * betas[5] + 0 * betas[6] 
  
  # 5c.) Matrix of marginal effects at mean.
  effect_class1_x1_at_mean <- - nabla_x1_g_at_mean * (dlogis(zetas[1] - g_function_pop_at_mean)) 
  effect_class1_x2_at_mean <- plogis(zetas[1] - g_function_pop_at_mean_ceil_x2) - plogis(zetas[1] - g_function_pop_at_mean_floor_x2)
  effect_class1_x3_at_mean <- - nabla_x3_g_at_mean * (dlogis(zetas[1] - g_function_pop_at_mean)) 
  effect_class1_x4_at_mean <- plogis(zetas[1] - g_function_pop_at_mean_ceil_x4) - plogis(zetas[1] - g_function_pop_at_mean_floor_x4)
  effect_class1_x5_at_mean <- - nabla_x5_g_at_mean * (dlogis(zetas[1] - g_function_pop_at_mean)) 
  effect_class1_x6_at_mean <- plogis(zetas[1] - g_function_pop_at_mean_ceil_x6) - plogis(zetas[1] - g_function_pop_at_mean_floor_x6)
  
  effect_class2_x1_at_mean <- - nabla_x1_g_at_mean * (dlogis(zetas[2] - g_function_pop_at_mean) - dlogis(zetas[1] - g_function_pop_at_mean)) 
  effect_class2_x2_at_mean <- plogis(zetas[2] - g_function_pop_at_mean_ceil_x2) - plogis(zetas[1] - g_function_pop_at_mean_ceil_x2) - (plogis(zetas[2] - g_function_pop_at_mean_floor_x2) - plogis(zetas[1] - g_function_pop_at_mean_floor_x2))
  effect_class2_x3_at_mean <- - nabla_x3_g_at_mean * (dlogis(zetas[2] - g_function_pop_at_mean) - dlogis(zetas[1] - g_function_pop_at_mean)) 
  effect_class2_x4_at_mean <- plogis(zetas[2] - g_function_pop_at_mean_ceil_x4) - plogis(zetas[1] - g_function_pop_at_mean_ceil_x4) - (plogis(zetas[2] - g_function_pop_at_mean_floor_x4) - plogis(zetas[1] - g_function_pop_at_mean_floor_x4))
  effect_class2_x5_at_mean <- - nabla_x5_g_at_mean * (dlogis(zetas[2] - g_function_pop_at_mean) - dlogis(zetas[1] - g_function_pop_at_mean)) 
  effect_class2_x6_at_mean <- plogis(zetas[2] - g_function_pop_at_mean_ceil_x6) - plogis(zetas[1] - g_function_pop_at_mean_ceil_x6) - (plogis(zetas[2] - g_function_pop_at_mean_floor_x6) - plogis(zetas[1] - g_function_pop_at_mean_floor_x6))
  
  effect_class3_x1_at_mean <- nabla_x1_g_at_mean * (dlogis(zetas[2] - g_function_pop_at_mean)) 
  effect_class3_x2_at_mean <- 1 - plogis(zetas[2] - g_function_pop_at_mean_ceil_x2) - (1 - plogis(zetas[2] - g_function_pop_at_mean_floor_x2))
  effect_class3_x3_at_mean <- nabla_x3_g_at_mean * (dlogis(zetas[2] - g_function_pop_at_mean)) 
  effect_class3_x4_at_mean <- 1 - plogis(zetas[2] - g_function_pop_at_mean_ceil_x4) - (1 - plogis(zetas[2] - g_function_pop_at_mean_floor_x4))
  effect_class3_x5_at_mean <- nabla_x5_g_at_mean * (dlogis(zetas[2] - g_function_pop_at_mean)) 
  effect_class3_x6_at_mean <- 1 - plogis(zetas[2] - g_function_pop_at_mean_ceil_x6) - (1 - plogis(zetas[2] - g_function_pop_at_mean_floor_x6))
  
  marginal_effects_at_mean <- matrix(c(effect_class1_x1_at_mean, effect_class2_x1_at_mean, effect_class3_x1_at_mean, 
                                       effect_class1_x2_at_mean, effect_class2_x2_at_mean, effect_class3_x2_at_mean,
                                       effect_class1_x3_at_mean, effect_class2_x3_at_mean, effect_class3_x3_at_mean,
                                       effect_class1_x4_at_mean, effect_class2_x4_at_mean, effect_class3_x4_at_mean,
                                       effect_class1_x5_at_mean, effect_class2_x5_at_mean, effect_class3_x5_at_mean,
                                       effect_class1_x6_at_mean, effect_class2_x6_at_mean, effect_class3_x6_at_mean), 
                                     ncol = 3, byrow = TRUE)
  colnames(marginal_effects_at_mean) <- paste("P'(Y=", 1:3, ")", sep = "")
  rownames(marginal_effects_at_mean) <- colnames(X_tr)
  
  ## 7.) Generate observed outcomes.
  Y_obs_tr <- rep(NA, n) 

  Y_obs_tr[Y_star_tr <= zetas[1]] <- 1 
  Y_obs_tr[zetas[1] < Y_star_tr & Y_star_tr <= zetas[2]] <- 2
  Y_obs_tr[zetas[2] < Y_star_tr] <- 3
  
  ## 8.) Output.
  dta_tr <- data.frame("Y" = Y_obs_tr, X_tr) 

  return(list("sample" = dta_tr,
              "true_probs" = real_probs_tr,
              "me_at_mean" = marginal_effects_at_mean))
}


#' Renaming Variables for LATEX Usage
#'
#' Renames variables where the character "_" is used, which causes clashes in LATEX. Useful for the \code{phased} print method.
#'
#' @param names string vector.
#'
#' @return
#' The renamed string vector. Strings where "_" is not found are not modified by \code{rename_latex}.
#' 
#' @keywords internal
rename_latex <- function(names) {
  ## Locating variables that need renaming.
  idx <- grepl("_", names, fixed = TRUE)
  
  if (sum(idx) == 0) return(names)
  
  ## Renaming variables.
  split_names <- stringr::str_split(string = names[idx], pattern = "_", simplify = TRUE)
  attach_names <- paste(split_names[, 1], split_names[, 2], sep = "\\_")
  
  ## Replacing.
  names[idx] <- attach_names
  
  ## Output.
  return(names)
}
