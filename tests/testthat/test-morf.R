test_that("morf splits and predicts as expected with continuos covariates", {
  ## Generating data.
  set.seed(rnorm(1, sd = 1000)) # Random seed.
  
  n <- 1000
  m <- sample(c(1, 2, 3), size = 1) # Class to be tested.
  
  y <- sample(c(1, 2, 3), size = n, replace = TRUE)
  x <- data.frame("x1" = rnorm(n))
  
  y_m <- ifelse(y <= m, 1, 0)
  y_m_1 <- ifelse(y <= m-1, 1, 0)
  
  alpha <- runif(1, 0, 0.5)
  
  ## Fitting a "stump".
  morf <- morf(x = x, y = y, n.trees = 1, max.depth = 1, replace = FALSE, sample.fraction = 1, min.node.size = 1, 
               honesty = FALSE, alpha = alpha)
  
  avg_split <- tree_info(morf[[m]])$splitval[1] 
  predictions <- tree_info(morf[[m]])$prediction[-1]
  split_values <- combn(x[, 1], 2)[, which(avg_split == combn(x[, 1], 2, mean))]
  
  ## R splitting criterion.
  modified_split <- function(x, y_m, y_m_1, alpha) {
    splits <- sort(unique(x))
    mse <- rep(NA, length(splits))
    
    ## Scanning all split points x.
    for (i in seq_along(splits)) {
      ## Skip this split value if alpha-regularity would be violated.
      if (sum(x < splits[i]) < length(x) * alpha | sum(x > splits[i]) > length(x) - (length(x) * alpha)) next
      
      split <- splits[i]
      
      mse_m <- sum(sum(y_m[ x < split])^2 / sum(x < split), sum(y_m[x >= split])^2 / sum(x >= split), na.rm = TRUE)
      mse_m_1 <- sum(sum(y_m_1[ x < split])^2 / sum(x < split), sum(y_m_1[x >= split])^2 / sum(x >= split), na.rm = TRUE)
      
      mce <- sum(mean(y_m[x < split] * y_m_1[x < split]), -mean(y_m[x < split]) * mean(y_m_1[x < split]),
                 mean(y_m[x >= split] * y_m_1[x >= split]), -mean(y_m[x >= split]) * mean(y_m_1[x >= split]), na.rm = TRUE)
      
      mse[i] <- mse_m + mse_m_1 - 2 * mce
    }
    
    ## Best split.
    best_split <- splits[which.max(mse)]
    
    left_prediction <- mean(y_m[x < best_split]) - mean(y_m_1[x < best_split])
    right_prediction <- mean(y_m[x >= best_split]) - mean(y_m_1[x >= best_split])
    
    predictions <- c(left_prediction, right_prediction)
    
    return(list("best_split" = best_split,
                "predictions" = predictions))
  }
  
  ## Comparing.
  treeR <- modified_split(x[, 1], y_m, y_m_1, alpha)
  
  check_split <- treeR$best_split %in% split_values
  
  expect_true(check_split)
  expect_setequal(treeR$predictions, predictions)
})


test_that("morf splits and predicts as expected with categorical covariates", {
  ## Generating data.
  set.seed(rnorm(1, sd = 1000)) # Random seed.
  
  n <- 1000
  m <- sample(c(1, 2, 3), size = 1) # Class to be tested.
  
  y <- sample(c(1, 2, 3), size = n, replace = TRUE)
  x <- data.frame("x1" = sample(c(1, 2, 3, 4, 5), size = n, replace = TRUE))
  
  y_m <- ifelse(y <= m, 1, 0)
  y_m_1 <- ifelse(y <= m-1, 1, 0)
  
  alpha <- runif(1, 0, 0.5)
  
  ## Fitting a "stump".
  morf <- morf(x = x, y = y, n.trees = 1, max.depth = 1, replace = FALSE, sample.fraction = 1, min.node.size = 1, 
               honesty = FALSE, alpha = alpha)
  
  avg_split <- tree_info(morf[[m]])$splitval[1] 
  predictions <- tree_info(morf[[m]])$prediction[-1]
  split_values <- combn(x[, 1], 2)[, which(avg_split == combn(x[, 1], 2, mean))]
  
  ## R splitting criterion.
  modified_split <- function(x, y_m, y_m_1, alpha) {
    splits <- sort(unique(x))
    mse <- rep(NA, length(splits))
    
    ## Scanning all split points x.
    for (i in seq_along(splits)) {
      ## Skip this split value if alpha-regularity would be violated.
      if (sum(x < splits[i]) < length(x) * alpha | sum(x > splits[i]) > length(x) - (length(x) * alpha)) next
      
      split <- splits[i]
      
      mse_m <- sum(sum(y_m[ x < split])^2 / sum(x < split), sum(y_m[x >= split])^2 / sum(x >= split), na.rm = TRUE)
      mse_m_1 <- sum(sum(y_m_1[ x < split])^2 / sum(x < split), sum(y_m_1[x >= split])^2 / sum(x >= split), na.rm = TRUE)
      
      mce <- sum(mean(y_m[x < split] * y_m_1[x < split]), -mean(y_m[x < split]) * mean(y_m_1[x < split]),
                 mean(y_m[x >= split] * y_m_1[x >= split]), -mean(y_m[x >= split]) * mean(y_m_1[x >= split]), na.rm = TRUE)
      
      mse[i] <- mse_m + mse_m_1 - 2 * mce
    }
    
    ## Best split.
    best_split <- splits[which.max(mse)]
    
    left_prediction <- mean(y_m[x < best_split]) - mean(y_m_1[x < best_split])
    right_prediction <- mean(y_m[x >= best_split]) - mean(y_m_1[x >= best_split])
    
    predictions <- c(left_prediction, right_prediction)
    
    return(list("best_split" = best_split,
                "predictions" = predictions))
  }
  
  ## Comparing.
  treeR <- modified_split(x[, 1], y_m, y_m_1, alpha)
  
  check_split <- treeR$best_split %in% split_values
  
  expect_true(check_split)
  expect_setequal(treeR$predictions, predictions)
})


test_that("Standard predictions and weight-based predictions are the same", {
  ## Generating data.
  set.seed(rnorm(1, sd = 1000)) # Random seed.
  
  n <- 1000

  y <- sample(c(1, 2, 3), size = n, replace = TRUE)
  x <- data.frame("x1" = rnorm(n))
  
  ## Fitting morf objects.
  set.seed(1986) # Set seed to get same honest split.
  morf <- morf(x = x, y = y, inference = FALSE)
  set.seed(1986)
  morf2 <- morf(x = x, y = y, inference = TRUE)
  
  ## Comparing.
  # expect_equal(sum(round(morf$predictions, 3) == round(morf2$predictions, 3)), n * 3)
  expect_setequal(round(morf$predictions, 3), round(morf2$predictions, 3))
})
