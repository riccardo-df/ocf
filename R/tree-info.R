#' Tree Information in Readable Format
#' 
#' Extracts tree information from a \code{morf.forest} object. 
#' 
#' @param object \code{morf.forest} object.
#' @param tree Number of the tree of interest.
#' 
#' @return A \code{data.frame} with the following columns:
#'    \item{\code{nodeID}}{The node's ID.} 
#'    \item{\code{leftChild}}{ID of the left child node.} 
#'    \item{\code{rightChild}}{ID of the right child node.} 
#'    \item{\code{splitvarID}}{ID of the splitting variable.}
#'    \item{\code{splitvarName}}{Name of the splitting variable.}
#'    \item{\code{splitval}}{The splitting value.} 
#'    \item{\code{terminal}}{Logical, TRUE for terminal nodes.} 
#'    \item{\code{prediction}}{One column with the predicted conditional class probabilities.} 
#'   
#' @details 
#' Nodes' and variables' IDs are 0-indexed, i.e., node 0 is the root node. \cr
#' 
#' All values smaller or equal than \code{splitval} go to the left and all values larger go to the right. 
#' 
#' @seealso \code{\link{morf}}
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
tree_info <- function(object, tree = 1) {
  ## Handling inputs and checks.
  if (!inherits(object, "morf.forest")) stop("Invalid class of input object.", call. = FALSE)
  
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
  result
}
