#' Renaming Variables for LATEX Usage (Internal Use)
#'
#' Renames variables where the character "_" is used, which causes clashes in LATEX. Useful for the \code{phased} print method.
#'
#' @param names string vector.
#'
#' @return
#' The renamed string vector. Strings where "_" is not found are not modified by \code{rename_latex}.
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
