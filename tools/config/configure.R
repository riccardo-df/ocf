# Prepare your package for installation here.
# Use 'define()' to define configuration variables.
# Use 'configure_file()' to substitute configuration values.

library(tools)

for (std in c('CXX20', 'CXX17', 'CXX14')) {
  cat('Checking if R knows a', std, 'compiler... ')
  out <- suppressWarnings(Rcmd(paste('config', std)))
  if (out == 0) break
}

if (out != 0) stop("Couldn't find a C++ >= 14 compiler")

f <- file(file.path('src', 'Makevars'), 'wb')
writeLines(c(paste('CXX_STD =', std)), f)
close(f)

