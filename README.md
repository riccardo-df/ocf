# Modified Ordered Random Forest for Estimating the Ordered Choice Model
 
 'morf' modifies the splitting criterion and prediction method of a standard random forest to estimate conditional choice probabilities. It also implements a nonparametric estimators of the marginal effects. Standard errors are computed via a weight-based approach that, under some conditions such as honesty and subsampling, allows the user to conduct approximate inference.
 
 ## Installation  
The current development version the package can be installed using the `devtools` package:

```
devtools::install_github("riccardo-df/morf")
library(morf)
```

## Usage Examples
This section demonstrates how to use `morf` to estimate conditional choice probabilities and marginal effects, and conduct approximate inference about these statistical targets. Let us generate some data:

```
## Generate data.
set.seed(1986)

n <- 100
k <- 3

X <- data.frame(matrix(rnorm(n)))
y <- sample(c(1, 2, 3), n, replace = TRUE)
```
