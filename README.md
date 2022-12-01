# Modified Ordered Random Forest for Estimating the Ordered Choice Model
 
`morf` modifies the splitting criterion and prediction method of a standard random forest to estimate conditional choice probabilities. It also implements a nonparametric estimators of the marginal effects. Standard errors are computed via a weight-based approach that, under particular conditions such as honesty and subsampling, allows the user to conduct approximate inference.
 
 ## Installation  
The current development version the package can be installed using the `devtools` package:

```
devtools::install_github("riccardo-df/morf")
library(morf)
```

## Usage Examples
This section demonstrates how to use `morf` to estimate conditional choice probabilities and marginal effects, and conduct approximate inference about these statistical targets. In the examples, I use a synthetic data set provided in the [orf](https://github.com/okasag/orf) package:

```
## Generate data.
set.seed(1986)

library(orf)
data(odata)

X <- as.matrix(odata[, -1])
y <- as.numeric(odata[, 1])
```

First, I estimate conditional choice probabilities and marginal effects. Note the argument `honesty` set to `TRUE`. This is necessary to conduct valid inference later. For the marginal effect estimation, I set `eval = "atmean"` for the marginal effect at the mean. Alternatively, we could target mean marginal effect and marginal effect at the median.

```
## Build forests.
forests <- morf(X, y, honesty = TRUE)
print(forests)

## Estimate marginal effects at mean.
atmean_effects <- marginal_effects(forests, eval = "atmean", inference = TRUE)

```
