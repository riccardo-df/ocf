# Modified Ordered Random Forest for Estimating the Ordered Choice Model
 
`morf` modifies the splitting criterion and prediction method of a standard random forest to estimate conditional choice probabilities. It also implements a nonparametric estimators of the marginal effects. Standard errors are computed via a weight-based approach that, under particular conditions such as honesty and subsampling, allows the user to conduct approximate inference.
 
 ## Installation  
The current development version of the package can be installed using the `devtools` package:

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

The following chunk of code estimates conditional choice probabilities and marginal effects. Note the argument `honesty` set to `TRUE`. This is necessary to conduct valid inference. If inference is not of interest, one can improve the prediction performance by setting `honesty = FALSE`. For the marginal effect estimation, I set `eval = "atmean"` for the marginal effect at the mean. Alternatively, we could target mean marginal effect and marginal effect at the median.  The output of the `morf` function can be used with the usual generic methods, such as `summary` and `predict`.

```
## Build forests.
forests <- morf(X, y, honesty = TRUE)
print(forests)

## Estimate marginal effects at mean.
atmean_effects <- marginal_effects(forests, eval = "atmean", inference = TRUE)
print(atmean_effects, latex = TRUE)
```

## References

Lechner, M., & Mareckova, J. (2022). 
<b>Modified Causal Forest.</b>
<i>arXiv preprint arXiv:2209.03744</i>, 2016.
[<a href="https://arxiv.org/abs/2209.03744">paper</a>]

Wright, M. N. & Ziegler, A. (2017).
<b>ranger: A fast implementation of random forests for high dimensional data in C++ and R.</b>
<i>Journal of Statistical Software</i>, 2017, 77(1).
[<a href="https://www.jstatsoft.org/article/view/v077i01">paper</a>]
