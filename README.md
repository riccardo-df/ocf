# Modified Ordered Random Forest
 
Nonparametric estimator of the ordered choice model using random forests. 

`morf` modifies a standard random forest splitting criterion to build a collection of forests, each estimating the conditional probability of a single class. Under an honesty condition, the predicted conditional probabilities are asymptotically normal and consistent. This allows us to conduct valid inference using standard methods, e.g., by constructing conventional confidence intervals. However, honesty generally comes at the expense of a larger mean squared error. Thus, if inference is not of interest, adaptive (i.e., non-honest) estimation is recommended.

The package also implements a nonparametric estimator of the covariates’ marginal effects. Following the approach of Lechner & Mareckova (2022), we can use the weights induced by each forest to estimate the variance of the marginal effects.

 ## Installation  
The current development version of the package can be installed using the `devtools` package:

```
devtools::install_github("riccardo-df/morf")
library(morf)
```

## Usage Examples
This section demonstrates how to use `morf` to estimate conditional choice probabilities and marginal effects, and conduct inference about these statistical targets. In the examples, I use a synthetic data set provided in the [orf](https://github.com/okasag/orf) package:

```
## Load data from orf package.
set.seed(1986)

library(orf)
data(odata)

y <- as.numeric(odata[, 1])
X <- as.matrix(odata[, -1])

## Training-test split.
train_idx <- sample(seq_len(length(y)), length(y)/2)

y_tr <- y[train_idx]
X_tr <- X[train_idx, ]

y_test <- y[-train_idx]
X_test <- X[-train_idx, ]
```

The `morf` function constructs a collection of forests, one for each category of `y` (three in this case). We can decide whether to grow honest forests by setting the argument `honesty`. Honesty is necessary to obtain asymptotically normal and consistent predictions. However, it generally comes at the expense of a larger mean squared error. Thus, if inference is not of interest, we recommend to set `honesty = FALSE` (the default).

```
## Fit morf on training sample.
forests <- morf(y_tr, X_tr)

## We have compatibility with generic S3-methods.
print(forests)
summary(forests)
predictions <- predict(forests, X_test)
head(predictions$probabilities)
table(y_test, predictions$classification)

## Compute standard errors. This requires honest forests.
honest_forests <- morf(y_tr, X_tr, honesty = TRUE, inference = TRUE)
honest_forests$predictions$standard.errors
```

To compute marginal effects, we can use the `marginal_effects` functions. This allows us to compute mean marginal effects, marginal effects at the mean, or marginal effects at the median, according to the `eval` argument. We recommend to grow a larger number of trees to obtain accurate estimates.

```
## Fit morf. Use large number of trees.
forests <- morf(y, X, n.trees = 4000)

## Marginal effects at the mean.
me <- marginal_effects(forests, eval = "atmean")
print(me)
summary(me)

## LATEX.
print(me, latex = TRUE)

## Compute standard errors. This requires honest forests.
honest_forests <- morf(y, X, n.trees = 4000, honesty = TRUE)
honest_me <- marginal_effects(honest_forests, eval = "atmean", inference = TRUE)
honest_me$standard.errors
honest_me$p.values # These are not corrected for multiple hypotheses testing!

print(honest_me, latex = TRUE)
```

## References

- Athey, S., Tibshirani, J., & Wager, S. (2019).
<b>Generalized Random Forests.</b> <i>Annals of Statistics</i>, 47(2).
[<a href="https://projecteuclid.org/euclid.aos/1547197251">paper</a>]

- Lechner, M., & Mareckova, J. (2022). 
<b>Modified Causal Forest.</b>
<i>arXiv preprint arXiv:2209.03744</i>.
[<a href="https://arxiv.org/abs/2209.03744">paper</a>]

- Lechner, M., & Okasa, G. (2019). 
<b>Random Forest Estimation of the Ordered Choice Model.</b>
<i>arXiv preprint arXiv:1907.02436</i>.
[<a href="https://arxiv.org/abs/1907.02436">paper</a>]

- Wager, S., & Athey, S. (2018).
<b>Estimation and Inference of Heterogeneous Treatment Effects using Random Forests.</b>
<i>Journal of the American Statistical Association</i>, 113(523).
[<a href="https://www.tandfonline.com/eprint/v7p66PsDhHCYiPafTJwC/full">paper</a>]

- Wright, M. N. & Ziegler, A. (2017).
<b>ranger: A fast implementation of random forests for high dimensional data in C++ and R.</b>
<i>Journal of Statistical Software</i>, 77(1).
[<a href="https://www.jstatsoft.org/article/view/v077i01">paper</a>]
