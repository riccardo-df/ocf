# Modified Ordered Random Forest <a href="https://riccardo-df.github.io/morf/"><img src="man/figures/logo.svg" align="right" height="130" /></a>
 
 <!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/morf)](https://CRAN.R-project.org/package=morf)
<!-- badges: end -->
 
R package to implement modified ordered random forests (MORF), a nonparametric estimator of the ordered choice model. Additionally, the package implements a nonparametric estimator of the marginal effects.

MORF modifies a standard random forest splitting criterion to build a collection of forests, each estimating the conditional probabilities of a single class. The estimator inherits the asymptotic properties of random forests. Thus, under an honesty condition (i.e., that different observations are used to place the splits and compute leaf predictions) the predicted conditional probabilities are asymptotically normal and consistent. The particular honesty implementation used by MORF allows for a weight-based estimation of the variance of both predicted probabilities and marginal effects.

To get started, please check the online [vignette](https://riccardo-df.github.io/morf/articles/morf-vignette.html) for a short tutorial.

## Installation  
The package can be downloaded from CRAN:

```
install.packages("morf")
```

Alternatively, the current development version of the package can be installed using the `devtools` package:

```
devtools::install_github("riccardo-df/morf") # run install.packages("devtools") if needed.
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
