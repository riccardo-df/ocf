# Ordered Correlation Forest <a href="https://riccardo-df.github.io/ocf/"><img src="man/figures/logo.svg" align="right" height="130" /></a>

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/ocf)](https://CRAN.R-project.org/package=ocf)
![CRAN Downloads overall](http://cranlogs.r-pkg.org/badges/grand-total/ocf)
<!-- badges: end -->
 
R package to implement ordered correlation forests (OCF), a nonparametric estimator specifically optimized for handling ordered non-numeric outcomes. 

OCF modifies a standard random forest splitting criterion to build a collection of forests, each estimating the conditional probabilities of a single class. Under an \open honesty" condition, the estimator inherits the asymptotic properties of random forests, namely the consistency and asymptotic normality of their predictions. The particular honesty implementation used by OCF allows us to obtain standard errors for the covariates' marginal effects. The estimated standard errors can then be used to construct conventional confidence intervals.

To get started, please check the online [vignette](https://riccardo-df.github.io/ocf/articles/ocf-short-tutorial.html) for a short tutorial.

## Installation  
The package can be downloaded from CRAN:

```
install.packages("ocf")
```

Alternatively, the current development version of the package can be installed using the `devtools` package:

```
devtools::install_github("riccardo-df/ocf") # run install.packages("devtools") if needed.
```

## References

- Athey, S., Tibshirani, J., & Wager, S. (2019).
<b>Generalized Random Forests.</b> <i>Annals of Statistics</i>, 47(2).
[<a href="https://projecteuclid.org/euclid.aos/1547197251">paper</a>]

- Di Francesco, R. (2023). 
<b>Ordered Correlation Forest.</b>
<i>arXiv preprint arXiv:2309.08755</i>.
[<a href="https://arxiv.org/abs/2309.08755">paper</a>]

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
