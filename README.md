# Ordered Correlation Forest <a href="https://riccardo-df.github.io/ocf/"><img src="man/figures/logo.svg" align="right" height="130" /></a>

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/ocf)](https://CRAN.R-project.org/package=ocf)
![CRAN Downloads overall](http://cranlogs.r-pkg.org/badges/grand-total/ocf)
<!-- badges: end -->
 
R package to implement *ordered correlation forest*, a machine learning estimator specifically optimized for predictive modeling of ordered non-numeric outcomes. 

`ocf` provides forest-based estimation of the conditional choice probabilities and the covariates’ marginal effects. Under an "honesty" condition, the estimates are consistent and asymptotically normal and standard errors can be obtained by leveraging the weight-based representation of the random forest predictions. Please reference the use as [Di Francesco (2023)](https://arxiv.org/abs/2309.08755).

To get started, please check the online [short tutorial](https://riccardo-df.github.io/ocf/articles/ocf-short-tutorial.html).

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
<b>Generalized Random Forests.</b> 
<i>Annals of Statistics</i>, 47(2).
[<a href="https://doi.org/10.1214/18-AOS1709">paper</a>]

- Di Francesco, R. (forthcoming). 
<b>Ordered Correlation Forest.</b>
<i>Econometric Reviews</i>.
[<a href="https://doi.org/10.48550/arXiv.2309.08755">paper</a>]

- Lechner, M., & Mareckova, J. (2022). 
<b>Modified Causal Forest.</b>
<i>arXiv preprint arXiv:2209.03744</i>.
[<a href="https://doi.org/10.48550/arXiv.2209.03744">paper</a>]

- Lechner, M., & Okasa, G. (2019). 
<b>Random Forest Estimation of the Ordered Choice Model.</b>
<i>arXiv preprint arXiv:1907.02436</i>.
[<a href="https://doi.org/10.48550/arXiv.1907.02436">paper</a>]

- Peracchi, F. (2014). 
<b>Econometric methods for ordered responses: Some recent developments.</b>
<i>In Econometric methods and their applications in finance, macro and related fields(pp. 133–165). World Scientific</i>.
[<a href="https://doi.org/10.1142/9789814513470_0006">paper</a>]

- Wager, S., & Athey, S. (2018).
<b>Estimation and Inference of Heterogeneous Treatment Effects using Random Forests.</b>
<i>Journal of the American Statistical Association</i>, 113(523).
[<a href="https://doi.org/10.1080/01621459.2017.1319839">paper</a>]

- Wright, M. N. & Ziegler, A. (2017).
<b>ranger: A fast implementation of random forests for high dimensional data in C++ and R.</b>
<i>Journal of Statistical Software</i>, 77(1).
[<a href="https://doi.org/10.18637/jss.v077.i01">paper</a>]
