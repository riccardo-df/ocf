# Ordered Correlation Forest <a href="https://riccardo-df.github.io/ocf/"><img src="man/figures/logo.svg" align="right" height="130"/></a>

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT) [![CRAN status](https://www.r-pkg.org/badges/version/ocf)](https://CRAN.R-project.org/package=ocf) [![Downloads](https://cranlogs.r-pkg.org/badges/ocf)](https://CRAN.R-project.org/package=ocf)

## Overview

`ocf` implements the *ordered correlation forest* estimator, a machine learning estimator specifically designed for predictive modeling of ordered non-numeric outcomes.

The package delivers:

✔ **Forest-based estimation** of conditional choice probabilities.\
✔ **Marginal effects** of covariates on the outcome probabilities.\
✔ **Standard error estimation**, leveraging the weight-based structure of random forest predictions.

------------------------------------------------------------------------

## Why Use `ocf`?

| Feature                            | Benefit                                                                              |
|------------------------|------------------------------------------------|
| **Optimized for Ordered Outcomes** | Unlike traditional random forests, `ocf` correctly handles ordered categorical data. |
| **Interpretable Marginal Effects** | Understand how covariates correlate with outcome probabilities.                      |
| **Easy to Use**                    | Integrates seamlessly into existing machine learning workflows.                      |
| **Active Development & Support**   | Open-source and actively maintained.                                                 |

------------------------------------------------------------------------

## 🚀 Installation

To install the latest stable version from CRAN:

```         
install.packages("ocf")
```

Alternatively, the current development version of the package can be installed using the `devtools` package:

```         
devtools::install_github("riccardo-df/ocf") # run install.packages("devtools") if needed.
```

------------------------------------------------------------------------

## Contributing

We welcome contributions! If you encounter issues, have feature requests, or want to contribute to the package, please follow the guidelines below.

📌 **Report an Issue:** If you encounter a bug or have a suggestion, please open an issue on GitHub:\
[Submit an Issue](https://github.com/riccardo-df/ocf/issues)

📌 **Contribute Code:** We encourage contributions via pull requests. Before submitting, please:\
1. Fork the repository and create a new branch.\
2. Ensure that your code follows the existing style and documentation conventions.\
3. Run tests and check for package integrity.\
4. Submit a pull request with a clear description of your changes.

📌 **Feature Requests:** If you have ideas for new features or extensions, feel free to discuss them by opening an issue.

------------------------------------------------------------------------

## Citation

If you use `ocf` in your research, please cite the corresponding paper:

> Di Francesco, R. (2025). Ordered Correlation Forest. Econometric Reviews 44(4).
