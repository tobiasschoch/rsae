# rsae: Robust Small Area Estimation<img src="inst/varia/logo.svg" align="right" width=120 height=139 alt="" />

[![CRAN](https://www.r-pkg.org/badges/version/rsae)](https://cran.r-project.org/package=rsae)


## Summary

The `rsae` package is an `R` ([R Core Team, 2022](#references)) package that provides functions to estimate the parameters of the basic unit-level model in small area estimation (also known as model type "B" in [Rao, 2003](#references),  or nested-error regression model in [Battese et al., 1988](#references)). 

In step 1, the model is fitted by one of the methods:

* maximum likelihood (see e.g., [Rao, 2003](#references), chapter 7.2),
* M-estimation, which is robust against outliers; see [Schoch (2012)](#references).

In step 2, the area-specific means are predicted using the empirical best linear unbiased predictor (EBLUP) or a robust prediction method due to [Copt and Victoria-Feser (2009)](#references). In addition, the mean square prediction error of the area-specific means can be computed by a parametric bootstrap.


## Installation

The package can be installed from CRAN using
```
install.packages("rsae")
```

## Building

Make sure that the R package `devtools` is installed. Then, the `rsae` package can be pulled from this GitHub repository and installed by
```
devtools::install_github("tobiasschoch/rsae")
```

## Community guidelines

#### Submitting an issue

If you have any suggestions for feature additions or any problems with the software that you would like addressed with the development community, please submit an issue on the Issues tab of the project GitHub repository. You may want to search the existing issues before submitting, to avoid asking a question or requesting a feature that has already been discussed.

#### How to contribute

If you are interested in modifying the code, you may fork the project for your own use, as detailed in the FreeBSD License we have adopted for the project. In order to contribute, please contact the developer by Tobias Schoch at gmail dot com (the names are separated by a dot) after making the desired changes.

#### Asking for help

If you have questions about how to use the software, or would like to seek out collaborations related to this project, you may contact Tobias Schoch (see contact details above).

## References

BATTESE, G. E., R. M. HARTER, AND W. A. FULLER (1988). An error component model for prediction of county crop areas using, *Journal of the American Statistical Association* **83**, 28–36.

COPT, S. AND M.-P. VICTORIA-FESER (2009). *Robust prediction in mixed linear models*, Tech. report, University of Geneva.

R CORE TEAM (2022). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

RAO, J. (2003). *Small Area Estimation*, Hoboken (NJ): John Wiley and Sons.

SCHOCH, T. (2012). Robust Unit-Level Small Area Estimation: A Fast Algorithm for Large Data, *Austrian Journal of Statistics* **41**, 243–265.
