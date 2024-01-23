\name{rsae-package}
\alias{rsae-package}
\alias{rsae}
\docType{package}
\title{Robust Small Area Estimation}
\description{
    The package implements methods to fit the basic unit-level model (also
    known as model type "B" in Rao, 2003, or nested-error regression model in
    Battese et al., 1988), to predict area-specific means by the empirical
    best linear unbiased predictor (EBLUP) or a robust prediction method
    (Copt and Victoria-Feser, 2009; Heritier et al., 2011), and to compute
    the mean square prediction error of the predicted area-level means
    by a parametric bootstrap (Sinha and Rao, 2009; see also Hall and Maiti,
    2006; Lahiri, 2003).

    The methods are discussed in Schoch (2012).
}
\details{
    \subsection{Implemented methods}{
        \itemize{
            \item maximum likelihood estimator
            \item Huber-type M-estimator (RML II of Richardson and Welsh,
                1995, not the method proposed in Sinha and Rao, 2009); see
                Schoch (2012) for details
        }
    }
    \subsection{How to}{
        Data analysis involves the following steps:
        \enumerate{
            \item prepare the data/ specify the model for estimation;
                see \code{\link[=saemodel]{saemodel()}}
            \item fit the model by various (robust) methods; see
                \code{\link[=fitsaemodel]{fitsaemodel()}}
            \item (robustly) predict the random effects and the area
                means; see \code{\link[=robpredict]{robpredict()}}
        }
    }
}
\references{
Battese, G. E., Harter, R. M., and W.A. Fuller (1988).
    An error component model for prediction of county crop areas using.
    \emph{Journal of the American Statistical Association} \bold{83}, 28--36.
    \doi{https://doi.org/10.2307/2288915}

Copt, S. and M.-P. Victoria-Feser (2009).
    \emph{Robust Predictions in Mixed Linear Models},
    Research Report, University of Geneva.

Lahiri, P. (2003).
    On the impact of bootstrap in survey sampling and small area estimation.
    \emph{Statistical Science} \bold{18}, 199--210.
    \doi{https://doi.org/10.1214/ss/1063994975}

Hall, P. and T. Maiti (2006).
    On parametric bootstrap methods for small area prediction.
    \emph{Journal of the Royal Statistical Society. Series B} \bold{68},
    221--238. \doi{https://doi.org/10.1111/j.1467-9868.2006.00541.x}

Heritier, S., Cantoni, E., Copt, S., and M.-P. Victoria-Feser (2009).
    \emph{Robust methods in Biostatistics}, New York: John Wiley and Sons.

Rao, J.N.K. (2003).
    \emph{Small Area Estimation}, New York: John Wiley and Sons.

Richardson, A.M. and A.H. Welsh (1995).
    Robust restricted maximum likelihood in mixed linear model.
    \emph{Biometrics} \bold{51}, 1429--1439.
    \doi{https://doi.org/10.2307/2533273}

Schoch, T. (2012). Robust Unit-Level Small Area Estimation: A Fast Algorithm
    for Large Datasets. \emph{Austrian Journal of Statistics} \bold{41},
    243--265. \doi{https://doi.org/10.17713/ajs.v41i4.1548}

Sinha, S.K. and J.N.K. Rao (2009).
    Robust small area estimation.
    \emph{Canadian Journal of Statistics} \bold{37}, 381--399.
    \doi{https://doi.org/10.1002/cjs.10029}
}
