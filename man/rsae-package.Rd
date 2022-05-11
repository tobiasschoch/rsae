\name{rsae-package}
\alias{rsae-package}
\alias{rsae}
\docType{package}
\title{Robust Small Area Estimation}
\description{
Computes robust basic unit- and area-level and predicts area-specific means
}
\details{
    \subsection{Implemented methods}{
        \itemize{
            \item maximum likelihood
            \item Huber-type M-estimation
        }
    }
    \subsection{How to}{
        Data analysis with \code{rsae} involves the following steps:
        \enumerate{
            \item prepare the data/ set up the model for estimation;
                see \code{\link{saemodel}}
            \item fit the model by various (robust) methods; see
                \code{\link{fitsaemodel}}
            \item (robustly) predict the random effects and the area
                means; see \code{\link{robpredict}}
        }
    }
}
\references{
Anderson, E., Bai, Z., Bischof, C., Blackford, L. S., Demmel, J.,
    Dongarra, J., et al. (2000).
    \emph{LAPACK users' guide} (3rd ed.).
    Philadelphia: Society for Industrial and Applied Mathematics (SIAM).

Battese, G. E., Harter, R. M., and W.A. Fuller (1988).
    An error component model for prediction of county crop areas using.
    \emph{Journal of the American Statistical Association} \bold{83}, 28-36.
    \doi{https://doi.org/10.2307/2288915}

Blackford, L.S., Petitet, A., Pozo, R., Remington, K., Whaley, R.C.,
    Demmel, J., et al. (2002).
    An updated set of basic linear algebra subprograms (BLAS).
    \emph{ACM Transactions on Mathematical Software} \bold{28}, 135-151.
    \doi{https://doi.org/10.1145/567806.567807}

Brent, R.P. (1973).
    \emph{Algorithms for minimization without derivatives}.
    Englewood Cliffs, NJ: Prentice-Hall.

Lahiri, P. (2003).
    On the impact of bootstrap in survey sampling and small area estimation.
    \emph{Statistical Science} \bold{18}, 199-210.
    \doi{https://doi.org/10.1214/ss/1063994975}

Hall, P. and T. Maiti (2006).
    On parametric bootstrap methods for small area prediction.
    \emph{Journal of the Royal Statistical Society. Series B} \bold{68},
    221-238. \doi{https://doi.org/10.1111/j.1467-9868.2006.00541.x}

Heritier, S., Cantoni, E., Copt, S., and M.-P. Victoria-Feser (2009).
    \emph{Robust methods in Biostatistics}, New York: John Wiley and Sons.

Maronna, R.A., Martin, D., and V.J. Yohai (2006).
    \emph{Robust statistics: Theory and methods}. Chichester: John Wiley.

Rao, J.N.K. (2003).
    \emph{Small Area Estimation}, New York: John Wiley and Sons.

Richardson, A.M. and A.H. Welsh (1995).
    Robust restricted maximum likelihood in mixed linear model.
    \emph{Biometrics} \bold{51}, 1429-1439.
    \doi{https://doi.org/10.2307/2533273}

Rousseeuw, P. J. and K. Van Driessen (2006).
    Computing LTS regression for large data sets.
    \emph{Data Mining and Knowledge Discovery} \bold{12}, 29-45.
    \doi{https://doi.org/10.1007/s10618-005-0024-4}

Schoch, T. (2012). Robust Unit-Level Small Area Estimation: A Fast Algorithm
    for Large Datasets. \emph{Austrian Journal of Statistics} \bold{41},
    243-265. \url{http://www.stat.tugraz.at/AJS/ausg124/124Schoch.pdf}

Sinha, S.K. and J.N.K. Rao (2009).
    Robust small area estimation.
    \emph{Canadian Journal of Statistics} \bold{37}, 381-399.
    \doi{https://doi.org/10.1002/cjs.10029}

Stahel, W. A. and A. Welsh (1997).
    Approaches to robust estimation in the simplest variance components model.
    \emph{Journal of Statistical Planning and Inference} \bold{57}, 295-319.
    \doi{https://doi.org/10.1016/S0378-3758(96)00050-X}
}
