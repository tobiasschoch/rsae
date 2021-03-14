\name{rsae-package}
\alias{rsae-package}
\alias{rsae}
\docType{package}
\title{
Robust Small Area Estimation
}
\description{
Computes robust basic unit- and area-level and predicts area-specific means
}
\details{
\tabular{ll}{
Package: \tab rsae\cr
Type: \tab Package\cr
Version: \tab 0.1-5\cr
Date: \tab 2014-02-12\cr
Suggests: \tab robustbase, nlme\cr
License: \tab GPL (>=2) | FreeBSD\cr
LazyLoad: \tab yes\cr
}

    \subsection{Implemented methods}{
        \itemize{
            \item maximum likelihood (as reference)
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

\author{
Tobias Schoch

Maintainer: Tobias Schoch <tobias.schoch@gmail.com>
}
\references{
Anderson, E., Bai, Z., Bischof, C., Blackford, L. S., Demmel, J.,
        Dongarra, J., et al. (2000).
    \emph{LAPACK users' guide} (3rd ed.).
    Philadelphia: Society for Industrial and Applied Mathematics (SIAM).

Battese, G. E., Harter, R. M., and W.A. Fuller (1988).
    An error component model for prediction of county crop areas using.
    \emph{Journal of the American Statistical Association} 83, 28--36.

Blackford, L.S., Petitet, A., Pozo, R., Remington, K., Whaley, R.C.,
        Demmel, J., et al. (2002).
    An updated set of basic linear algebra subprograms (BLAS).
    \emph{ACM Transactions on Mathematical Software}, 28, 135--151.

Brent, R.P. (1973).
    \emph{Algorithms for minimization without derivatives}.
    Englewood Cliffs, NJ: Prentice-Hall.

Lahiri, P. (2003).
    On the impact of bootstrap in survey sampling and small area estimation.
    \emph{Statistical Science} 18, 199--210.

Hall, P. and T. Maiti (2006).
    On parametric bootstrap methods for small area prediction.
    \emph{Journal of the Royal Statistical Society. Series B}, 68, 221--238.

Heritier, S., Cantoni, E., Copt, S., and M.-P. Victoria-Feser (2009).
    \emph{Robust methods in Biostatistics}, New York: John Wiley and Sons.

Maronna, R.A., Martin, D., and V.J. Yohai (2006).
    \emph{Robust statistics: Theory and methods}. Chichester: John Wiley.

Rao, J.N.K. (2003).
    \emph{Small Area Estimation}, New York: John Wiley and Sons.

Richardson, A.M. and A.H. Welsh (1995).
    Robust restricted maximum likelihood in mixed linear model.
    \emph{Biometrics} 51, 1429--1439.

Rousseeuw, P. J. and K. Van Driessen (2006).
    Computing LTS regression for large data sets.
    \emph{Data Mining and Knowledge Discovery} 12, 29--45.

Schoch, T. (2012 under revision).
    Robust Unit-Level Small Area Estimation: A Fast Algorithm for Large
        Datasets,
    \emph{Austrian Journal of Statistics}.

Sinha, S.K. and J.N.K. Rao (2009).
    Robust small area estimation.
    \emph{Canadian Journal of Statistics} 37, 381--399.

Stahel, W. A. and A. Welsh (1997).
    Approaches to robust estimation in the simplest variance components model.
    \emph{Journal of Statistical Planning and Inference} 57, 295--319.
}
