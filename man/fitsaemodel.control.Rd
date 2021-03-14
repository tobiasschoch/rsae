\name{fitsaemodel.control}
\alias{fitsaemodel.control}
\title{
Tuning parameters of fitsaemodel
}
\description{
This function carries global settings and parameter definitions that
are used by fitsaemodel (and its derivatives). Modifications of the
parameters can be delivered as additional arguments in the fitsaemodel call.
}
\usage{
fitsaemodel.control(niter = 40, iter = c(200, 200), acc = 1e-05,
                    dec=0, decorr=0, init="default", ...)
}
\arguments{
    \item{niter}{
        integer, defining the maximum number of outer-loop iterations
        (default: \code{niter=40})
    }
    \item{iter}{
        integer or vector of size 2, defining the maximum loops of the
        inner loops (default: \code{iter=c(200, 200)}; element 1 refers
        to beta; element 2 refers to v; note that d has an implicitly
        defined maxiter of 100 and cannot be modified)
    }
    \item{acc}{
        scalar or vector of size 4, defining the numeric tolerance used
        in the termination rule of the loops (default: \code{acc=1e-05};
        the positions of elements in the vector of size 4 are: 1=acc
        outer-loop; 2=acc inner-loop beta; 3=acc inner-loop v; 4=acc
        inner-loop d).
    }
    \item{dec}{
        type of matrix square root (decomposition); \code{dec=0} for
        SVD (default) or \code{dec=1} for Cholesky decomposition
    }
    \item{decorr}{
        in order to robustly decorrelate the residuals, one has to chose
        \code{decorr=1}  (default \code{decorr=0})
    }
    \item{init}{
        a character string; specifies by what method the main algorithm
        is initialized; by default: \code{init="default"}; alternatively,
        (and provided that the \pkg{robustbase} package is installed) one
        may choose a high-breakdown-point initial estimate: either "lts"
        (fast LTS regression) or "s" (S-estimate of regression). For
        datasets with more than 100,000 observations, the former is rather
        slow. For more details on the initializing methods see the
        documentation of \pkg{robustbase} (\code{"ltsReg"} and
        \code{"lmrob.S"}).
    }
    \item{\dots}{
        (will be used in the future)
    }
}
\details{
Caution! Modifying the default values of the parameters may result in
onvergence failure and/or loss of convergence speed.
}
\value{
(an object used by the robust methods)
}
\author{
Tobias Schoch
}
\seealso{
   \code{\link{fitsaemodel}}
}
