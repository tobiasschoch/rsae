%FIXME:
\name{robpredict}
\alias{robpredict}
\alias{print.meanssaemodel}
\alias{plot.meanssaemodel}
\alias{residuals.meanssaemodel}
\title{
Robust Prediction of Random Effects, Fixed Effects, and Area-Specific Means
}
\description{
The function \code{robpredict} robustly predicts the random effects,
fixed effects, and area-specific means under the model. As concerned with
robustly predicting the realizations of the random effects, we rely on
the method of Copt and Victoria-Feser (cf. Heritier et al., 2009, 113-114);
not the method of Sinha and Rao (2009).
}
\usage{
robpredict(fit, areameans = NULL, k = NULL, reps = NULL)

\method{print}{meanssaemodel}(x, digits = max(3L, getOption("digits") - 3L),
    ...)
\method{plot}{meanssaemodel}(x, y = NULL, type = "e", sort = NULL, ...)
\method{residuals}{meanssaemodel}(object, ...)
}
\arguments{
    \item{fit}{a fitted SAE model; object of class \code{fitsaemodel}.}
    \item{areameans}{numeric matrix (typically, with area-level means);
        the no. of rows must be equal to the no. of areas; the no. of
        columns must be equal to the no. of fixed-effects coefficients
        (incl. intercept). By default: \code{areadata = NULL}, i.e.,
        predictions are based on those data that have been used to
        estimate the model.}
    \item{k}{robustness tuning constant (of the Huber psi-function) for robust
        prediction. Notice that \code{k} does not necessarily be the same
        as the \code{k} that has been used in \code{fitsaemodel}. By
        default, \code{k} is equal to the tuning constant used in
        estimating the model parameters.}
    \item{reps}{number (integer) of bootstrap replicates for mean squared
        prediction error; default: \code{reps = NULL}.}
    \item{x}{object of the class \code{"meanssaemodel"}; this argument is
        only used in the \code{print} method.}
    \item{digits}{integer, defining the number of decimal places to be
        shown in the \code{print} method (default: \code{digits = 4}).}
    \item{y}{has no meaning, yet! (default: \code{y = NULL}; needs to
        included in the args list, because it is part of plot's generic
        arg definition).}
    \item{type}{character specifying the \code{plot} method; either \code{"e"}
        (error bars; default) or \code{"l"} (lines).}
    \item{sort}{only used in the \code{plot} method; if \code{sort = "means"},
        the predicted means are ploted in ascending order (default:
        \code{sort = NULL}); similarly, with \code{sort = "fixef"} and
        \code{sort = "ranef"} the predicted means are sorted along the
        fixed effects or the random effects, respectively.}
    \item{object}{object of the class \code{fitsaemodel}; a fitted model
        used in the \code{residuals} method.}
    \item{\dots}{not used.}
}
\details{
The \code{robpredict} function enables the following modes of prediction:

\itemize{
    \item if \code{areameans = NULL}, then the predictions are exclusively
        based on the sample values,
    \item if \code{robpredict} is called with \code{areameans} (i.e.,
        matrix with area-specific means of the auxiliary data of
        conformable size), then the fixed-effect predictions and thus
        also the predictions of the area-specific means are based on
        the auxiliary data,
    \item if, in addition to specifying \code{areameans}, one specifies
        also the number of bootstrap replications (i.e., \code{reps};
        some positive integer), the function computes area-specific mean
        square prediction error (MSPE) estimates for the area-level means.
        The MSPE is obtained, in line with Sinha and Rao (2009), from a
        (robust) parametric bootstrap; see Lahiri (2003) and Hall and
        Maiti (2006) for more details.
    }

The tuning constant \code{k} regulates the degree of robustness
(i.e., degree of winsorization of the Huber psi-function) when predicting
the random effects. If \code{k} is sufficiently large (ideally, if
\code{k} is equal to infinity), the predictions correspond to the EBLUP.
}
\value{
Instance of the S3 class \code{meanssaemodel}
}
\references{
Copt, S. and M.-P. Victoria-Feser (2009).
    \emph{Robust Predictions in Mixed Linear Models},
    Research Report, University of Geneva.

Lahiri, P. (2003).
    On the impact of bootstrap in survey sampling and small area estimation.
    \emph{Statistical Science} 18, 199-210.
    \doi{https://doi.org/10.1214/ss/1063994975}

Hall, P. and T. Maiti (2006).
    On parametric bootstrap methods for small area prediction.
    \emph{Journal of the Royal Statistical Society. Series B}, 68, 221-238.
    \doi{https://doi.org/10.1111/j.1467-9868.2006.00541.x}

Heritier, S., Cantoni, E., Copt, S., and M.-P. Victoria-Feser (2009).
    \emph{Robust methods in biostatistics}. New York: John Wiley and Sons.

Schoch, T. (2012). Robust Unit-Level Small Area Estimation: A Fast Algorithm
    for Large Datasets. \emph{Austrian Journal of Statistics} 41 (4), 243-265.
    \url{http://www.stat.tugraz.at/AJS/ausg124/124Schoch.pdf}

Sinha, S.K. and J.N.K. Rao (2009).
    Robust small area estimation.
    \emph{Canadian Journal of Statistics} 37, 381-399.
    \doi{https://doi.org/10.1002/cjs.10029}

}
\examples{
# generate synthetic data
mymodel <- makedata()
# Huber M-estimate with robustness tuning constant k = 2
myfittedmodel <- fitsaemodel("huberm", mymodel, k = 2)
myfittedmodel
# summarize the fitted model
summary(myfittedmodel)
# robust prediction of the random effects and the area-level means; here,
# we choose the robustness tuning constant k equal to 1.8
mypredictions <- robpredict(myfittedmodel, k = 1.8)
mypredictions
}
