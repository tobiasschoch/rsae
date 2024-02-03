\name{robpredict}
\alias{robpredict}
\alias{print.pred_model_b}
\alias{plot.pred_model_b}
\alias{residuals.pred_model_b}
\alias{as.matrix.pred_model_b}
\alias{head.pred_model_b}
\alias{tail.pred_model_b}

\title{Robust Prediction of Random Effects, Fixed Effects, and
    Area-Specific Means}
\description{
    Function \code{robpredict()} predicts the area-level means by (1) the
    empirical best linear unbiased predictor (EBLUP) or (2) a robust
    prediction method which is due to Copt and Victoria-Feser (2009).
    In addition, the function computes the mean square prediction
    error (MSPE) of the predicted area-level means by a parametric
    bootstrap method.
}
\usage{
robpredict(fit, areameans = NULL, k = NULL, reps = NULL, seed = 1024,
    progress_bar = TRUE)

\method{print}{pred_model_b}(x, digits = max(3L, getOption("digits") - 3L),
    ...)
\method{plot}{pred_model_b}(x, type = "e", sort = NULL, ...)
\method{residuals}{pred_model_b}(object, ...)
\method{as.matrix}{pred_model_b}(x, ...)
\method{head}{pred_model_b}(x, n = 6L, ...)
\method{tail}{pred_model_b}(x, n = 6L, keepnums = TRUE, addrownums, ...)
}
\arguments{
    \item{fit}{an object of class \code{fit_model_b}; a fitted SAE model.}
    \item{areameans}{\code{[matrix}]  or \code{[NULL]} area-level means of
        dimension \code{(g, p)}, where \code{g} and \code{p} denote,
        respectively, the number of areas and number of fixed-effects terms
        in the regression model (incl. intercept). By default,
        \code{areadata = NULL} which implies that the predictions are
        base on the data used in fitting the model (not new data).}
    \item{k}{\code{[numeric]} or \code{[NULL]} robustness tuning constant
        (of the Huber psi-function) for robust prediction. Note that
        \code{k} does not necessarily have to be the same as the \code{k}
        that has been used in \code{\link[=fitsaemodel]{fitsaemodel()}}.
        By default, \code{k = NULL} which implies that the tuning constant
        specified in \code{fitsaemodel()} is used.}
    \item{reps}{\code{[integer}] or \code{[NULL]} number of bootstrap
        replicates for the computation of the mean squared prediction
        error (MSPE).  If \code{reps = NULL} the MSPE is not computed.}
    \item{seed}{\code{[integer]} a positive integer used as argument
        \code{seed} in \code{\link[=set.seed]{set.seed()}} to specify
        the random seed.}
    \item{progress_bar}{\code{[logical]} whether a progress bar is displayed
        for the parametric bootstrap; see \bold{NOTE} below.}
    \item{x}{an object of class \code{"pred_model_b"}.}
    \item{digits}{\code{[integer]} number of digits to be printed by.}
    \item{type}{\code{[character]} type of \code{plot} method: \code{"e"}
        (error bars; default) or \code{"l"} (lines).}
    \item{sort}{\code{[character]} or \code{[NULL]} if \code{sort = "means"},
        the predicted means are plotted in ascending order (default:
        \code{sort = NULL}); similarly, with \code{sort = "fixef"} and
        \code{sort = "ranef"} the predicted means are sorted along the
        fixed effects or the random effects, respectively.}
    \item{object}{an object of class \code{fit_model_b}.}
    \item{n}{\code{[integer]} vector of length up to \code{dim(x)}, i.e.,
        number of areas.}
    \item{keepnums}{in each dimension, if no names in that dimension are
        present, create them using the indices included in that dimension.
        Ignored if \code{dim(x)} is \code{NULL} or its length 1.}
    \item{addrownums}{deprecated - \code{keepnums} should be used instead.}
    \item{\dots}{additional arguments (not used).}
}
\details{
    Function \code{robpredict()} computes predictions of the area-level means
    and---if required---an estimate of the area-specific mean square
    prediction error (MSPE).
    \describe{
        \item{Prediction of area-level means}{
            \itemize{
                \item Case 1: If \code{areameans} is a \code{matrix} with
                    area-level means of the explanatory variables, then
                    the computation of the fixed effects effects are
                    based on \code{areameans}.
                \item Case 2: If \code{areameans = NULL}, then the predictions
                    are based on the sample data that have been used to
                    fit the model.
            }
        }
        \item{Mean square prediction error}{
            \itemize{
                \item If \code{reps = NULL}, the number of bootstrap
                    replicates is not specified; hence, MSPE is not computed.
                \item If \code{reps} is a positive integer and
                \code{areameans} is not \code{NULL} (see Case 1 above), then
                a (robust) parametric bootstrap estimate of MSPE is computed
                as proposed by Sinha and Rao (2009); see also Lahiri (2003)
                and Hall.
            }
        }
        \item{Robustness}{
            \itemize{
                \item The EBLUP obtains if \code{k = NULL}, i.e., if the
                    robustness tuning constant \code{k} is unspecified.
                \item Robust predictions of the area-level means are
                    computed if \code{k} is a nonnegative real number.
                    Small values of \code{k} imply that outliers are
                    heavily downweighted; formally, the EBLUP corresponds
                    to choosing the tuning constant \code{k} equal to
                    infinity. The value of the tuning constant \code{k}
                    specified in \code{robpredict()} can be different
                    from the tuning constant \code{k} used in fitting
                    the model. The robust prediction method is due to Copt
                    and Victoria-Feser (2009); see also Heritier et al.
                    (2009, 113-114) and differs from the method in Sinha
                    and Rao (2009).
            }
        }
    }
}
\section{NOTE}{
    Users of \code{Rgui.exe} on Windows are recommended to call
    \code{robpredict()} with argument \code{progress_bar = FALSE}
    because \code{Rgui.exe} does not handle calls to
    \code{\link[=txtProgressBar]{txtProgressBar()}} well (the
    execution time of the same job increases and it tends to stall the
    execution of R). Users of \code{R-Studio} and \code{Rterm.exe}
    are \bold{not} affected.
}
\value{
    An instance of the S3 class \code{pred_model_b}
}
\references{
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
    \emph{Robust methods in biostatistics}. New York: John Wiley and Sons.

Schoch, T. (2012). Robust Unit-Level Small Area Estimation: A Fast Algorithm
    for Large Datasets. \emph{Austrian Journal of Statistics} \bold{41},
    243--265. \doi{https://doi.org/10.17713/ajs.v41i4.1548}

Sinha, S.K. and J.N.K. Rao (2009).
    Robust small area estimation.
    \emph{Canadian Journal of Statistics} \bold{37}, 381--399.
    \doi{https://doi.org/10.1002/cjs.10029}
}
\seealso{
    \code{\link[=saemodel]{saemodel()}}, \code{\link[=makedata]{makedata()}},
    \code{\link[=fitsaemodel]{fitsaemodel()}}
}
\examples{
# use the landsat data
head(landsat)

# set up the model
model <- saemodel(formula = HACorn ~ PixelsCorn + PixelsSoybeans,
    area = ~CountyName,
    data = subset(landsat, subset = (outlier == FALSE)))

# summary of the model
summary(model)

# Huber M-estimate with robustness tuning constant k = 2
m <- fitsaemodel("huberm", model, k = 2)
m

# summary of the fitted model/ estimates
summary(m)

# robust prediction of the random effects and the area-level means (robust
# EBLUP) using the counts-specific means (landsat_means)
head(landsat_means)

# for robust prediction, we use the robustness tuning constant 'k = 1.8'
m_predicted <- robpredict(m, landsat_means, k = 1.8)
head(m_predicted)

# extract prediction as matrix
as.matrix(m_predicted)

# extract residuals from the predictions
head(residuals(m_predicted))

# prediction incl. MSPE; parametric bootstrap with only 'reps = 10'
# replications (for demonstration purposes; in practice, 'reps' should be
# considerably larger)
m_predicted_mspe <- robpredict(m, landsat_means, k = 1.8, reps = 10,
                               progress_bar = FALSE)
head(m_predicted_mspe)
}
