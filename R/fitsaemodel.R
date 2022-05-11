# workhorse function
fitsaemodel <- function(method, model, ...)
{
    if (!inherits(model, "saemodel"))
        stop("Argument 'model' must be an object of class 'saemodel'\n",
             call. = FALSE)

    method <- match.arg(method, c("ml", "huberm", "tukeys"))
    if (method == "tukeys")
        stop("Tukey S-estimator not implemented yet!\n")
    else
        .fitsaemodel.huberm(method, model, ...)

}
# control function used in fitsaemodel
fitsaemodel.control <- function(niter = 40, iter = c(200, 200), acc = 1e-5,
    dec = 0, decorr = 0, init = "default", k_Inf = 20000, ...)
{
    stopifnot(all(acc > 0), all(iter > 0), length(niter) == 1, niter > 0,
        dec %in% c(0, 1), decorr %in% c(0, 1), k_Inf > 0, is.finite(k_Inf))

    # define acc
    if (length(acc) != 4)
        acc = rep(acc, 4)

    if (length(iter) != 2)
        iter = rep(iter[1], 2)

    # make them all positive
    init <- switch(match.arg(init, c("default", "lts", "s")),
        "default" = 0,
        "lts" = 1,
        "s" = 2)
    list(niter = niter, iter = iter, acc = acc, k_Inf = k_Inf, init = init,
        dec = dec, decorr = decorr, add = list(...))
}
# S3 print method
#FIXME:
print.fitsaemodel <- function (x, digits = 6, ...)
{
    saemodel <- attr(x, "saemodel")
    # check whether the model converged
    if (x$converged == 1) {
        yname <- attr(saemodel, "yname")
        xnames <- attr(saemodel, "xnames")
        areadef <- attr(saemodel, "areadef")
        # retrieve the estimating method
        method <- attr(x, "method")
        cat("ESTIMATES OF SAE-MODEL (model type B) \n")
        cat("Method: ", method$type, "\n")
        # branch: robust vs. non-robust methods
        if (length(method) > 1) {
            tuning <- method$tuning$k
            if (length(tuning) == 1)
                cat(paste("Robustness tuning constant: k = ", tuning, "\n"))
            else
                cat(paste("Robustness tuning constants: k_beta = ", tuning[1],
                    ", k_v = ", tuning[2], ", k_d = ", tuning[3],"\n",
                    sep = ""))
        }
        cat("---\nFixed effects\n")
        cat(paste("Model: ", yname, " ~ ", paste(xnames, collapse = " + "),
            sep = ""), "\n")
        cat("  Coefficients: \n")
        beta <- x$beta
        names(beta) <- xnames
        if (saemodel$intercept == 1)
            names(beta)[1] <- "(Intercept)"

        print.default(format(beta, digits = digits), print.gap = 2,
            quote = FALSE)
        cat("--- \nRandom effects \n")
        cat(paste("  Model: ~1| ", areadef, sep = ""), "\n")
        # warn if the raneff variance is almost zero
        if (x$theta[2] <= .Machine$double.eps^(1 / 4)) {
            cat("---\n")
            cat("NOTE THAT THE VARIANCE OF THE AREA-LEVEL RANDOM\n")
            cat("EFFECT IS ALMOST ZERO! DO YOU REALLY NEED THE\n")
            cat("RANDOM EFFECT? IF SO, GO AHEAD. HOWEVER, YOU\n")
            cat("SHOULD CONSIDER FITTING A (ROBUST) GLS MODEL.\n")
            cat("---\n")
        }
        theta <- sqrt(x$theta)
        # change the order of theta
        theta <- c(theta[2], theta[1])
        names(theta) <- c("(Intercept)", "Residual")
        theta <- as.matrix(theta)
        colnames(theta) <- "Std. Dev."
        print.default(format(t(theta), digits = digits), print.gap = 2,
            quote = FALSE)
        cat("--- \n")
        cat("Number of Observations: ", saemodel$n, "\n")
        cat("Number of Areas: ", saemodel$g, "\n\n")
    } else {
        #FIXME: take this case first (possibly using an early return instruction)
        # not converged
        cat("THE METHOD DID NOT CONVERGE!\n---\n")
        cat("  1) use convergence() of your fitted model to learn more\n")
        cat("  2) study the documentation using the command ?fitsaemodel\n")
        cat("  3) you may call fitsaemodel with 'init' equal to (either) 'lts'\n")
        cat("     or 's' (this works also for ML, though it may no be very efficient)\n")
        cat("  4) if it still does not converge, the last resort is to modify\n")
        cat("     'acc' and/or 'niter' (and hope and pray)\n")
    }
    invisible(x)
}
# S3 summary method
#FIXME: + digits ; see https://github.com/wch/r-source/blob/trunk/src/library/stats/R/lm.R
summary.fitsaemodel <- function (object, digits = 6, ...)
{
    # failure of convergence
    if (object$converged != 1) {
        cat("ALGORITHM FAILED TO CONVERGE! use convergence() to learn more\n")
        return()
    }

    # otherwise
    cat("ESTIMATION SUMMARY\n")
    method <- attr(object, "method")
    cat("Method: ", method$type, "\n")
    # branch: robust vs. non-robust methods
    if (length(method) > 1) {
        tuning <- method$tuning
        if (length(tuning) == 1) {
            cat(paste("Robustness tuning constant: ", names(tuning), " = ",
                as.numeric(tuning), "\n"))
        } else {
            for (i in 1:length(tuning))
                cat(tuning[i], "\n")
        }
    }
    #----------------------
    # fixed effects table
    saemodel <- attr(object, "saemodel")
    df <- saemodel$n - saemodel$g - saemodel$p + 1
    fixed <- object$beta
    stdfixed <- sqrt(diag(object$vcovbeta))
    tTable <- cbind(fixed, stdfixed, fixed / stdfixed, df, fixed)
    colnames(tTable) <- c("Value", "Std.Error", "t-value", "df", "p-value")
    rownames(tTable) <- attr(saemodel, "xnames")
    # p-value
    tTable[, 5] <- 2 * pt(-abs(tTable[, 3]), tTable[, 4])
    cat("---\nFixed effects\n")
    printCoefmat(tTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
    #----------------------
    # robustness properties (only for robust methods)
    wgt <- attr(object, "robustness")$wgt
    if (!is.null(wgt)) {
        wgt <- t(wgt) / saemodel$n
        colnames(wgt) <- c("fixeff", "residual var", "area raneff var")
        rownames(wgt) <- "sum(wgt)/n"
        cat("---\nDegree of downweighting/winsorization:\n\n")
        print.default(format(t(wgt), digits = digits), print.gap = 2,
            quote = FALSE)
    }
    invisible(tTable)
}
# S3 coef method to extract the coefficients
coef.fitsaemodel <- function(object, type = "both", ...)
{
    model <- attr(object, "saemodel")
    raneff <- t(as.matrix(object$theta))
    colnames(raneff) <- c("ResidualVar", "AreaVar")
    rownames(raneff) <- "raneff"

    fixeff <- t(as.matrix(object$beta))
    colnames(fixeff) <- colnames(model$X)
    rownames(fixeff) <- "fixeff"

    type <- match.arg(type, c("both", "raneff", "fixeff"))
    if (type == "both")
        type <- c("fixeff", "raneff")

    res <- NULL
    if ("fixeff" %in% type)
        res$fixeff <- fixeff
    if ("raneff" %in% type)
        res$raneff <- raneff
    res
}
