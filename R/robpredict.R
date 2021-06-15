robpredict <- function(fit, areameans = NULL, k = NULL, reps = NULL)
{
    if (!inherits(fit, "fitsaemodel"))
        stop("fit must be of class 'fitsaemodel'")
    # get the modelk
    modelk <- attr(fit, "method")$tuning$k
    # get the decomposition
    dec <- attr(fit, "dec")
    # for ml, modelk = 20000
    if (is.null(modelk))
        modelk <- 20000
    # prediction k (which is not necessarily equal to modelk)
    if (is.null(k))
        k <- modelk # here: k == modelk
    else
        if (k <= 0 )
            stop("Robustness tuning constant k must be > 0!\n")
    kappa <- .computekappa(k)
    # initialize
    mspe <- NULL
    # from the model definitions; used in order to compute the random effects
    model <- attr(fit, "saemodel")
    areaNames <- attr(model, "areaNames")
    x <- model$X
    y <- model$y
    n <- model$n
    p <- model$p
    g <- model$g
    nsize <- model$nsize
    # from the fitted model
    beta <- fit$beta
    v <- fit$theta[1]
    d <- fit$theta[2] / v
    # tau
    tau <- c(beta, v, d)
    # preparations for fortran-call
    predre <- rep(0, g)
    predfe <- rep(0, g)
    tmp <- .Fortran("drsaehubpredict", n = as.integer(n), p = as.integer(p),
        g = as.integer(g), nsize = as.integer(nsize), k = as.double(k),
        kappa = as.double(kappa), d = as.double(d), v = as.double(v),
        beta = as.matrix(beta), yvec = as.matrix(y), xmat = as.matrix(x),
        predfe = as.matrix(predfe), predre = as.matrix(predre),
        dec = as.integer(dec), PACKAGE = "rsae")
    # retrieve the area-level random effects; it is used whether new data
    # is present or not
    raneff <- tmp$predre
    # branch: old vs new data
    if (is.null(areameans)) {
        fixeff <- as.matrix(tmp$predfe)
    } else {
        # check whether the new data are proper
        if (!is.matrix(areameans))
            areameans <- as.matrix(areameans)
        # check the dimensions
        if (dim(areameans)[1] != g)
            stop("'areameans' is not of conformable size! \n")
        if (dim(areameans)[2] != p)
            stop("'areameans' is not of conformable size! \n")
        # compute the fixed-effect spredictions (at the area level)
        fixeff <- areameans %*% beta
        # compute mspe, given the fitted model (this option is only valid
        #if areameans != NULL
        if (!is.null(reps))
            mspe <- .mspe(fit, abs(round(reps)), areameans, fixeff)
    }
    means <- raneff + fixeff
    rownames(fixeff) <- areaNames
    rownames(raneff) <- areaNames
    rownames(means) <- areaNames
    # compute the residuals of the model (i.e. e_ij = y_ij - X_ij*beta - u_i)
    vn <- numeric(n)
    getres <- .Fortran("drsaeresid", n = as.integer(n), p = as.integer(p),
        g = as.integer(g), nsize = as.integer(nsize), k = as.double(modelk),
        tau = as.matrix(tau), u = as.matrix(raneff), xmat = as.matrix(x),
        yvec = as.matrix(y), res = as.matrix(vn), stdres = as.matrix(vn),
        wgt = as.matrix(vn), dec = as.integer(dec), PACKAGE = "rsae")
    result <- list(fixeff = fixeff, raneff = raneff, means = means,
        res = getres$res, stdres = getres$stdres, wgt = getres$wgt,
        mspe = mspe)
    attr(result, "robustness") <- k
    attr(result, "fit") <- fit
    attr(result, "mspe") <- reps
    class(result) <- "meanssaemodel"
    return(result)
}
# S3 plot method
plot.meanssaemodel <- function(x, y = NULL, type = "e", sort = NULL, ...)
{
    # y is part of the generic (later, we will use y to allow plot comparison)
    # extract robustness tuning constant
    k <- attr(x, "robustness")
    fe <- x$fixeff
    re <- x$raneff
    means <- x$means
    g <- length(re)
    areaNames <- rownames(re)
    mspe <- x$mspe
    # sorting
    if (!is.null(sort)) {
        ord <- switch(sort, raneff = order(re), ranef = order(re),
            fixeff = order(fe), fixef = order(fe), means = order(means))
        re <- re[ord]
        fe <- fe[ord]
        means <- means[ord]
        areaNames <- areaNames[ord]
        if (!is.null(mspe))
            mspe <- mspe[ord]
    }
    # prepare the plot
    at <- 1:g
    # without mspe
    if (is.null(mspe)) {
        ra <- range(means)
        ra[1] <- min(fe)
        # add an null-line to the plot (for the legend)
        g <- g + 1
        re <- c(re, NA)
        fe <- c(fe, NA)
        means <- c(means, NA)
        areaNames <- c(areaNames, "")
        op <- par(mfcol = c(1, 1), mar = c(4, 8, 2, 4))
        plot(means, at, type = "b", col = 2, lwd = 2, axes = FALSE,
            xlab = "predicted mean", ylab = "", xlim = ra,
            main = "Predicted means")
        lines(fe, at, type = "b", col = 1, lwd = 2, xlim = ra)
        axis(2, seq(1, g), labels = areaNames, las = 1)
        axis(1)
        grid(col = "gray65", lty = 2, lwd = 1)
        box()
        legend("top", pch = c(1, 1), lty = c(1, 1), col = c(1, 2),
            legend = c("fixeff prediction", "full prediction"),
            bg = "white", ncol = 2)
    } else {
        # with mspe
        # check type
        if (is.na(match(type, c("e", "l"))))
            stop("plot type must be either 'e' or 'l' \n")
        # compute plotting range
        ra <- range(means + sqrt(mspe), means - sqrt(mspe))
        op <- par(mfcol =c (1, 1), mar = c(8, 4, 2, 4))
        plot(at, means, type = "b", col = 2, lwd = 2, axes = FALSE,
                xlab = "", ylab = "predicted area mean", ylim = ra,
                main = "Predicted means (+/- SQRT[MSPE])")
        if (type == "l") {
            lines(at, means + sqrt(mspe), type = "b", col = 1, lwd = 2,
                lty = 2, xlim = ra)
            lines(at, means - sqrt(mspe), type = "b", col = 1, lwd = 2,
                lty = 2, xlim = ra)
        }
        if (type == "e") {
            arrows(at, means - sqrt(mspe), at, means + sqrt(mspe), angle = 90,
                code = 3, ...)
        }
        axis(1, seq(1, g), labels = abbreviate(areaNames, minlength = 12),
            las = 2)
        axis(2)
        grid(col = "gray65", lty = 2, lwd = 1)
        box()
    }
}
# S3 print method
print.meanssaemodel <- function(x, digits = 4, ...)
{
    cat("Robustly Estimated/Predicted Area-Level Means:\n")
    hasmspe <- attr(x, "mspe")
    if (is.null(hasmspe)) {
        all <- cbind(x$raneff, x$fixeff, x$means)
        colnames(all) <- c("raneff", "fixeff", "predicted mean")
    } else {
        all <- cbind(x$raneff, x$fixeff, x$means, x$mspe)
        colnames(all) <- c("raneff", "fixeff", "area mean", "MSPE")
    }
    print.default(format(all, digits = digits), print.gap = 2, quote = FALSE)
    if (!is.null(hasmspe))
        cat(paste("(MSPE: ", hasmspe," boostrap replicates)\n", sep = ""))
}
# S3 residual method to extract residuals
residuals.meanssaemodel <- function(object, ...)
{
    res <- object$res
    return(res)
}
