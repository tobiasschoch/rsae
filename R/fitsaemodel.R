fitsaemodel <- function(method, model, ...)
{
    thecall <- match.call()
    # check if model is appropriate
    if (inherits(model, "saemodel")) {
        # check if method exists
        # FIXME: match.arg()
        m <- match(method, c("ml", "huberm", "tukeys"))
        if (is.na(m))
            stop(paste("Method: ", method, " is not supported! \n", sep = ""))
        # ml and huberm
        # FIXME: switch
        if (m == 1 | m == 2)
            tmp <- .fitsaemodel.huberm(method, model, ...)
        # tukeys
        if (m == 3)
            stop("Tukey S-estimator not implemented, yet!\n")
    } else {
        stop(paste("Model", thecall[3],
            "must be an instance of 'saemodel' class.\n"))
    }
    tmp
}
# control function used in fitsaemodel
#FIXME: add tuning consant for robust initialization as an argument
fitsaemodel.control <- function(niter = 40, iter = c(200, 200), acc = 1e-5,
    dec = 0, decorr = 0, init = "default", ...)
{
    # define acc
    if (length(acc) != 4)
        acc = rep(acc, 4)

    if (length(iter) != 2)
        iter = rep(iter[1], 2)

    # implicitly check for postitivity
    #FIXME: move abs checks down to the return list
    acc = abs(acc)
    iter = abs(iter)
    niter = abs(niter[1])
    # define maxk (define ml method)
    maxk = 20000
    # machine eps
    eps <- .Machine$double.eps^(1 / 4)
    # make them all positive
    init <- switch(init, "default" = 0, "lts" = 1, "s" = 2)
    # define decomposition of the matrix-squareroot (0=SVD; 1=Cholesky)
    dec <- ifelse(dec == 0, 0, 1)
    # robustly decorrelate (center by median instead of the mean)
    decorr <- ifelse(decorr == 0, 0, 1)
    res = list(niter = niter, iter = iter, acc = acc, maxk = maxk, init = init,
        dec = dec, decorr = decorr, add = list(...))
    return(res)
}
# S3 print method
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
summary.fitsaemodel <- function (object, digits = 6, ...)
{
    saemodel <- attr(object, "saemodel")
    n <- saemodel$n
    g <- saemodel$g
    p <- saemodel$p
    fixeffnames <- attr(saemodel, "xnames")
    # retrieve the estimating method
    method <- attr(object, "method")
    # check whether the model converged
    converged <- object$converged
    if (converged != 1) {
        cat("ALGORITHM FAILED TO CONVERGE! use convergence() to learn more \n")
    } else {
        cat("ESTIMATION SUMMARY \n")
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
        df <- n - g - p + 1
        fixed <- object$beta
        stdfixed <- sqrt(diag(object$vcovbeta))
        tTable <- cbind(fixed, stdfixed, fixed / stdfixed, df, fixed)
        colnames(tTable) <- c("Value", "Std.Error", "t-value", "df", "p-value")
        rownames(tTable) <- fixeffnames
        tTable[, 5] <- 2 * pt(-abs(tTable[, 3]), tTable[, 4])
        cat("---\n")
        cat("Fixed effects\n")
        printCoefmat(tTable, digits = digits, P.values = TRUE,
            has.Pvalue = TRUE)

        #----------------------
        # robustness properties
        wgt <- attr(object, "robustness")$wgt
        # branch robust vs non-robust
        if (!is.null(wgt)) {
            wgt <- t(wgt) / n
            colnames(wgt) <- c("fixeff", "residual var", "area raneff var")
            rownames(wgt) <- "sum(wgt)/n"
            cat("---\n")
            cat("Degree of downweighting/winsorization:\n")
            cat("\n")
            print.default(format(t(wgt), digits = digits), print.gap = 2,
                quote = FALSE)
        }
    }
}
# S3 coef method to extract the coefficients
coef.fitsaemodel <- function(object, type = "both", ...)
{
    model <- attr(object, "saemodel")
    theta <- object$theta
    beta <- object$beta
    # FIXME: extract raneff and fixeff globaly; then, return the ones that
    # are required
    if (type == "both") {
        raneff <- t(as.matrix(object$theta))
        colnames(raneff) <- c("ResidualVar", "AreaVar")
        rownames(raneff) <- "raneff"
        fixeff <- t(as.matrix(object$beta))
        colnames(fixeff) <- colnames(model$X)
        rownames(fixeff) <- "fixeff"
        res <- list(fixeff = fixeff, raneff = raneff)
        print(fixeff)
        cat("\n")
        print(raneff)
    }
    if (type == "raneff") {
        res <- t(as.matrix(object$theta))
        colnames(res) <- c("ResidualVar", "AreaVar")
        rownames(res) <- "raneff"
        print(res)
    }
    if (type == "fixeff") {
        res <- t(as.matrix(object$beta))
        colnames(res) <- colnames(model$X)
        rownames(res) <- "fixeff"
        print(res)
    }
    invisible(res)
}
