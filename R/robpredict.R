robpredict <- function(fit, areameans = NULL, k = NULL, reps = NULL,
        seed = 1024, progress_bar = TRUE)
{
    if (!inherits(fit, "fit_model_b"))
        stop("fit must be of class 'fit_model_b'", call. = FALSE)
    if (fit$converged == 0)
        stop("Prediction is not possible because algorithm\nof fitted model did not converge\n\n",
             call. = FALSE)

    model <- attr(fit, "saemodel")          # sae model
    dec <- attr(fit, "dec")                 # type of decomposition
    k_fit <- attr(fit, "method")$tuning$k   # tuning constant of fitted model
    if (is.null(k_fit))
        k_fit <- fitsaemodel.control()$k_Inf

    # tuning constant k for prediction (not necessarily equal to k_fit)
    if (is.null(k))
        k <- k_fit
    else
        stopifnot(is.numeric(k), k > 0, is.finite(k))
    kappa <- .computekappa(k)               # consistency correction
    v <- fit$theta[1]                       # variance components
    d <- fit$theta[2] / v

    # preparations for fortran-call
    predre <- predfe <- rep(0, model$g)
    tmp <- .Fortran("drsaehubpredict", n = as.integer(model$n),
        p = as.integer(model$p), g = as.integer(model$g),
        nsize = as.integer(model$nsize), k = as.double(k),
        kappa = as.double(kappa), d = as.double(d), v = as.double(v),
        beta = as.matrix(fit$beta), yvec = as.matrix(model$y),
        xmat = as.matrix(model$X), predfe = as.matrix(predfe),
        predre = as.matrix(predre), dec = as.integer(dec), PACKAGE = "rsae")
    # predicted area-level random effects
    raneff <- tmp$predre
    # branch: old vs new data
    mspe <- NULL
    if (is.null(areameans)) {
        fixeff <- as.matrix(tmp$predfe)
    } else {
        # check whether the new data are proper
        if (!is.matrix(areameans))
            areameans <- as.matrix(areameans)
        # check the dimensions
        if (!all(dim(areameans) == c(model$g, model$p)))
            stop("'areameans' is not of conformable size! \n")
        # compute the fixed-effect predictions (at the area level)
        fixeff <- areameans %*% fit$beta
        # compute mean square prediction error (bootstrap)
        if (!is.null(reps)) {
            stopifnot(is.numeric(reps), reps > 0, is.numeric(seed))
            set.seed(seed)
            mspe <- .mspe(fit, as.integer(reps), areameans, fixeff,
                progress_bar)
        }
    }
    rownames(fixeff) <- attr(model, "areaNames")
    rownames(raneff) <- rownames(fixeff)
    means <- raneff + fixeff

    # compute the residuals of the model (i.e. e_ij = y_ij - X_ij*beta - u_i)
    tau <- c(fit$beta, v, d)
    vn <- numeric(model$n)
    getres <- .Fortran("drsaeresid", n = as.integer(model$n),
        p = as.integer(model$p), g = as.integer(model$g),
        nsize = as.integer(model$nsize), k = as.double(k_fit),
        tau = as.matrix(tau), u = as.matrix(raneff), xmat = as.matrix(model$X),
        yvec = as.matrix(model$y), res = as.matrix(vn), stdres = as.matrix(vn),
        wgt = as.matrix(vn), dec = as.integer(dec), PACKAGE = "rsae")

    result <- list(fixeff = fixeff, raneff = raneff, means = means,
        res = getres$res, stdres = getres$stdres, wgt = getres$wgt,
        mspe = mspe)
    attr(result, "robustness") <- k
    attr(result, "fit") <- fit
    attr(result, "mspe") <- reps
    class(result) <- "pred_model_b"
    result
}
# Bootstrap
.mspe <- function(fit, reps, areameans, fixeff, progress_bar)
{
    if (progress_bar)
        p_bar <- txtProgressBar(max = reps, style = 3)

    theta <- sqrt(fit$theta)
    model <- attr(fit, "saemodel")
    Xbeta <- as.matrix(model$X) %*% fit$beta
    predicts <- matrix(NA_real_, reps, model$g)
    for (i in 1:reps) {
        # draw model error, e
        e <- rnorm(model$n, 0, theta[1])
        # draw raneff, v
        v <- unlist(sapply(model$nsize, function(u) rep(rnorm(1, 0,
            theta[2]), u), simplify = TRUE))
        # modify pred (add random effec)
        predrf <- fixeff + unique(v)
        # generate bootstrap samples (and fix it to the model)
        model$y <- Xbeta + e + v
        # compute the model parameters using ml
        tmp <- fitsaemodel("ml", model)
        # predict
        tmp <- robpredict(tmp, areameans, k = 20000, reps = NULL,
            progress_bar = FALSE)
        predicts[i, ] <- t(tmp$means - predrf)

        # update progress bar
        if (progress_bar)
            setTxtProgressBar(p_bar, i)
    }
    if (progress_bar)
        close(p_bar)

    colMeans(predicts^2)
}
# S3 plot method
plot.pred_model_b <- function(x, type = "e", sort = NULL, ...)
{
    fe <- x$fixeff
    re <- x$raneff
    means <- x$means
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
    g <- length(re)
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
        plot(means, 1:g, type = "b", col = 2, lwd = 2, axes = FALSE,
            xlab = "predicted mean", ylab = "", xlim = ra,
            main = "Predicted means")
        lines(fe, 1:g, type = "b", col = 1, lwd = 2, xlim = ra)
        axis(2, seq(1, g), labels = areaNames, las = 1)
        axis(1)
        grid(col = "gray65", lty = 2, lwd = 1)
        box()
        legend("top", pch = c(1, 1), lty = c(1, 1), col = c(1, 2),
            legend = c("fixeff prediction", "full prediction"),
            bg = "white", ncol = 2)
    } else {
        # with mspe
        if (is.na(match(type, c("e", "l"))))
            stop("plot type must be either 'e' or 'l' \n")
        # compute plotting range
        ra <- range(means + sqrt(mspe), means - sqrt(mspe))
        op <- par(mfcol =c (1, 1), mar = c(8, 4, 2, 4))
        plot(1:g, means, type = "b", col = 2, lwd = 2, axes = FALSE,
                xlab = "", ylab = "predicted area mean", ylim = ra,
                main = "Predicted means (+/- SQRT[MSPE])")
        if (type == "l") {
            lines(1:g, means + sqrt(mspe), type = "b", col = 1, lwd = 2,
                lty = 2, xlim = ra)
            lines(1:g, means - sqrt(mspe), type = "b", col = 1, lwd = 2,
                lty = 2, xlim = ra)
        }
        if (type == "e") {
            arrows(1:g, means - sqrt(mspe), 1:g, means + sqrt(mspe),
                angle = 90, code = 3, ...)
        }
        axis(1, seq(1, g), labels = abbreviate(areaNames, minlength = 12),
            las = 2)
        axis(2)
        grid(col = "gray65", lty = 2, lwd = 1)
        box()
    }
}
# S3 print method
print.pred_model_b <- function(x, digits = max(3L, getOption("digits")
    - 3L), ...)
{
    cat("Robustly Estimated/Predicted Area-Level Means:\n")

    all <- cbind(x$raneff, x$fixeff, x$means, x$mspe)
    colnames(all) <- if (is.null(x$mspe))
        c("raneff", "fixeff", "area mean")
    else
        c("raneff", "fixeff", "area mean", "MSPE")

    print.default(format(all, digits = digits), print.gap = 2, quote = FALSE)
    if (!is.null(x$mspe))
        cat("(MSPE:", attr(x, "mspe"), "boostrap replicates)\n")
}
# S3 residual method to extract residuals
residuals.pred_model_b <- function(object, ...)
{
    object$res
}
