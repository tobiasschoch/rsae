# consistency correction term
.computekappa <- function(k)
{
    2 * (k^2 * (1 - pnorm(k)) + pnorm(k) - 0.5  - k * dnorm(k))
}
# workhorse function
.fit_model_b_huberm <- function(method, model, k,
    control = fitsaemodel.control(...), ...)
{
    # methods: 'ml' vs. 'huberm'
    if (method == "ml") {
        k <- rep(control$k_Inf, 3)
        kappa <- rep(1, 2)
        methodName <- list(type = "Maximum likelihood estimation")
    } else {
        if (missing(k))
            stop("Robustness tuning constant is missing! \n", call. = FALSE)
        if (is.list(k)) {
            if (any(is.na(match(names(k), c("beta", "v", "d")))))
                stop("k can be list with named entries 'beta', 'v', and 'd'\n",
                    call. = FALSE)
            k <- c(k$beta, k$v, k$d)
            stopifnot(all(is.numeric(k)), all(k) > 0)
            kappa <- c(.computekappa(k[2]), .computekappa(k[3]))
            k.report <- k
        } else {
            stopifnot(length(k) == 1, is.numeric(k), k > 0)
            kappa <- rep(.computekappa(k), 2)
            k <- rep(k, 3)
            k.report <- k[1]
        }
        methodName <- list(type = "Huber-type M-estimation", tuning =
            list(k = k.report))
    }
    # initialize the full parameter vector
    init <- if (length(control$add) == 0)
        .initmethod(model, control$init)
    else
        .initmethod(model, control$init, control$add)
    # compute estimates
    taurecord <- matrix(0, control$niter, (model$p + 2))
    eps <- .Machine$double.eps
    tmp <- .Fortran("drsaehub", n = as.integer(model$n),
        p = as.integer(model$p), g = as.integer(model$g),
        niter = as.integer(control$niter), nsize = as.integer(model$nsize),
        iter = as.integer(control$iter), iterrecord = as.matrix(matrix(0,
        control$niter, 3)), allacc = as.double(control$acc[1]),
        acc = as.matrix(control$acc[2:4]), sumwgt = as.matrix(rep(0, 3)),
        xmat = as.matrix(model$X), yvec = as.matrix(model$y),
        k = as.matrix(k), kappa = as.matrix(kappa),
        epsd = as.double(eps^(1 / 4)), tau = as.matrix(init),
        taurecord = as.matrix(taurecord), converged = as.integer(0),
        dec = as.integer(control$dec), decorr = as.integer(control$decorr),
        PACKAGE = "rsae")
    # check for cycling an choose the parameter-vector estimate whose
    # estimate of v is closer to the (robust) init (i.e., either the "lts" or
    # "s" estimate) value. This method is not supported for init=default or ml
    converged <- tmp$converged
    p <- model$p
    if (control$init > 0 & converged == 0) {
        taurecord <- tmp$taurecord
        # take the second diff, take the mean over each parameter vector,
        # and take the last that fullfiled the criterion
        u <- max(which(rowMeans(diff(taurecord, 2)) <= eps^(1 / 2)))
        # take the difference from init for at u
        uat <- abs(taurecord[u, (p + 1)] - init[(p + 1)])
        # take the difference from init for one before u
        ubefore <- abs(taurecord[(u - 1), (p + 1)] - init[(p + 1)])
        if (uat <= ubefore) {
            tau <- taurecord[u, ]
            converged <- 1
        } else {
            tau <- taurecord[(u-1), ]
            converged <- 1
        }
    } else {
        tau <- tmp$tau
    }
    if (converged != 1)
        tau <- rep(NA_real_, p + 2)

    res <- list(beta = tau[1:p], theta = c(tau[p+1], tau[p + 1] * tau[p + 2]),
        converged = converged)
    # additional attributes
    attr(res, "optim") <- list(acc = control$acc, niter = c(control$niter,
        control$iter), usediter = tmp$iterrecord, tau = tmp$taurecord,
        kappa = kappa)
    if (method == "huberm")
        attr(res, "robustness") <- list(wgt = tmp$sumwgt)
    attr(res, "init") <- init
    attr(res, "method") <- methodName
    attr(res, "saemodel") <- model
    attr(res, "dec") <- control$dec
    class(res) <- "fit_model_b"
    res
}
# initialization
.initmethod <- function(model, init, ...)
{
    n <- model$n
    p <- model$p
    intercept <- model$intercept
    #FIXME: use a switch?
    #-------------
    # default (i.e., robust fixed-effects estimator; see AJS2012)
    if (init == 0) {
        # FIXME: make it an argument of .initmethod
        k <- 1.345
        # retrieve all the model characteristics
        y <- model$y
        X <- as.data.frame(model$X)
        g <- model$g
        areaID <- model$areaID
        # center y by the area-specific median of y
        y.list <- split(y, areaID)
        y.centered.list <- lapply(y.list, function(u) u - median(u))
        y.centered <- unsplit(y.centered.list, areaID)
        # center X by the area-specific mean of x
        X.list <- split(X, areaID)
        X.centered.list <- lapply(X.list, function(u)
            as.data.frame(sweep(as.matrix(u), 2, colMeans(u))))
        X.centered <- unsplit(X.centered.list, areaID)
        if (intercept == 1) {
            X.centered <- X.centered[, -1]
            p <- p - 1
        }
        # prepare the model.frame
        mm <- model.matrix(~ -1 + as.factor(areaID))
        x <- cbind(X.centered, mm)
        # compute the robust fixed-effects estimator
        initbeta <- rep(1, (p + g))
        #FIXME: remove magic numbers niter = 20 and s = 1.2, acc = 0.000001
        tmp <- .Fortran("drlm", n = as.integer(n), p = as.integer(p + g),
            xmat = as.matrix(x), yvec = as.matrix(y.centered),
            k = as.double(k), beta = as.matrix(initbeta), s = as.double(1.2),
            info = as.integer(1), niter = as.integer(20),
            acc = as.double(0.00001), PACKAGE = "rsae")
        result <- c(0, tmp$beta[1:p], tmp$s^2, 100)
    }
    #-------------
    # check whether robustbase must be loaded
    if (init > 0) {
        check <- requireNamespace("robustbase", quietly = TRUE)
        if (!check)
            stop("You cannot use 'lts' or 's', because the \n
                'robustbase' package is not installed! \n", call. = FALSE)
    }
    #FIXME: use a switch?
    #-------------
    # lts
    if (init == 1) {
        x <- as.matrix(model$X)
        # check if it has an intercept
        if (intercept == 1) {
            x <- as.matrix(x[, (2:p)])
            intercept <- TRUE
        } else {
            cat
            intercept <- FALSE
        }
        y <- model$y
        tmp <- robustbase::ltsReg(x = x, y = y, intercept=intercept, ...)
        # compute (robust) variance bound of d
        # repare return value (beta, v, d)
        result <- as.numeric(c(tmp$coefficients, tmp$raw.scale^2, 1))
    }
    #-------------
    # lmrob.S
    if (init == 2) {
        x <- as.matrix(model$X)
        y <- model$y
        control <- robustbase::lmrob.control(...)
        tmp <- robustbase::lmrob.S(x = x, y = y, control=control)
        result <- as.numeric(c(tmp$coefficients, tmp$scale^2, 1))
    }
    return(result)
}
# S3 print method
print.fit_model_b <- function (x, digits = max(3L, getOption("digits") - 3L),
    ...)
{
    # failure of convergence
    if (x$converged != 1) {
        cat("THE ALGORITHM DID NOT CONVERGE!\n---\n")
        cat("  1) use convergence() of your fitted model to learn more\n")
        cat("  2) study the documentation using the command ?fitsaemodel\n")
        cat("  3) you may call fitsaemodel with 'init' equal to (either) 'lts'\n")
        cat("     or 's' (this works also for ML, though it may no be very efficient)\n")
        cat("  4) if it still does not converge, the last resort is to modify\n")
        cat("     'acc' and/or 'niter' (and hope and pray)\n\n")
    }
    # otherwise
    saemodel <- attr(x, "saemodel")
    method <- attr(x, "method")
    cat("ESTIMATES OF SAE-MODEL (model type B) \n")
    cat("Method: ", method$type, "\n")
    # branch: robust vs. non-robust methods
    if (length(method) > 1) {
        tuning <- method$tuning$k
        if (length(tuning) == 1)
            cat("Robustness tuning constant: k =", tuning, "\n")
        else
            cat(paste0("Robustness tuning constants: k_beta = ", tuning[1],
                ", k_v = ", tuning[2], ", k_d = ", tuning[3],"\n"))
    }
    cat("---\nFixed effects\n")
    cat(paste0("Model: ", attr(saemodel, "yname"), " ~ ",
        paste0(attr(saemodel, "xnames"), collapse = " + ")), "\n")
    cat("\nCoefficients: \n")
    beta <- x$beta
    names(beta) <- attr(saemodel, "xnames")
    if (saemodel$intercept == 1)
        names(beta)[1] <- "(Intercept)"

    print.default(format(beta, digits = digits), print.gap = 2, quote = FALSE)
    cat("---\nRandom effects\n")
    cat("\nModel: ~1|", attr(saemodel, "areadef"), "\n")

    # warn if the raneff variance is almost zero
    if (!is.na(x$theta[2])){
        if (x$theta[2] <= .Machine$double.eps^(1 / 4)) {
            cat("---\n")
            cat("NOTE THAT THE VARIANCE OF THE AREA-LEVEL RANDOM\n")
            cat("EFFECT IS ALMOST ZERO! DO YOU REALLY NEED THE\n")
            cat("RANDOM EFFECT? IF SO, GO AHEAD. HOWEVER, YOU\n")
            cat("SHOULD CONSIDER FITTING A (ROBUST) GLS MODEL.\n")
            cat("---\n")
        }
    }

    theta <- sqrt(x$theta)
    theta <- c(theta[2], theta[1])  # change the order of theta
    names(theta) <- c("(Intercept)", "Residual")
    theta <- as.matrix(theta)
    colnames(theta) <- "Std. Dev."
    print.default(format(t(theta), digits = digits), print.gap = 2,
        quote = FALSE)
    cat("---\nNumber of Observations: ", saemodel$n, "\nNumber of Areas: ",
        saemodel$g, "\n\n")
    invisible(x)
}
# S3 summary method
summary.fit_model_b <- function (object, ...)
{
    # covariance matrix
    model <- attr(object, "saemodel")
    vcovbeta <- if (object$converged == 1)
        .Fortran("drsaehubvariance", n = as.integer(model$n),
            p = as.integer(model$p), g = as.integer(model$g),
            nsize = as.integer(model$nsize), kappa = as.double(attr(object,
            "optim")$kappa), v = as.double(object$theta[1]),
            d = as.double(object$theta[2] / object$theta[1]),
            xmat = as.matrix(model$X), vcovbeta = as.matrix(matrix(0,
            model$p, model$p)), dec = as.integer(attr(object, "dec")),
            PACKAGE = "rsae")$vcovbeta
    else
        matrix(NA_real_, model$p, model$p)

    method <- attr(object, "method")
    res <- list(converged = object$converged, method = method,
        vcovbeta = vcovbeta)
    #----------------------
    # fixed effects table
    saemodel <- attr(object, "saemodel")
    df <- saemodel$n - saemodel$g - saemodel$p + 1
    fixed <- object$beta
    stdfixed <- sqrt(diag(vcovbeta))
    tTable <- cbind(fixed, stdfixed, fixed / stdfixed, df, fixed)
    colnames(tTable) <- c("Value", "Std.Error", "t-value", "df", "p-value")
    rownames(tTable) <- attr(saemodel, "xnames")
    # p-value
    tTable[, 5] <- 2 * pt(-abs(tTable[, 3]), tTable[, 4])
    res$tTable <- tTable
    #----------------------
    wgt <- attr(object, "robustness")$wgt
    if (!is.null(wgt) && object$converged) {
        wgt <- t(wgt) / saemodel$n
        colnames(wgt) <- c("fixeff", "residual var", "area raneff var")
        rownames(wgt) <- "sum(wgt)/n"
        res$wgt <- wgt
    }
    class(res) <- "summary_fit_model_b"
    res
}
# S3 print method for model summary
print.summary_fit_model_b <- function(x, digits = max(3L, getOption("digits")
    - 3L), ...)
{
    # failure of convergence
    if (x$converged != 1) {
        cat("ALGORITHM FAILED TO CONVERGE! use convergence() to learn more\n\n")
    }
    cat("ESTIMATION SUMMARY\n")
    cat("Method:", x$method$type, "\n")
    # branch: robust vs. non-robust methods
    if (length(x$method) > 1) {
        tuning <- x$method$tuning$k
        if (length(tuning) == 1)
            cat("Robustness tuning constant: k =", tuning, "\n")
        else
            cat(paste0("Robustness tuning constants: k_beta = ", tuning[1],
                ", k_v = ", tuning[2], ", k_d = ", tuning[3],"\n"))
    }
    cat("---\nFixed effects\n")
    printCoefmat(x$tTable, digits = digits, P.values = TRUE,
        has.Pvalue = TRUE)
    if (!is.null(x$wgt)) {
        cat("---\nDegree of downweighting/winsorization:\n\n")
        print.default(format(t(x$wgt), digits = digits), print.gap = 2,
            quote = FALSE)
    }
}
# S3 coef method to extract the coefficients
coef.fit_model_b <- function(object, type = "both", ...)
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
