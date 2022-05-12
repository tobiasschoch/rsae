# consistency correction term
.computekappa <- function(k)
{
    2 * (k^2 * (1 - pnorm(k)) + pnorm(k) - 0.5  - k * dnorm(k))
}
# workhorse function
.fitsaemodel_huberm <- function(method, model, k,
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
