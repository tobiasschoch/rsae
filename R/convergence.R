convergence <- function(object)
{
    if (!inherits(object, "fitsaemodel"))
        stop("Method 'convergence()' is not applicable\n", call. = FALSE)
    saemodel <- attr(object, "saemodel")
    # retrieve the estimating method
    method <- attr(object, "method")
    cat("CONVERGENCE REPORT\n")
    # check whether the model converged
    if (object$converged != 1)
        cat("NOTE: ALGORITHM DID NOT CONVERGE!\n")

    #----------------------
    # niter and acc specification
    optim <- attr(object, "optim")
    acc <- optim$acc

    # 100 is the max value defined in estimation of "d"
    niter <- c(optim$niter, 100)
    together <- cbind(niter, acc)
    colnames(together) <- c("niter", "acc")
    rownames(together) <- c("overall loop", "fixeff", "residual var",
        "area raneff var")
    cat("---\nUser specified number of iterations (niter) and\nnumeric precision (acc):\n\n")
    print(together)

    #----------------------
    # used iters
    iters <- optim$usediter
    iters <- iters[rowSums(iters) != 0, ]
    colnames(iters) <- c("fixeff", "residual var", "area raneff var")
    rownames(iters) <- 1:NROW(iters)
    cat("---\nNumber of EE-specific iterations in each call\n(given the user-defined specs), reported for\neach of the",
        NROW(iters), "overall iterations:\n\n")
    print(iters)
}
