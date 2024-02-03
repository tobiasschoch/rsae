# workhorse function
fitsaemodel <- function(method, model, ...)
{
    if (!inherits(model, "saemodel"))
        stop("Argument 'model' must be an object of class 'saemodel'\n",
             call. = FALSE)

    method <- match.arg(method, c("ml", "huberm", "tukeys"))

    res <- if (method == "tukeys")
        stop("Tukey S-estimator not implemented yet!\n")
    else
        .fit_model_b_huberm(method, model, ...)
    attr(res, "call") <- match.call()
    res
}
# control function used in fitsaemodel
fitsaemodel.control <- function(niter = 40, iter = c(200, 200), acc = 1e-5,
                                dec = 0, decorr = 0, init = "default",
                                k_Inf = 20000, ...)
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
        "default" = 0, "lts" = 1, "s" = 2)
    list(niter = niter, iter = iter, acc = acc, k_Inf = k_Inf, init = init,
         dec = dec, decorr = decorr, add = list(...))
}
