# function to specify a SAE model
saemodel <- function(formula, area, data, type = "b", na.omit = FALSE)
{
    if (!is.language(formula))
        stop("Argument 'formula' must be formula object\n", call. = FALSE)
    if (!is.language(area))
        stop("Argument 'area' must be formula object\n", call. = FALSE)
    # variables in 'formula' must be different from variables in 'area'
    if (length(intersect(all.vars(formula), all.vars(area))) != 0)
        stop("Variables in 'formula' and 'area' must not be the same \n",
             call. = FALSE)
    # only one variables
    if (length(all.vars(area)) != 1)
        stop("Only one variable is allowed in the specification of 'area'\n",
             call. = FALSE)
    # unit- vs. area-level model
    model <- switch(type,
        "a" = .saemodel_area(),
        "b" = .saemodel_unit(formula, area, data, na.omit),
         stop(paste0("Argument 'type = ", type,  "' is unknown\n"),
              call. = FALSE))
    attr(model, "call") <- match.call()
    class(model) <- "saemodel"
    model
}
.saemodel_area <- function()
{
    .NotYetImplemented()
}
.saemodel_unit <- function(formula, area, data, na.omit)
{
    mf <- model.frame(formula, data, na.action = ifelse(na.omit, "na.omit",
                                                        "na.fail"))
    mt <- attr(mf, "terms")
    y <- as.numeric(model.response(mf, "numeric"))
    X <- model.matrix(mt, mf)
    p <- NCOL(X)

    if (qr(X)$rank < p)                         # stop if X is rank deficient
        stop("Rank(design matrix) < p! Choose another model!\n", call. = FALSE)

    areaID <- data[, all.vars(area)]            # area identifiers
    if (!is.factor(areaID))                     # areaID as factor
        areaID <- as.factor(areaID)

    # order the data along areaID so that within-area units form blocks
    mod <- data.frame(y, X, areaID)
    mod <- mod[order(areaID), ]
    nsize <- unname(table(areaID))

    model <- list(X = mod[, 2:(p + 1)], y = mod[, 1], nsize = nsize,
                  areaID = mod[, (p + 2)], g = length(nsize), p = p,
                  n = length(y), intercept = attr(mt, "intercept"))
    colnames(model$X) <- colnames(X)
    attr(model, "areaNames") <- levels(areaID)
    attr(model, "yname") <- paste(formula)[[2]]
    attr(model, "xnames") <- colnames(X)
    attr(model, "areadef") <- paste(area[[2]])
    model
}
# S3 print method
print.saemodel <- function(x, ...)
{
    #FIXME: the print method should only print; do preparations of strings
    # in saemodel() and makedata()

    if (is.null(attr(x, "contam"))) {
        cat("SAE MODEL TYPE: B (J.N.K. Rao's classification)\n---\n")
        cat(paste0("FIXED EFFECTS: ", attr(x, "yname"), " ~ ",
                   paste0(attr(x, "xnames"), collapse = " + ")), "\n")
        cat(paste0("AREA-SPECIFIC RANDOM EFFECTS: ", attr(x, "areadef")), "\n")
    } else {
        skeleton <- attr(x, "contam")$skeleton
        cat("SAE MODEL TYPE: B (J.N.K. Rao's classification)\n")
        cat("DATA: Synthetic, simulated data \nMODEL:\n")
        beta <- skeleton$beta
        intercept <- skeleton$intercept
        p <- length(beta)
        lhs <- "   y_ij = "
        if (p > 1) {
            rhs <- "sum_k[ beta_k * x_kij] + v_i + e_ij"
            if (is.null(intercept))
                cat(paste0(lhs, rhs, "\n with\n"))
            else
                cat(paste0(lhs, "intercept + ", rhs, "\n with\n"))
            cat("   each x_kij ~ N(0, 1), k = 1, ...,", p, "\n")
        } else {
            rhs <- "beta * x_ij + v_i + e_ij"
            if (is.null(intercept))
                cat(paste0(lhs, rhs, "\n with\n"))
            else
                cat(paste0(lhs, "intercept + ", rhs, "\n with\n"))
                cat("   x_ij ~ N(0, 1)\n")
        }
        # random effects
        if(skeleton$vu.epsilon == 0)
            cat("   v_i ~  N(0, ", skeleton$vu,  ")\n", sep = "")
        else
            cat("   v_i ~ (", 1 - skeleton$vu.epsilon, ")*N(0, ", skeleton$vu,
                ") + ", skeleton$vu.epsilon, "*N(0, ", skeleton$vu.contam,
                ") \n", sep = "")

        if (skeleton$ve.epsilon == 0)
            cat("   e_ij ~ N(0, ", skeleton$ve,  ")\n", sep = "")
        else
            cat("   e_ij ~ (", 1 - skeleton$ve.epsilon, ")*N(0, ", skeleton$ve,
                ") + ", skeleton$ve.epsilon, "*N(0, ", skeleton$ve.contam,
                ") \n", sep = "")
    }
}
# S3 summary method
summary.saemodel <- function(object, ...)
{
    cat("Model summary: \n")
    print(attr(object, "call"))
    cat("---\n")
    nsize <- object$nsize
    cat("No. of areas:", length(nsize), "\n")
    cat("No. of obs.: ", sum(nsize), "\n")
    if (length(unique(nsize)) == 1) {
        cat("Balanced data, each area has", nsize[1], "units \n")
    } else {
        cat("Smallest area:", min(nsize), "units \n")
        cat("Largest area: ", max(nsize), "units \n")
    }
}
# S3 method to extract the data as matrix
as.matrix.saemodel <- function(x, ...)
{
    m <- as.matrix(cbind(x$y, x$X, x$areaID))
    colnames(m)[c(1, NCOL(m))] <- c("y", "areaID")
    m
}
