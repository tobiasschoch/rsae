saemodel <- function(formula, area, data, type = "b", na.omit = FALSE)
{
    # get all variables in the formula object
    variables <- all.vars(formula)
    if (!is.language(area))
        stop("area-specific random effect must be defined as formula \n")
    # Rao's model type "b" (basic unit-level model)
    if (toupper(type) == "B") {
        # get the variable(s) that define the area-specifc random effect
        areavariable <- all.vars(area)
        if (!is.na(match(areavariable, variables)))
            stop("NOTE: '", areavariable,
                "' can be either in the formula or defining the area-specific random effect, but NOT both \n")
        # limit the number of area-defining variables to one
        if (length(areavariable) != 1)
            stop("area-specific random effect not properly specified \n")
        # generate the model frame (in doing so, we can specify the behavior
        # in the case of NAs)
        na.action <- ifelse(na.omit == TRUE, "na.omit", "na.fail")
        mf <- model.frame(formula, data, na.action = na.action)
        # extract the terms
        mt <- attr(mf, "terms")
        # register if it has an intercept
        hasintercept <- attr(mt, "intercept")
        # extract the response
        y <- as.numeric(model.response(mf, "numeric"))
        # extract the design matrix
        X <- model.matrix(mt, mf)
        # dim of X
        p <- dim(X)[2]
        # check if X has full rank, else: stop
        qr <- qr(X)
        if (qr$rank < p)
            stop("Rank(design matrix) < p! Choose another model!\n")
        # generate area identifiers and sort the data
        areaID <- data[, areavariable]
        # reshape areaID into a factor
        if (!is.factor(areaID))
            areaID <- as.factor(areaID)
        # order the data along areaID so that within-area units form blocks
        ord <- order(areaID)
        mod <- data.frame(y, X, areaID)
        mod <- mod[ord, ]
        # get the area names
        areaNames <- levels(areaID)
        # vector of area size
        nsize <- unname(table(areaID))
        model <- list(X = mod[, 2:(p + 1)], y = mod[, 1], nsize = nsize,
            areaID = mod[, (p + 2)], g = length(nsize), p = p, n = length(y),
            intercept = hasintercept)
        attr(model, "areaNames") <- areaNames
        # additional stuff
        pf <- paste(formula)
        areadef <- paste(area[[2]])
        xnames <- colnames(X)
        yname <- pf[[2]]
    }
    # Rao's model type a (Fay-Herriot model, known variances must be given
    # as argument of area)

    # FIXME: function signature type = c("b", "a"); then match.arg() and switch

    if (toupper(type) == "A")
        .NotYetImplemented()
    # build the model
    attr(model, "yname") <- yname
    attr(model, "xnames") <- xnames
    attr(model, "areadef") <- areadef
    attr(model, "call") <- match.call()
    class(model) <- "saemodel"
    return(model)
}
# S3 print method
print.saemodel <- function(x, ...)
{
    # distinguish synthetic (potentially contaminated) from ordinary data

    #FIXME: this distinction must be done at the stage of model generation
    #       makedata and saemodel, respectively; the print method should only
    #       print

    if (is.null(attr(x, "contam"))) {
        cat("SAE MODEL TYPE: B (J.N.K. Rao's classification)\n---\n")
        cat(paste("FIXED EFFECTS: ", attr(x, "yname"), " ~ ", paste(attr(x,
            "xnames"), collapse = " + "), sep = ""), "\n")
        cat(paste("AREA-SPECIFIC RANDOM EFFECTS: ", attr(x, "areadef"),
            sep = ""), "\n")
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
    cat("No. of areas: ", length(nsize), "\n")
    cat("No. of obs.: ", sum(nsize), "\n")
    if (length(unique(nsize)) == 1) {
        cat("Balanced data, each area has ", nsize[1], " units \n")
    } else {
        cat("Smallest area: ", min(nsize), " units \n")
        cat("Largest area: ", max(nsize), " units \n")
    }
}
# S3 method to extract the data as matrix
as.matrix.saemodel <- function(x, ...)
{
    as.matrix(cbind(x$y, x$X, x$areaID))
}
