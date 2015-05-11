#' Estimate functional q-values based on p-values and an informative factor
#' z0
#' 
#' Estimate functional q-values based on a vector of p-values and
#' z0, where z0 is known to be related to either power or the likelihood of
#' the null hypothesis.
#' 
#' @param p.value A vector of p-values
#' @param z0 a vector providing information either on the power of each
#' test or the likelihood of the null
#' @param pi0.method Method for estimating pi0(z): either "gam" (default)
#' @param lambda Lambda parameter, or choices of parameter, to use in estimating pi0(z)
#' @param fixed.pi0 Whether the true pi0 (null hypothesis rate) is believed to be
#' the same for all z.
#' @param monotone.window Parameter used to force estimated densities to
#' decrease with increasing p-values- higher means densities are smoothed
#' more aggressively. If NULL, perform no such smoothing.
#' @param ... Extra arguments to be passed to kernelUnitInterval for estimating the density
#' 
#' @return An object of S3 class "fqvalue"
#' 
#' @import data.table
#' @importFrom Rcpp evalCpp
#' 
#' @useDynLib fFDR
#' 
#' @export
fqvalue = function(p.value, z0, pi0.method = "gam", lambda = seq(.4, .9, .1),
                   fixed.pi0 = FALSE, monotone.window = .01, ...) {
    dt = data.table(p.value = p.value, z = rank(z0) / length(z0), z0 = z0)

    # calculate functional pi0
    if (fixed.pi0 == TRUE) {
        # assume the true pi0 doesn't actually change with Z, only power
        print("Fixed pi0")
        dt$fpi0 = fixed_pi0(p.value, dt$z, ...)
        fpi0 = NULL
    }
    else if (fixed.pi0 == FALSE) {
        fpi0 = estFPi0(dt$p.value, dt$z, method = pi0.method, ...)
        dt$fpi0 = as.numeric(fpi0)
    } else {
        # given a specific pi0 value to use in estimation
        dt$fpi0 = fixed.pi0
        fpi0 = NULL
    }

    # calculate the density with a kernel estimator, on a transformed scale
    kd = kernelUnitInterval(cbind(dt$z, dt$p.value), ...)
    
    # if monotone.window is given, force density to be monotonically decreasing
    # as p-value increases
    if (!is.null(monotone.window)) {
        dt$original.fx = kd$fx
        orderer = rank(dt$p.value)
        dt = dt[order(p.value), ]
        # call C++ function monoSmooth to monotonically constrain
        dt[, fx := monoSmooth(p.value, z, original.fx, monotone.window)]
        dt = dt[orderer, ]
    } else {
        dt$fx = kd$fx
    }

    dt$lfdr = dt$fpi0 / dt$fx

    # fqvalue is the cumulative mean of the functional lfdr
    # each lFDR must be <= 1, though ordered using the unconstrained density
    dt[, fq.value := (cumsum(pmin(sort(lfdr), 1)) / seq_along(lfdr))[rank(lfdr)]]
    dt[, lfdr := pmin(lfdr, 1)]

    
    ret = list(table = dt, fPi0 = fpi0, density = kd)
    class(ret) <- "fqvalue"
    ret
}


#' Convert an fqvalue object to a vector of fq-values
#' 
#' @param x "fqvalue" object
#' @param ... extra arguments, not used
#' 
#' @export
as.double.fqvalue <- function(x, ...) {
    x$table$fq.value
}


#' Convert an fqvalue object to a data frame
#' 
#' @param x "fqvalue" object
#' @param ... extra arguments, not used
#' 
#' @export
as.data.frame.fqvalue <- function(x, ...) {
    x$table
}

#' Print an fqvalue object
#' 
#' @param x "fqvalue" object
#' @param ... Extra arguments, not used
#' 
#' @export
print.fqvalue <- function(x, ...) {
    cat("Estimation of functional qvalue on", nrow(x$table), "pvalues\n")
    cat("Use plot() on this object to construct a z vs pvalue scatterplot.",
        "Use as.numeric() on this object to access the vector of fqvalues,",
        "or as.data.frame() to extract a table comparing",
        "p-values, z, and fqvalue\n")
    
    print(x$table)
}
