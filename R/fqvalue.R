#' Estimate functional q-values based on p-values and an informative factor
#' z0
#' 
#' Estimate functional q-values based on a vector of p-values and
#' z0, where z0 is known to be related to either power or the likelihood of
#' the null hypothesis.
#' 
#' @param pvalue A vector of p-values
#' @param z0 a vector providing information either on the power of each
#' test or the likelihood of the null
#' @param fixed.pi0 Whether the true pi0 (null hypothesis rate) is believed to be the same for all z.
#' @param monotone.window Parameter used to force estimated densities to
#' decrease with increasing p-values- higher means densities are smoothed
#' more aggressively. If NULL, perform no such smoothing.
#' @param transformation Transformation of p-values before the kernel density is estimated: either "cloglog" (default)=log(-log(p)) or "probit", the inverse of the cumulative density of the normal distribution.
#' @param ... Extra arguments to be passed to estFPi0
#' 
#' @return An object of class \linkS4class{fqvalue}
#' 
#' @import data.table
#' @importFrom Rcpp evalCpp
#' 
#' @useDynLib fFDR
#' 
#' @export
fqvalue = function(pvalue, z0, fixed.pi0=FALSE, monotone.window=.01, 
                   transformation="cloglog", ...) {
    dt = data.table(pvalue=pvalue, z=rank(z0) / length(z0), z0=z0)

    # calculate functional pi0
    if (fixed.pi0) {
        # assume the true pi0 doesn't actually change with Z, only power
        stop("Not currently functional")
        dt$fpi0 = fixedPi0(pvalue, dt$z, ...)
    }
    else {
        fpi0 = estFPi0(dt$pvalue, dt$z)
        dt$fpi0 = as.numeric(fpi0)
    }

    # calculate the density with a kernel estimator, on a transformed scale
    kd = kernelUnitInterval(cbind(dt$z, dt$pvalue), transformation = "cloglog")
    
    # if monotone.window is given, force density to be monotonically decreasing
    # as p-value increases
    if (!is.null(monotone.window)) {
        dt$original.fx = kd$fx
        orderer = rank(dt$pvalue)
        dt = dt[order(pvalue), ]
        # call C++ function monoSmooth to monotonically constrain
        dt[, fx:=monoSmooth(pvalue, z, original.fx, monotone.window)]
        dt = dt[orderer, ]
    } else {
        dt$fx = kd$fx
    }

    dt$lfdr = dt$fpi0 / dt$fx

    # fqvalue is the cumulative mean of the functional lfdr
    # each lFDR must be <= 1, though ordered using the unconstrained density
    dt[, fqvalue:=(cumsum(pmin(sort(lfdr), 1)) / seq_along(lfdr))[rank(lfdr)]]
    dt[, lfdr:=pmin(lfdr, 1)]

    ret = new("fqvalue", table=dt, fPi0=fpi0, density=kd)
    ret
}
