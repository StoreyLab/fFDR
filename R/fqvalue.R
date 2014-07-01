#' Estimate functional q-values based on p-values and a power surrogate
#' 
#' Estimate functional q-values based on a vector of p-values and
#' a "power surrogate" Z, where Z is a known value that tends to increase
#' with the power of the test.
#' 
#' @param p A vector of p-values
#' @param Z a vector of the power surrogate variable, where a higher
#' Z generally represents a higher power of the hypothesis test.
#' @param fixed.pi0 Whether the true pi0 (null hypothesis rate) is believed to be the same for all Z.
#' @param lambda The value of the tuning parameter to estimate pi0.
#' @param pi0.method Method used to estimate pi0 as a function of Z:
#' either "kernel" (default), "spline", "logistic", or "LOESS". Ignored if
#' fixed.pi0 is TRUE
#' @param bandwidth The scaling bandwidth to use for the kernel density estimation: either a single value or a vector of bandwidths for smoothing in the qZ and p-value directions, respectively.
#' @param grid.points The number of grid points to use in each direction in kernel smoothing, either a single value or a vector of integers for smoothing in the qZ and p-value directions, respectively.
#' @param transformation Transformation of p-values before the kernel density is estimated: either "cloglog" (default)=log(-log(p)) or "probit", the inverse of the cumulative density of the normal distribution.
#' @param pi0.bandwidth If pi0.method is "kernel" (ignored otherwise), the scaling bandwidth to use for the kernel estimation of pi0
#' @param smooth.df If pi0.method is "spline" (ignored otherwise),
#' the number of degrees of freedom to use. If NULL (default), use
#' cross validation to decide
#' @param only.pi0 Return as soon as (f)pi0 is calculated, without calculating
#' fqvalue. Useful in simulation
#' @param ... Extra arguments to be passed to either
#' 
#' @return A data.table with columns:
#' 
#' \item{pvalue}{The original pvalues given to the function (in the same order)}
#' \item{Z}{The original Z given to the function}
#' \item{qZ}{The quantiles of Z, which are used for the nonparametric calculations of pi0 and fqvalue.}
#' \item{pi0}{Functional pi0, which depends on qZ.}
#' \item{fqvalue}{Functional q-value: the estimated false discovery rate of rejecting all hypotheses with this f-qvalue or lower}
#' 
#' @import data.table
#' @import MASS
#' @import Iso
#' @import fields
#' 
#' @export
fqvalue = function(pvalue, z0, oracle=NULL, fixed.pi0=FALSE,
                   p.bandwidth=NULL, z.bandwidth=NULL,
                   grid.points=250, transformation="cloglog",
                   min.pval=1e-30, ...) {
    # calculate functional pi0
    dt = data.table(pvalue=pvalue, z=rank(z0) / length(z0), z0=z0)
    
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

    # rows of k$z are each value of qZ, columns are each p-value
    # calculate local FDR grid
    #pi0qZ = approx(dt$qZ, dt$pi0, k$x)$y
    # lFDR is pi0 over the density
    dt$lfdr = dt$fpi0 / kd$fx

    # ensure lFDR is monotonically increasing with increasing p-values
    # using isotonic regression
    #iso.lfdr = t(apply(lfdr.grid, 1, function(z) pava(z)))
    #iso.lfdr = lfdr.grid

    # match hypotheses back to the grid with bilinear interpolation
    #iso.k = list(x=k$x, y=k$y, z=iso.lfdr)
    #dt$lfdr = interp.surface(iso.k, cbind(dt$qZ, dt$transformed.pvalue))

    # fqvalue is the cumulative mean of the functional lfdr
    # each lFDR must be <= 1, though ordered using the unconstrained density
    dt[, fqvalue:=(cumsum(pmin(sort(lfdr), 1)) / seq_along(lfdr))[rank(lfdr)]]
    dt[, lfdr:=pmin(lfdr, 1)]

    # if there was an oracle column in the input, save it
    if (!is.null(oracle)) {
        dt$oracle = oracle
    }
    ret = new("fqvalue", table=dt, fPi0=fpi0, density=kd)
    ret
}
