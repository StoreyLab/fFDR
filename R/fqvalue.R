#' Estimate functional q-value based on p-value and an informative variable
#' 
#' Estimate functional q-values based on p-values and realizations of
#' the informative variable z, where z may affect either the power of a statistical test or the likelihood of a true null hypothesis. 
#' 
#' @param p.value A vector of p-values
#' @param z0 A vector of observations from the informative variable, of the same length as \code{p}
#' @param pi0.method Method for estimating the functional proportion pi0(z); either "gam" (default), "glm", "kernel" or "bin"
#' @param lambda The parameter Lambda or its values to use in estimating pi0(z)
#' @param fixed.pi0 Whether pi0(z) is believed to be independent of z
#' @param monotone.window Parameter used to force estimated densities to
#' decrease with increasing p-values; higher value of this parameter means densities are smoothed
#' more aggressively. If NULL, perform no such smoothing
#' @param ... Extra arguments to be passed to kernelUnitInterval for estimating the density
#' 
#' @details Assume the random variable z0 may affect the power of a statistical test (that induces the p-values) or the likelihood of a true null hypothesis. The m observations z0_i, i=1,...,m of z0 are quantile transformed into z_i, i=1,...,m such that z_i = rank(z0_i) / m, where rank(z0_i) is the rank of z0_i among z0_i, i=1,...,m. Consequently, z_i, i=1,...,m are approximately uniformly distributed on the interval [0,1]. When z_i, i=1,...,m are regarded as observations from the random variable z, then z is approximately uniformly distributed on [0,1]. Namely, z0 has been quantile transformed into z, and they are equivalent. Further, z or z0 is referred to as the informative variable.
#' 
#' The likelihood of a true null hypothesis is regarded as a function of z, referrred to as the functional proportion, and denoted by pi0(z), and the fFDR methodology is applied to the m paired observations (p_i,z_i), i=1,...,m of the p-value p and z. 
#' 
#' @return An object of S3 class "fqvalue"
#' 
#' @importFrom dplyr %>%
#' @importFrom Rcpp evalCpp
#' @import broom
#' @importFrom graphics plot
#' @importFrom stats approx binom.test binomial cor dnorm  fitted.values glm lm optimize plogis pnorm predict qnorm quantile rbinom rlnorm rnorm runif t.test
#' 
#' @useDynLib fFDR
#' 
#' @export
fqvalue = function(p.value, z0, pi0.method = "gam", lambda = seq(.4, .9, .1),
                   fixed.pi0 = FALSE, monotone.window = .01, ...) {
    dt <- tibble::tibble(p.value = p.value, z = rank(z0) / length(z0), z0 = z0)

    # calculate functional pi0
    if (fixed.pi0 == TRUE) {
        # assume the true pi0 doesn't actually change with Z, only power
        dt$fpi0 <- estimate_fixed_pi0(p.value, dt$z, ...)
        fpi0 <- NULL
    }
    else if (fixed.pi0 == FALSE) {
        fpi0 <- estimate_fpi0(dt$p.value, dt$z, method = pi0.method, ...)
        dt$fpi0 <- as.numeric(fpi0)
    } else {
        # given a specific pi0 value to use in estimation
        dt$fpi0 <- fixed.pi0
        fpi0 <- NULL
    }

    # calculate the density with a kernel estimator, on a transformed scale
    kd <- kernelUnitInterval(cbind(dt$z, dt$p.value), ...)
    
    # if monotone.window is given, force density to be monotonically decreasing
    # as p-value increases
    if (!is.null(monotone.window)) {
        dt$original.fx <- kd$fx
        orderer <- rank(dt$p.value)
        # call C++ function monoSmooth to monotonically constrain
        dt <- dt %>%
            dplyr::arrange(p.value) %>%
            dplyr::mutate(fx = monoSmooth(p.value, z, original.fx, monotone.window)) %>%
            dplyr::slice(orderer)
    } else {
        dt$fx <- kd$fx
    }

    dt$lfdr <- dt$fpi0 / dt$fx

    # fqvalue is the cumulative mean of the functional lfdr
    # each lFDR must be <= 1, though ordered using the unconstrained density
    dt <- dt %>%
        mutate(fq.value = (cumsum(pmin(sort(lfdr), 1)) / seq_along(lfdr))[rank(lfdr)]) %>%
        mutate(lfdr = pmin(lfdr, 1))

    ret <- list(table = dt, fPi0 = fpi0, density = kd)
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
    cat("Use plot() and grid::grid.draw() on this object to construct a z vs pvalue scatterplot.",
        "Use as.numeric() on this object to access the vector of fqvalues,",
        "or as.data.frame() to extract a table comparing",
        "p-values, z, and fqvalue\n")
    
    print(x$table)
}
