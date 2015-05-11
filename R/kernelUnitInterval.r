#' estimate density on the unit interval; either 1d or 2d, using
#' local regression desity estimation
#' 
#' Used in the fpi0 and fqvalue density estimation steps
#' 
#' @param x either a vector or a 2-column matrix
#' @param transformation either probit (default), complementary log-log, or
#' identity (not recommended)
#' @param eval.points Points at which to evaluate the estimate, default x
#' @param subsample Number of points that are randomly subsampled for
#' computing the fit- useful for computational efficiency and for ensuring
#' the density estimation does not run out of memory. NULL means no the
#' fit is performed on all points
#' @param cv Whether to use generalized cross-validation to choose the
#' nn (nearest neighbor) smoothing parameter
#' @param epsilon How close values are allowed to come to 0
#' @param epsilon.max How close values are allowed to come to 1
#' @param maxk maxk argument passed to locfit
#' @param trim In one-dimensional fitting, the very edges often have high
#' variance. This parameter fixes the estimate on the intervals
#' (0, trim) and (1 - trim, 1).
#' @param ... additional arguments to be passed to lp in locfit, used only
#' if cv=FALSE
#' 
#' @import locfit
#' 
#' @export
kernelUnitInterval <- function(x, transformation = "probit",
                              eval.points = x, subsample = 1e5,
                              cv = TRUE, epsilon = 1e-15, epsilon.max = .999,
                              maxk = 100, trim = .02, ...) {
    transformation <- match.arg(as.character(transformation),
                               c("ident", "cloglog", "probit"))
    trans <- switch(transformation,
                    "ident" = identity,
                    "cloglog" = function(x) -log(-log(x)),
                    "probit" = qnorm)
    inv <- switch(transformation,
                  "ident" = identity,
                  "cloglog" = function(x) exp(-exp(-x)),
                  "probit" = pnorm)
    dens <- switch(transformation,
                   "ident" = function(x) 1,
                   "cloglog" = function(x) exp(-x) * exp(-exp(-x)),
                   "probit" = dnorm)

    # remove 0s and 1s by replacing them with the the most extreme
    # non 0/1 value. Further threshold them above a minimum
    process.vals <- function(vals) {
        vals <- pmax(vals, epsilon)
        vals <- pmin(vals, epsilon.max)
    }
    x <- process.vals(x)
    eval.points <- process.vals(eval.points)

    s <- trans(x)
    if (!is.null(subsample)) {
        if (is.matrix(s) && subsample < nrow(s)) {
            s <- s[sample(nrow(s), subsample), ]
        } else if (!is.matrix(s) && subsample < length(s)) {
            s <- s[sample(length(s), subsample)]
        }
    }
    if (is.matrix(s)) {
        fitfunc <- function(...) locfit(~ lp(s[, 1], s[, 2], ...), maxk = maxk)
    } else {
        fitfunc <- function(...) locfit(~ lp(s, ...), maxk = maxk)
    }

    if (cv) {
        # try fitting at several nearest-neighbor bandwidths
        # it's possible some will fail (esp for many hypotheses) due to
        # insufficient memory
        nns <- seq(.1, .9, .1)
        gcvs <- sapply(nns, function(nn) {
            tryCatch(gcv(fitfunc(nn = nn))["gcv"], error = function(e) NA)
            })

        if (all(is.na(gcvs))) {
            stop(paste("Fitting fails at all attempted bandwidths; try", 
                       "increasing maxk or working with a subsample"))
        }

        opt.nn <- nns[which.min(gcvs)]
        lfit <- fitfunc(nn = opt.nn)
    }
    else {
        lfit <- fitfunc(...)
    }
    
    # evaluate at the given points
    eval.s <- trans(eval.points)

    fs.hat <- predict(lfit, newdata=eval.s)
    if (is.matrix(eval.points)) {
        corrector <- apply(dens(eval.s), 1, prod)
    } else {
        corrector <- dens(eval.s)
    }
    fx.hat <- fs.hat / corrector
    
    if (is.matrix(x)) {
        colnames(eval.points) <- c("x1", "x2")
        colnames(eval.s) <- c("s1", "s2")
    }
    ret <- cbind(x = eval.points, fx = fx.hat, eval.s, fs = fs.hat) %>%
        as.data.frame() %>%
        dplyr::tbl_df()

    if (trim && !is.matrix(x)) {
        ret$fx[x < trim] <- ret %>%
            dplyr::slice(which.min(abs(x - trim))) %>%
            .$fx
        ret$fx[x > 1 - trim] <- ret %>%
            dplyr::slice(which.min(abs(x - (1 - trim)))) %>%
            .$fx
    }

    attr(ret, "lfit") = lfit

    return(ret)
}
