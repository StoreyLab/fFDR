#' estimate density on the unit interval; either 1d or 2d, using
#' local regression desity estimation
#' 
#' Used in the fpi0 and fqvalue density estimation steps
#' 
#' @param x either a vector or a 2-column matrix
#' @param transformation either probit (default), complementary log-log, or
#' identity (not recommended)
#' @param eval.points Points at which to evaluate the estimate, default x
#' @param cv Whether to use generalized cross-validation to choose the
#' nn (nearest neighbor) smoothing parameter. Currently not permitted for
#' 2d estimation
#' @param epsilon How close are allowed to come to 0 or 1; values are
#' thresholded to be in the range [epsilon, 1-epsilon]
#' @param maxk maxk argument passed to locfit
#' @param ... additional arguments to be passed to lp in locfit, used only
#' if cv=FALSE
#' 
#' @import locfit
#' 
#' @export
kernelUnitInterval = function(x, transformation="probit",
                              eval.points=x, cv=TRUE,
                              epsilon=1e-15,
                              maxk=100, grid=FALSE, ...) {
    transformation = match.arg(as.character(transformation),
                               c("ident", "cloglog", "probit"))
    trans = switch(transformation,
                   "ident"=identity,
                   "cloglog"=function(x) -log(-log(x)),
                   "probit"=qnorm)
    inv = switch(transformation,
                 "ident"=identity,
                 "cloglog"=function(x) exp(-exp(-x)),
                 "probit"=pnorm)
    dens = switch(transformation,
                     "ident"=function(x) 1,
                     "cloglog"=function(x) exp(-x) * exp(-exp(-x)),
                     "probit"=dnorm)

    # remove 0s and 1s by replacing them with the the most extreme
    # non 0/1 value. Further threshold them above a minimum
    process.vals = function(vals) {
        vals = pmax(vals, epsilon)
        vals = pmin(vals, 1 - epsilon)
    }
    x = process.vals(x)
    eval.points = process.vals(eval.points)

    s = trans(x)

    if (is.matrix(s)) {
        fitfunc = function(...) locfit(~ lp(s[, 1], s[, 2], ...), maxk=maxk)
    } else {
        fitfunc = function(...) locfit(~ lp(s, ...), maxk=maxk)
    }

    if (cv && !is.matrix(s)) {
        opt.nn = optimize(function(nn) {
            gcv(fitfunc(nn=nn))["gcv"]
            }, interval=c(0, 1))$minimum
        lfit = fitfunc(nn=opt.nn)
    }
    else {
        lfit = fitfunc(...)
    }
    
    if (grid & is.matrix(s)) {
        # set up a grid in the desired space
        gridsize = 50
        s1 = seq(min(s[, 1]), max(s[, 1]), diff(range(s[, 1])) / gridsize)
        s2 = seq(min(s[, 2]), max(s[, 2]), diff(range(s[, 2])) / gridsize)
        eval.s = as.matrix(expand.grid(s1, s2))
        eval.points = inv(eval.s)
    }
    else {
        # evaluate at the given points
        eval.s = trans(eval.points)
    }
    fs.hat = predict(lfit, newdata=eval.s)
    if (is.matrix(eval.points)) {
        corrector = apply(dens(eval.s), 1, prod)
    } else {
        corrector = dens(eval.s)
    }
    fx.hat = fs.hat / corrector
    
    if (is.matrix(x)) {
        colnames(eval.points) = c("x1", "x2")
        colnames(eval.s) = c("s1", "s2")
    }
    ret = as.data.frame(cbind(x=eval.points, fx=fx.hat, eval.s, fs=fs.hat))

    attr(ret, "lfit") = lfit

    return(ret)
}
