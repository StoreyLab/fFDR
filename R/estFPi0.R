#' estimate pi0 as a function of z
#' 
#' \code{fpi0} estimates a pi0 function that depends on another
#' parameter z.
#' 
#' @param p a vector of p-values
#' @param z0 values on which pi0 is believed to depend
#' @param lambda Possible choices of lambda, by default {.4, .5, ..., .9}
#' @param method Either "gam" (default), "glm", "kernel", or "bin"
#' @param df Degrees of freedom to use in spline in "gam" method
#' @param breaks Either a number of (evenly spaced) break points in "bin" method,
#' or a vector of break points (from 0 to 1) to use for bins
#' @param ... Additional arguments to glm, gam or the kernel estimator
#' 
#' @details In short, the "glm", "gam", and "kernel" methods attempt
#' to estimate:
#' 
#' \code{pi_0(z) = Pr(p>lambda|z)/(1-lambda)}
#' 
#' The factor z0 that the user provides is transformed to produce z, according to
#' z = rank(z0) / length(z0). This ensures that z has an approximately uniform
#' distribution on the interval (0, 1].
#' 
#' The glm and gam approaches define a variable \code{phi=I{p>lambda}}, and
#' use a modification of logistic regression to fit \code{phi~f(z)}. The
#' kernel density estimate examines the density of z where p>lambda and computes
#' Pr(p>lambda) from that.
#' 
#' Binning simply computes the Storey estimate with the given lambda
#' within each bin.
#' 
#' @return an FPi0 object
#' 
#' @export
estFPi0 = function(p, z0, lambda=seq(.4, .9, .1), method="gam",
                   df=3, breaks=5, ...) {
    # check p-values, and assumptions of model
    if (min(p) < 0 || max(p) > 1) {
        stop("P-values not in valid range")
    }
    if (length(p) != length(z0)) {
        stop("Vectors of p-values and z0 must be same length")
    }

    # transform z0 into quantiles
    z = rank(z0) / length(z0)
    
    choices = c("glm", "gam", "kernel", "bin")
    method = match.arg(method, choices)

    pi0hat.func = function(lambda) {
        # set up a function for estimating the fpi0, which
        # returns a list containing fpi0 and the fit object
        phi = as.numeric(p > lambda)
        fit = NULL

        if (method == "glm") {
            fit = glm(phi ~ z, family=constrained.binomial(1 - lambda), ...)
            pi0 = fitted.values(fit) / (1 - lambda)
        } else if (method == "gam") {
            fit = gam::gam(phi ~ splines::ns(z, df),
                             family=constrained.binomial(1 - lambda), ...)
            pi0 = fitted.values(fit) / (1 - lambda)
        } else if (method == "kernel") {
            kd = kernelUnitInterval(z[phi == 1], transformation="probit",
                                    eval.points=z, ...)
            pi0 = kd$fx * mean(phi) / (1 - lambda)
        } else if (method == "bin") {
            if (length(breaks) == 1) {
                breaks = seq(0, 1, 1 / breaks)
            }
            b = findInterval(z, breaks, rightmost.closed=TRUE)
            freq = table(b, phi)
            b.pi0 = (freq[, 2] / rowSums(freq)) / (1 - lambda)
            pi0 = b.pi0[b]
        }
        pi0
    }
    
    # choose lambda
    fpi0s = data.table(lambda=lambda)
    fpi0s = fpi0s[, list(pvalue=p, z=z, fpi0=pi0hat.func(lambda)), by=lambda]
    
    # estimate \hat{phi} using the lowest lambda as reference 
    ref = fpi0s[lambda == min(lambda), fpi0]
    fpi0s[, k:=optimize(function(k) mean((ref-k*(1-ref)-fpi0)^2), interval=c(-1, 1))$minimum, by=lambda]
    fpi0s[, phi.hat:=ref-k*(1-ref)]

    pi0.S = qvalue(p)$pi0

    # correct for cases > 1
    fpi0s[, fpi0:=pmin(fpi0, 1)]
    fpi0s[, phi.hat:=pmin(phi.hat, 1)]
    
    stats = fpi0s[, list(omega=mean((fpi0-phi.hat)^2),
                         delta.sq=(max(mean(fpi0)-pi0.S, 0))^2),
                  by=lambda]
    stats[, MISE:=omega+delta.sq]

    # extract chosen lambda and fpi0 estimate
    lambda.hat = stats[which.min(MISE), lambda]
    fpi0s[, chosen := (lambda == lambda.hat)]
    fpi0 = fpi0s[chosen == TRUE, fpi0]

    tab = data.table(pvalue=p, z=z, z0=z0, fpi0=fpi0)
    ret = new("fPi0", table=tab, tableLambda=fpi0s, MISE=stats,
               lambda=lambda.hat, method=method)
    ret
}


#' constrained binomial family
#' 
#' This is a variation on the binomial family for GLMs and GAMs, where the
#' value of g(x) ranges from 0 to a maximum value (while for the logistic
#' link, it ranges from (0, 1).
#' 
#' @param maximum ceiling constraint on the value of g(x)
#' 
#' @return An object of class "family"
constrained.binomial = function(maximum) {
    link = structure(list(name = paste0("constrained.logit (0, ", maximum, ")"),
                          linkfun = function(mu) log((mu / maximum) / (1 - mu / maximum)),
                          linkinv = function(eta) maximum / (1 + exp(-eta)),
                          mu.eta = function(eta) maximum / (1 + exp(-eta)),
                          valideta = function(eta) TRUE), class = "link-glm")
    
    fam = binomial(link)
    fam$validmu = function(mu) all(mu > 0) && all(mu < maximum)
    fam$family = paste0("constrained.binomial (0, ", maximum, ")")
    
    # for mgcv
    fam$d2link <- function(mu) 1/(1 - (mu / maximum))^2 - 1/(mu / maximum)^2
    fam$d3link <- function(mu) 2/(1 - (mu / maximum))^3 + 2/(mu / maximum)^3
    fam$d4link <- function(mu) 6/(1 - (mu / maximum))^4 - 6/(mu / maximum)^4
    fam$dvar <- function(mu) rep.int(1, length(mu))
    fam$d3var <- fam$d2var <- function(mu) rep.int(0, length(mu))
    
    # new addition to initialization: mu cannot be greater than maximum
    new.line = substitute({mustart <- mustart * maximum}, list(maximum=maximum))
    fam$initialize <- c(fam$initialize, new.line)
    fam
}
