#' @name simtests
#' @title Simulated t-tests for the fqvalue package
#' @description Simulated p-values from 4000 one sample t-tests,
#' 2000 of them null and the rest with mean Norm(0, .6) and standard
#' deviation 1, with sample sizes drawn from 2 + round(rlnorm(num, 2, 2)).
#' @docType data
#' @usage simtests
#' @details Seed used was 010814.
#' @format a data frame, one row per hypothesis, with the p-value, sample size, and true mean of each sample.
#' @author David Robinson, 01-08-14
NULL


#' simulate one-sample t-tests with a distribution of means and sample sizes
#' 
#' \code{simulateTTests} simulates a number of one-sample t-tests on
#' normally distributed samples of varying mean and sample size.
#' 
#' @param m Number of hypotheses to simulate
#' @param pi0 Either a value for the fraction of hypotheses where mu=0
#' (hypothesis is null) or a vector containing one pi0 value for each
#' hypothesis
#' @param mu.sd Standard deviation used to generate the true means in the
#' alternative hypothesis cases
#' @param mu.min The minumum value of abs(mu) in an alternative hypothesis test
#' (allowing a distinct gap between the null and the smallest alternative)
#' @param n.lmean Mean of the sample size distribution on a log scale
#' @param n.lsd Standard deviation of the sample size distribution on a log scale
#' @param sample.sd Standard deviation for generating each sample
#' @param seed Optionally, random seed to set before simulating
#' @param ... Extra arguments, such as replication (may not be used)
#' 
#' @return A data frame, one row per hypothesis, with the columns
#' 
#' \item{pvalue}{Two-sided t-test p-value}
#' \item{n}{Sample size of each sample}
#' \item{mu}{True mean of each sample}
#' \item{oracle}{TRUE if the alternative hypothesis holds, FALSE for null}
#' 
#' @details Each sample and test is constructed according to the following
#' generative process.
#' 
#' The true mean of each sample, mu, is either 0, for null hypotheses, or
#' generated as (\code{mu.min} + abs(Norm(0, \code{mu.sd})), for alternative
#' hypotheses. This means the true mean is always either positive or 0.
#' \code{mu.min} can be used to ensure a gap between the null hypotheses and the
#' weakest alternative hypotheses- by default an alternative hypothesis can
#' be arbitrarily close to being a null.
#' 
#' The sample size of each sample, n, is generated independently of the mean, from
#' a log-normal distribution. Specifically it is generated as
#' 2 + round(LNorm(n.lmean, n.lsd)), such that the sample size is never less
#' than 2 (which would make a t-test impossible).
#' 
#' At this point, a random sample of size n is generated from the
#' Norm(mu, \code{sample.sd}) distribution. Increasing the sample standard
#' deviation decreases the power of all tests. A one-sample two-sided t-test
#' is performed and the p-value reported.
#' 
#' @export
simulateTTests = function(m=4000, pi0=.5, mu.sd=.3, mu.min=0,
                          n.lmean=2, n.lsd=2.5, sample.sd=1, seed=NULL,
                          ...) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    if (length(pi0) == 1) {
        neg = round(pi0 * m)
        pos = m - neg
        oracle = c(rep(TRUE, pos), rep(FALSE, neg))
    } else if (length(pi0) == m) {
        oracle = as.logical(rbinom(m, 1, 1 - pi0))
    } else {
        stop(paste("pi0 must either be a vector of length 1 or m"))
    }
    
    # generate sample sizes from a log-normal distribution
    # (though the smallest value is always 2)
    n = 2 + round(rlnorm(m, n.lmean, n.lsd))
    
    # generate mu as a mixture of 0 and a normal distribution, optionally
    # setting up a gap between 0 and the smallest value of the alternative
    mu = abs(rnorm(m, 0, mu.sd))
    mu = (mu + mu.min) * oracle

    # Generate random samples and perform one-sample t-tests
    # for null (oracle == FALSE), just use uniform
    simdata = lapply(1:sum(oracle), function(i) {
        rnorm(n[oracle][i], mu[oracle][i], sample.sd)
        })
    pvalue = runif(m)
    pvalue[oracle] = sapply(simdata, function(v) t.test(v)$p.value)

    return(data.frame(pvalue=pvalue, n=n, mu=mu, oracle=oracle))
}


#' simulation of t-tests with pi0 varying as a function of some variable z
#' 
#' Unlike simulateTTests, this fixes the sample size of each t-test, instead
#' varying the probability that the hypothesis is true.
#' 
#' @param shape True shape of pi0 relative to z0: either "Monotonic",
#' "Asymptotic", "Symmetric", "2 Step" or "3 Step"
#' @param m Number of hypotheses
#' @param n sample size of each test
#' @param ... Additional parameters passed on to simulateTTests
#' 
#' @details First a latent variable z0 is generated from a standard normal. pi0
#' is computed according to a function of z0.
#' 
#' z is calculated as rank(z0) / length(z0) (as is done by estFPi0).
#' 
#' @export
simulatefPi0TTests = function(shape="Monotonic", m=2500, n=30, ...) {
    pi0.funcs = list("Monotonic" = function(x) plogis(x),
                     "2 Step" = function(x) c(.4, .8)[cut(x, c(-Inf, 0, Inf))],
                     "3 Step" = function(x) c(.2, .5, .8)[cut(x, c(-Inf, -.5, .5, Inf))],
                     "Symmetric" = function(x) plogis(x^2),
                     "Asymptotic"= function(x) plogis(3 * z0+3) * .75)

    z0 = rnorm(m)
    z = rank(z0) / m
    pi0 = pi0.funcs[[shape]](z0)

    tt = simulateTTests(m, pi0=pi0, n.lmean=log(n - 2), n.lsd=0, ...)
    cbind(tt, z=z, pi0=pi0)
}


#' Factorial simulation on all combinations of simulateTTests and
#' fqvalue parameters
#' 
#' Given a list of possible parameters to be passed to \code{simulateTTests} and
#' parameters to be passed to \code{fqvalue}, perform a simulation that combines
#' those parameters in all possible ways, performs simulated t-tests, then uses
#' fqvalue to find functional q-values for each. Also use the traditional qvalue
#' package to find and report q-values.
#' 
#' @param sim.pars A list of vectors of parameters that can be passed to
#' simulateTTests
#' @param fq.pars A list of vectors of parameters that can be passed to
#' fqvalue, along with the pvalue and n from the simulations
#' @param replications Number of replications of each simulation
#' @param fpi0 Whether to perform an fpi0 specific simulation, where
#' pi0 is generated as a function of a latent variable z
#' 
#' @return A data.table, the first columns of which contain the simulation
#' parameters and the fqvalue parameters used in each simulation, followed by
#' 
#' \item{pvalue}{P-value from the simulation}
#' \item{n}{Sample size of sample from the simulation, used as z in fqvalue}
#' \item{qZ}{Quantiles of n (used as Z) in fqvalue package}
#' \item{mu}{True mean of randomly generated sample}
#' \item{oracle}{TRUE or FALSE on whether the sample was null (mu=0)}
#' \item{qvalue}{Traditional q-value, calculated using the qvalue package}
#' \item{qpi0}{Pi0 calculated with the qvalue package (constant across tests)}
#' \item{fqvalue}{Functional q-value calculated with the fqvalue function}
#' \item{fpi0}{Estimate of functional pi0 reported by fqvalue}
#' 
#' @examples
#' 
#' sim = factorialSim(sim.pars=list(m=c(2500, 5000), pi0=c(.5, .75)),
#'                    fq.pars=list(pi0.method=c("kernel", "spline")))
#' 
#' @export
factorialSim = function(sim.pars=list(), fq.pars=list(), replications=NULL,
                        fpi0=FALSE) {
    if (!is.null(replications)) {
        sim.pars = c(sim.pars, list(replication=1:replications))
    }
    par = as.data.table(expand.grid(sim.pars, stringsAsFactors=FALSE))

    if (fpi0) {
        simfunc = simulatefPi0TTests
    } else {
        simfunc = simulateTTests
    }
    # have to handle the case of no arguments separately
    if (length(sim.pars) > 0) {
        sim = par[, do.call(simfunc, mget(names(sim.pars))), by=names(sim.pars)]
    }
    else {
        sim = par[, simfunc(), by=names(par)]
    }

    if (fpi0) {
        sim = sim[, run.fq.params(pvalue, n, mu, oracle, fq.pars, z, pi0), by=names(sim.pars)]
    } else {
        sim = sim[, run.fq.params(pvalue, n, mu, oracle, fq.pars), by=names(sim.pars)]
    }
    
    # return a simulation object
    ret = new("Simulation", table=sim, parameters=c(names(sim.pars), names(fq.pars)))
    ret
}


run.fq.params = function(pvalue, n, mu, oracle, fq.pars, z=NULL, pi0=NULL) {
    # helper function for factorialSim
    if (length(fq.pars) == 0) {
        return(as.data.table(add.columns(pvalue, n, mu, oracle, z=z, pi0=pi0)))
    }
    innersim = as.data.table(expand.grid(fq.pars, stringsAsFactors=FALSE))
    pn = names(fq.pars)
    as.list(innersim[, add.columns(pvalue, n, mu, oracle, mget(pn), z=z, pi0=pi0), by=pn])
}


add.columns = function(pvalue, n, mu, oracle, pars=list(), z=NULL, pi0=NULL) {
    # helper function for run.fq.params
    # add columns to the data table for the qvalue, f-qvalue, and
    # quantile of Z
    q = qvalue(pvalue)
    if (!is.null(z)) {
        fp = do.call(estFPi0, c(list(pvalue, z), pars))
        res = fp@tableLambda
        as.list(data.frame(pvalue=pvalue, n=n, mu=mu, oracle=oracle, z=z, pi0=pi0, qpi0=q$pi0, fpi0=res$fpi0, lambda=res$lambda, phi.hat=res$phi.hat, chosen=res$chosen))
    }
    else {
        fq = do.call(fqvalue, c(list(pvalue, n), pars))
        tab = fq@table
        
        c(list(pvalue=pvalue, n=n, z=tab$z, mu=mu, oracle=oracle,
             qvalue=q$qvalue, qpi0=q$pi0),
          as.list(tab[, list(fpi0, lfdr, fqvalue)]))
    }
}


#' Summarize a factorial fqvalue simulation
#' 
#' @param object an fqvalueSimulation object
#' @param Desired confidence level
#' 
#' @return a data.table summarizing the simulation
#' 
#' @export
setMethod("summary", "Simulation", function(object, alpha=.05, ...) {
    parnames = object@parameters
    #if (!("fqvalue" %in% colnames(object))) {
    #    # pi0 only simulation
    #    object[, list(qpi0=qpi0[1], fpi0.min=min(fpi0)), by=parnames]
    #}
    subcols = object@table[, c("oracle", "qpi0", "fpi0", "qvalue", "fqvalue",
                               parnames), with=FALSE]
    mtab = melt(subcols, id=c("oracle", "qpi0", "fpi0", parnames))
    setnames(mtab, "variable", "method")
    ret = mtab[, list(power=mean(value[oracle] < alpha),
                      FDR=mean(!oracle[value < alpha]),
                      minpi0=ifelse(method == "qvalue", min(qpi0), min(fpi0))),
               by=c("method", parnames)]

    return(ret)
})
