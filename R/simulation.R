#' simulate one-sample t-tests
#' 
#' \code{simulateTTests} simulates a number of one-sample t-tests on
#' Normally distributed samples with varying mean mu0 and sample size.
#' 
#' @param m Number of hypotheses to simulate
#' @param pi0 Either a value for the fraction of true null hypotheses for which mu=0,
#' or a vector each of whose entries is the likelihood of a true null hypothesis
#' @param mu.sd Standard deviation used to generate the means of the Normal distributions under the
#' alternative hypothesis
#' @param mu.min The minumum value of abs(mu) for the means of the Normal distributions under the alternative hypothesis
#' @param n.lmean Mean of the sample sizes on a log scale
#' @param n.lsd Standard deviation of the sample sizes on a log scale
#' @param sample.sd Standard deviation for the Normal distributions
#' @param seed Optionally, random seed to set before simulating
#' @param ... Extra arguments, such as replication (may not be used)
#' 
#' @return A data frame, one row per hypothesis, with the columns
#' 
#' \item{p.value}{Two-sided t-test p-value}
#' \item{n}{Sample size of each sample}
#' \item{mu}{Mean of each sample}
#' \item{oracle}{TRUE if the alternative hypothesis holds, and FALSE for a true null hypothesis}
#' 
#' @details Each sample and test is constructed according to the following
#' generative process.
#' 
#' The mean of a Normal distribution, mu, is either 0 under the null hypothesis, or
#' generated as (\code{mu.min} + abs(Norm(0, \code{mu.sd})) under the alternative
#' hypothesis. So, the means are always either positive or 0.
#' \code{mu.min} can be used to ensure a minimal effect size between the null hypothesis and the
#' weakest alternative hypothesis (since by default an alternative hypothesis can
#' be arbitrarily close to being a true null hypothesis).
#' 
#' The sample size of each sample, n, is generated from a log-normal distribution, independent of the means of the Normal distributions. Specifically, it is generated as
#' 2 + round(LNorm(n.lmean, n.lsd)), such that the sample size is never less
#' than 2 (which would make a t-test impossible).
#' 
#' At this point, a random sample of size n is generated from the
#' Norm(mu, \code{sample.sd}) distribution. A one-sample two-sided t-test
#' is performed and the p-value obtained.
#' 
#' @export
simulate_t_tests = function(m = 4000, pi0 = .5, mu.sd = .3,
                            mu.min = 0, n.lmean = 2, n.lsd = 2.5,
                            sample.sd = 1, seed = NULL, ...) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    if (length(pi0) == 1) {
        neg <- round(pi0 * m)
        pos <- m - neg
        oracle <- c(rep(TRUE, pos), rep(FALSE, neg))
    } else if (length(pi0) == m) {
        oracle <- as.logical(rbinom(m, 1, 1 - pi0))
    } else {
        stop(paste("pi0 must either be a vector of length 1 or m"))
    }
    
    # generate sample sizes from a log-normal distribution
    # (though the smallest value is always 2)
    n = 2 + round(rlnorm(m, n.lmean, n.lsd))
    
    
    # setting up a gap between 0 and the smallest value of the alternative
    mu <- abs(rnorm(m, 0, mu.sd))
    mu <- (mu + mu.min) * oracle

    # Generate random samples and perform one-sample t-tests
    # for null (oracle == FALSE), just use uniform
    simdata <- lapply(1:sum(oracle), function(i) {
        rnorm(n[oracle][i], mu[oracle][i], sample.sd)
        })
    p.value <- runif(m)
    p.value[oracle] <- sapply(simdata, function(v) t.test(v)$p.value)

    return(data.frame(p.value = p.value, n = n, mu = mu, oracle = oracle))
}


#' simulation of t-tests with functional proportion
#' 
#' Unlike simulateTTests, this simulation fixes the sample size of each t-test but
#' models the likelihood that a null hypothesis is true as a function of a random variable z0.
#' 
#' @param shape True shape of pi0 relative to z0: either "Monotonic",
#' "Asymptotic", "Symmetric", "2 Step" or "3 Step"
#' @param m Number of hypotheses
#' @param n Sample size of each test
#' @param ... Additional parameters passed on to simulateTTests
#' 
#' @details First, the random variable z0 is generated from a standard Normal distribution. Then the likelihood of a true null hypothesis pi0 is modelled as a function of z0.
#' 
#' z is calculated as rank(z0) / length(z0).
#' 
#' @export
simulate_fPi0_t_tests = function(shape="Monotonic", m=2500, n=30, ...) {
    pi0.funcs <- list("Monotonic" = function(x) plogis(x),
                     "2 Step" = function(x) c(.4, .8)[cut(x, c(-Inf, 0, Inf))],
                     "3 Step" = function(x) c(.2, .5, .8)[cut(x, c(-Inf, -.5, .5, Inf))],
                     "Symmetric" = function(x) plogis(x ^ 2),
                     "Asymptotic" = function(x) plogis(3 * z0 + 3) * .75)

    z0 <- rnorm(m)
    z <- rank(z0) / m
    pi0 <- pi0.funcs[[shape]](z0)

    tt <- simulate_t_tests(m, pi0 = pi0, n.lmean = log(n - 2), n.lsd = 0, ...)
    cbind(tt, z = z, pi0 = pi0)
}

# 
# #' Factorial simulation on all combinations of simulateTTests and
# #' fqvalue parameters
# #' 
# #' Given a list of possible parameters to be passed to \code{simulateTTests} and
# #' parameters to be passed to \code{fqvalue}, perform a simulation that combines
# #' those parameters in all possible ways, performs simulated t-tests, then uses
# #' fqvalue to find functional q-values for each. Also use the qvalue
# #' package to find and report traditional q-values.
# #' 
# #' @param sim.pars A list of vectors of parameters that can be passed to
# #' simulateTTests
# #' @param fq.pars A list of vectors of parameters that can be passed to
# #' fqvalue, along with the p-value and n from the simulations
# #' @param replications Number of replications of each simulation
# #' @param fpi0 Whether to perform an fpi0 specific simulation, where
# #' pi0 is generated as a function of a latent variable z
# #' 
# #' @return A "Simulation" object, which is data.table, the first columns of
# #' which contain the simulation parameters and the fqvalue parameters used
# #' in each simulation, followed by
# #' 
# #' \item{p.value}{P-value from the simulation}
# #' \item{n}{Sample size of sample from the simulation, used as z in fqvalue}
# #' \item{qZ}{Quantiles of n (used as Z) in fqvalue package}
# #' \item{mu}{True mean of randomly generated sample}
# #' \item{oracle}{TRUE or FALSE on whether the sample was null (mu=0)}
# #' \item{qvalue}{Traditional q-value, calculated using the qvalue package}
# #' \item{qpi0}{Pi0 calculated with the qvalue package (constant across tests)}
# #' \item{fqvalue}{Functional q-value calculated with the fqvalue function}
# #' \item{fpi0}{Estimate of functional pi0 reported by fqvalue}
# #' 
# #' @examples
# #' 
# #' sim <- factorialSim(sim.pars = list(m = c(2500, 5000), pi0 = c(.5, .75)),
# #'                     fq.pars = list(pi0.method = c("kernel", "gam")))
# #' 
# #' @export
# factorialSim = function(sim.pars=list(), fq.pars=list(), replications=NULL,
#                         fpi0=FALSE) {
#     if (!is.null(replications)) {
#         sim.pars = c(sim.pars, list(replication = 1:replications))
#     }
#     par = as.data.table(expand.grid(sim.pars, stringsAsFactors = FALSE))
# 
#     if (fpi0) {
#         simfunc = simulate_fPi0_t_tests
#     } else {
#         simfunc = simulate_t_tests
#     }
#     # have to handle the case of no arguments separately
#     if (length(sim.pars) > 0) {
#         sim = par[, do.call(simfunc, mget(names(sim.pars))), by = names(sim.pars)]
#     } else {
#         sim = par[, simfunc(), by = names(par)]
#     }
# 
#     if (fpi0) {
#         sim = sim[, run.fq.params(p.value, n, mu, oracle, fq.pars, z, pi0), by = names(sim.pars)]
#     } else {
#         sim = sim[, run.fq.params(p.value, n, mu, oracle, fq.pars), by = names(sim.pars)]
#     }
#     
#     # return a simulation object
#     attr(sim, "parameters") <- c(names(sim.pars), names(fq.pars))
#     class(sim) <- "Simulation"
#     sim
# }
# 
# 
# run.fq.params = function(p.value, n, mu, oracle, fq.pars, z=NULL, pi0=NULL) {
#     # helper function for factorialSim
#     if (length(fq.pars) == 0) {
#         return(as.data.table(add.columns(p.value, n, mu, oracle, z = z, pi0 = pi0)))
#     }
#     innersim = as.data.table(expand.grid(fq.pars, stringsAsFactors = FALSE))
#     pn = names(fq.pars)
#     as.list(innersim[, add.columns(p.value, n, mu, oracle, mget(pn), z = z, pi0 = pi0), by = pn])
# }
# 
# 
# # #' Factorial simulation on all combinations of simulateTTests and
# # #' fqvalue parameters
# # #' 
# # #' Given a list of possible parameters to be passed to \code{simulateTTests} and
# # #' parameters to be passed to \code{fqvalue}, perform a simulation that combines
# # #' those parameters in all possible ways, performs simulated t-tests, then uses
# # #' fqvalue to find functional q-values for each. Also use the traditional qvalue
# # #' package to find and report q-values.
# # #' 
# # #' @param sim.pars A list of vectors of parameters that can be passed to
# # #' simulateTTests
# # #' @param fq.pars A list of vectors of parameters that can be passed to
# # #' fqvalue, along with the p-value and n from the simulations
# # #' @param replications Number of replications of each simulation
# # #' @param fpi0 Whether to perform an fpi0 specific simulation, where
# # #' pi0 is generated as a function of a latent variable z
# # #' 
# # #' @return A data.table, the first columns of which contain the simulation
# # #' parameters and the fqvalue parameters used in each simulation, followed by
# # #' 
# # #' \item{p.value}{P-value from the simulation}
# # #' \item{n}{Sample size of sample from the simulation, used as z in fqvalue}
# # #' \item{qZ}{Quantiles of n (used as Z) in fqvalue package}
# # #' \item{mu}{True mean of randomly generated sample}
# # #' \item{oracle}{TRUE or FALSE on whether the sample was null (mu=0)}
# # #' \item{qvalue}{Traditional q-value, calculated using the qvalue package}
# # #' \item{qpi0}{Pi0 calculated with the qvalue package (constant across tests)}
# # #' \item{fqvalue}{Functional q-value calculated with the fqvalue function}
# # #' \item{fpi0}{Estimate of functional pi0 reported by fqvalue}
# # #' 
# # #' @examples
# # #' 
# # #' sim = factorialSim(sim.pars=list(m=c(2500, 5000), pi0=c(.5, .75)),
# # #'                    fq.pars=list(pi0.method=c("kernel", "spline")))
# # #' 
# # #' @export
# # factorialSim2 = function(sim.pars=list(), fq.pars=list(), replications=NULL,
# #                         fpi0=FALSE) {
# #     if (!is.null(replications)) {
# #         sim.pars = c(sim.pars, list(replication=1:replications))
# #     }
# #     
# #     if (fpi0) {
# #         simfunc = simulatefPi0TTests
# #     } else {
# #         simfunc = simulateTTests
# #     }
# #     # have to handle the case of no arguments separately
# #     if (length(sim.pars) > 0) {
# #         sim = do.call(do_factorial, list(data.frame(), simulateTTests(), sim.pars))
# #     }
# #     else {
# #         sim = simulateTTests()
# #     }
# #     
# #     if (fpi0) {
# #         sim = sim[, run.fq.params(p.value, n, mu, oracle, fq.pars, z, pi0), by=names(sim.pars)]
# #     } else {
# #         sim = sim[, run.fq.params(p.value, n, mu, oracle, fq.pars), by=names(sim.pars)]
# #     }
# #     
# #     # return a simulation object
# #     ret = new("Simulation", table=sim, parameters=c(names(sim.pars), names(fq.pars)))
# #     ret
# # }
# 
# 
# add.columns = function(p.value, n, mu, oracle, pars=list(), z=NULL, pi0=NULL) {
#     # helper function for run.fq.params
#     # add columns to the data table for the qvalue, f-qvalue, and
#     # quantile of Z
#     # browser()
#     q = qvalue(p.value)
#     if (!is.null(z)) {
#         fp = do.call(estFPi0, c(list(p.value, z), pars))
#         res = fp$tableLambda
#         as.list(data.frame(p.value = p.value, n = n, mu = mu, oracle = oracle,
#                            z = z, pi0 = pi0, qpi0 = q$pi0, fpi0 = res$fpi0,
#                            lambda = res$lambda, phi.hat = res$phi.hat,
#                            chosen = res$chosen))
#     }
#     else {
#         fq = do.call(fqvalue, c(list(p.value, n), pars))
#         tab = fq$table
#         
#         # browser()
#         c(list(p.value=p.value, n=n, z=tab$z, mu=mu, oracle=oracle,
#              qvalue=q$qvalue, qpi0=q$pi0),
#           as.list(tab[, list(fpi0, lfdr, fq.value)]))
#     }
# }
# 
# 
# #' Summarize a factorial fqvalue simulation
# #' 
# #' @param object an fqvalueSimulation object
# #' @param alpha Desired confidence level
# #' @param ... Extra arguments (not used)
# #' 
# #' @return a data.table summarizing the simulation
# #' 
# #' @export
# summary.Simulation <- function(object, alpha = .05, ...) {
#     parnames = attr(object, "parameters")
#     #if (!("fqvalue" %in% colnames(object))) {
#     #    # pi0 only simulation
#     #    object[, list(qpi0=qpi0[1], fpi0.min=min(fpi0)), by=parnames]
#     #}
#     subcols = object[, c("oracle", "qpi0", "fpi0", "qvalue", "fqvalue",
#                          parnames), with = FALSE]
#     mtab = melt(subcols, id = c("oracle", "qpi0", "fpi0", parnames))
#     setnames(mtab, "variable", "method")
#     ret = mtab[, list(power = mean(value[oracle] < alpha),
#                       FDR = mean(!oracle[value < alpha]),
#                       AUC = wilcox.test(value[!oracle], value[oracle])$statistic / (sum(oracle) * sum(!oracle)),
#                       minpi0 = ifelse(method == "qvalue", min(qpi0), min(fpi0))),
#                by = c("method", parnames)]
# 
#     return(ret)
# }
