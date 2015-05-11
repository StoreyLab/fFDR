library(testthat)
library(ggplot2)
library(data.table)

set.seed(2014)
sim.ttests = simulate_t_tests(m = 1000)

context("fqvalue")

test_consistent_fqvalue = function(fq, p, z0) {
    # function for determining whether a the result of an fqvalue call
    # is consistent with the pvalues and z0 that were given to it
    expect_is(fq, "fqvalue")
    expect_is(fq$table, "tbl_df")
    expect_equal(fq$table$p.value, p)
    expect_true(!is.null(fq$table$p.value))
    expect_true(!is.null(fq$table$fq.value))
    expect_equal(rank(fq$table$z), rank(z0))
    expect_that(as.numeric(fq), equals(fq$table$fq.value))
    expect_that(fq$fPi0, is_a("fPi0"))
}

test_that("fqvalue returns an object with the right structure", {
    fq = fqvalue(sim.ttests$p.value, sim.ttests$n)
    
    test_consistent_fqvalue(fq, sim.ttests$p.value, sim.ttests$n)
    
    # test that summary can be performed
    s = summary(fq, sim.ttests$oracle)
    expect_less_than(s$q.value.power, s$fq.value.power)
    # should be less than .05, give it some breathing room
    expect_less_than(s$fq.value.fdr, .2)
    # expect the FDR is not significantly higher than you'd
    # expect by chance
    expect_less_than(.005, s$FDR.binom.pval)

    # test that plots can be built
    print(plot(fq))
    print(compare_qvalue(fq))
    print(plot_MISE(fq$fPi0))
})

test_that("Monotone smoothing function works", {
    fqm = fqvalue(sim.ttests$p.value, sim.ttests$n, monotone.window = .1)
    test_consistent_fqvalue(fqm, sim.ttests$p.value, sim.ttests$n)
    tab <- fqm$table
    
    for (i in 1:nrow(tab)) {
        # use .099 to avoid floating point error getting in the way
        cond = (tab$p.value < tab$p.value[i] &
                     abs(tab$z - tab$z[i]) < .099)
        if (!all(tab$fx[cond] >= tab$fx[i])) {
            browser()
        }
        expect_true(all(tab$fx[cond] >= tab$fx[i]))
    }
})


test_that("fqvalue works on null hypotheses", {
    set.seed(2014)
    nullpvals = runif(1000)
    nullz = runif(1000)
    fqn = fqvalue(nullpvals, nullz)
    test_consistent_fqvalue(fqn, nullpvals, nullz)

    # should be no false discoveries, allow 2 anyway
    expect_less_than(sum(as.numeric(fqn) < .1), 3)

    # check it can be plotted
    print(plot(fqn))
    print(compare_qvalue(fqn))
    print(plot_MISE(fqn$fPi0))
})


context("fPi0")

set.seed(2015)
simfPi0 = simulate_fPi0_t_tests()

test_that("estimate_fpi0 returns the right kind of object for all methods", {
    # test all methods
    methods <- c("kernel", "glm", "gam", "bin")
    fpi0s <- lapply(methods, function(m) {
        estimate_fpi0(simfPi0$p.value, simfPi0$z, method = m)
    })

    # test for consistency
    for (fp in fpi0s) {
        expect_is(fp, "fPi0")
        expect_equal(length(as.numeric(fp)), nrow(simfPi0))
        expect_equal(fp$table$z, simfPi0$z)  
        expect_equal(fp$table$z0, simfPi0$z)
        
        expect_true(all(fp$table$fpi0 <= 1 & fp$table$fpi0 >= 0))
    }

    # test that the fpi0 values are similar
    pi0.matrix = sapply(fpi0s, function(fp) fp$table$fpi0)
    expect_less_than(.85, min(cor(pi0.matrix, method = "spearman")))
    
    # FPi0 plots are built without errors
    for (fp in fpi0s) {
        print(plot(fp))
        print(plot_MISE(fp))
    }
})

test_that("Incorrect usage of estimate_fpi0 throws errors", {
    expect_that(estimate_fpi0(simfPi0$p.value, simfPi0$z, method = "nomethod"),
                throws_error("should be one of"))
    expect_that(estimate_fpi0(simfPi0$p.value + 1, simfPi0$z),
                throws_error("valid range"))
})

