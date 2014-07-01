library(ggplot2)
library(data.table)

set.seed(2014)
sim.ttests = simulateTTests()

context("fqvalue")

test_consistent_fqvalue = function(fq, p, z0) {
    # function for determining whether a the result of an fqvalue call
    # is consistent with the pvalues and z0 that were given to it
    expect_that(fq, is_a("fqvalue"))
    expect_that(fq@table$pvalue, equals(p))
    expect_that(colnames(fq@table), matches("pi0", all=FALSE))
    expect_that(colnames(fq@table), matches("fqvalue", all=FALSE))
    expect_that(rank(fq@table$z), equals(rank(z0)))
    expect_that(as.numeric(fq), equals(fq@table$fqvalue))
    expect_that(fq@fPi0, is_a("fPi0"))
}

test_that("fqvalue returns an object with the right structure", {
    fq = fqvalue(sim.ttests$pvalue, sim.ttests$n)
    
    test_consistent_fqvalue(fq, sim.ttests$pvalue, sim.ttests$n)
    
    # test that plots can be built
    p1 = plot(fq)
    p2 = compareQvalue(fq)
    p3 = plotMISE(fq)
})


context("fPi0")

set.seed(2014)
simfPi0 = simulatefPi0TTests()

test_that("estFPi0 returns the right kind of object for all methods", {
    # test all methods
    methods = c("kernel", "glm", "gam", "bin")
    fpi0s = lapply(methods, function(m) {
                    estFPi0(simfPi0$pvalue, simfPi0$z, method=m)
    })
    
    # test for consistency
    for (fp in fpi0s) {
        expect_that(fp, is_a("fPi0"))
        expect_that(length(as.numeric(fp)), equals(nrow(simfPi0)))
        expect_that(fp@table$z, equals(simfPi0$z))        
        expect_that(fp@table$z0, equals(simfPi0$z))
        
        expect_that(all(fp@table$fpi0 <= 1 & fp@table$fpi0 >= 0), is_true())
    }

    # test that the fpi0 values are similar
    pi0.matrix = sapply(fpi0s, function(fp) fp@table$fpi0)
    expect_that(min(cor(pi0.matrix, method="spearman")) > .85, is_true())
    
    # FPi0 plots are built without errors
    for (fp in fpi0s) {
        p1 = plot(fp)
        p2 = plotMISE(fp)
    }
})

test_that("Incorrect usage of estFPi0 throws errors", {
    expect_that(estFPi0(simfPi0$pvalue, simfPi0$z, method="nomethod"),
                throws_error("should be one of"))
    expect_that(estFPi0(simfPi0$pvalue + 1, simfPi0$z),
                throws_error("valid range"))
})

