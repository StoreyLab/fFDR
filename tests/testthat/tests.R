library(ggplot2)
library(data.table)

data(simtests)

# test_that("fqvalue returns a data table with the right structure", {
#     p = simtests$pvalue
#     x = simtests$sample.size
#     fq = fqvalue(p, x)
# 
#     expect_that(fq, is_a("data.table"))
#     expect_that(fq$pvalue, equals(p))
#     expect_that(fq$X, equals(x))
#     expect_that(rank(fq$X), equals(rank(fq$qX)))
# 
#     expect_that(colnames(fq), matches("pi0", all=FALSE))
#     expect_that(colnames(fq), matches("fqvalue", all=FALSE))
# 
#     expect_that(fqvalue(p, x, lambda=.5, method="nomethod"), throws_error())
# })

set.seed(2014)
simdata = simulatefPi0TTests()

test_that("estFPi0 returns the right kind of object for all methods", {
    # test all methods
    methods = c("kernel", "glm", "gam", "bin")
    fpi0s = lapply(methods, function(m)
                                estFPi0(simdata$pvalue, simdata$z, method=m))
    
    # test for consistency
    for (fp in fpi0s) {
        expect_that(fp, is_a("FPi0"))
        expect_that(length(as.numeric(fp)), equals(length(p)))
        expect_that(fp@table$z, equals(z))        
        expect_that(fp@table$z0, equals(z))
        
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


test_that("simulation of one-sample t-tests works", {
    
})

#fq = fqvalue(simtests$pvalue, simtests$sample.size, fixed.pi0=TRUE)
#print(plot.pi0(fq))
#fq2 = fqvalue(simtests2$pvalue, simtests2$sample.size)
#print(plot.pvalue.qX(fq2))
#print(g)
