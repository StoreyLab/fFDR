library(fFDR)
simtests = simulateTTests(4000, .5, seed=010914)
save(simtests, file="../data/simtests.rda")
