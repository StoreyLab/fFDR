\docType{data}
\name{simtests}
\alias{simtests}
\title{Simulated t-tests for the fqvalue package}
\format{a data frame, one row per hypothesis, with the p-value, sample size, and true mean of each sample.}
\usage{
simtests
}
\description{
Simulated p-values from 4000 one sample t-tests, 2000 of
them null and the rest with mean Norm(0, .6) and standard
deviation 1, with sample sizes drawn from 2 +
round(rlnorm(num, 2, 2)).
}
\details{
Seed used was 010814.
}
\author{
David Robinson, 01-08-14
}

