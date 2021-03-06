% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{simulate_fPi0_t_tests}
\alias{simulate_fPi0_t_tests}
\title{simulation of t-tests with functional proportion}
\usage{
simulate_fPi0_t_tests(shape = "Monotonic", m = 2500, n = 30, ...)
}
\arguments{
\item{shape}{True shape of pi0 relative to z0: either "Monotonic",
"Asymptotic", "Symmetric", "2 Step" or "3 Step"}

\item{m}{Number of hypotheses}

\item{n}{Sample size of each test}

\item{...}{Additional parameters passed on to simulateTTests}
}
\description{
Unlike simulateTTests, this simulation fixes the sample size of each t-test but
models the likelihood that a null hypothesis is true as a function of a random variable z0.
}
\details{
First, the random variable z0 is generated from a standard Normal distribution. Then the likelihood of a true null hypothesis pi0 is modelled as a function of z0.

z is calculated as rank(z0) / length(z0).
}
