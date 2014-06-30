\name{fqvalue}
\alias{fqvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate functional q-values}
\description{
Estimate the functional q-values for a set of p-values and a corresponding
vector X. A functional q-value takes the vector X into account by assuming that
the either the power of the test or the true effect size is positively
correlated with X.
}
\usage{
fqvalue(p, X, lambda = 0.5, method = "kernel", smooth.df = NULL, check.assumptions = TRUE, bandwidth = 0.01, pi0.bandwidth = 0.2, fixed.pi0 = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
A vector of p-values.
}
  \item{X}{
Some surrogate for power, where higher X indicates greater power. Must be the
same length as p.
}
  \item{lambda}{
The value of the tuning parameter to estimate pi_0 as a function of X.
}
  \item{method}{
The method used to estimate pi_0 as a function of X: either "kernel",
"spline", "logistic", or "LOESS".
}
  \item{smooth.df}{
If a spline is used to estimate the pi_0 function, use this degrees of freedom.
If not given, use cross-validation to determine the degrees of freedom.
}
  \item{check.assumptions}{
If TRUE (default), check that the p-values actually do tend to be smaller
with larger X, using a nonparametric spearman correlation test. Give a warning
if it does not appear to be the case (p > .05), in which case using functional
q-values might not be appropriate. 
}
  \item{bandwidth}{
The standard deviation of the Gaussian smoother for estimating the cumulative
density function.
}
  \item{pi0.bandwidth}{
If the pi_0 "kernel" estimation method is used, use a Gaussian smoother with
this standard deviation to smooth pi0.
}
  \item{fixed.pi0}{
If TRUE, assume that pi0 does not actually change as a function of X, only the
power of the test. If TRUE, estimate pi0 based on the largest values of X and
assume it holds for all tests. FALSE (default) is the conservative assumption.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a 
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
David Robisnon
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (p, X, lambda = 0.5, method = "kernel", smooth.df = NULL, 
    check.assumptions = TRUE, bandwidth = 0.01, pi0.bandwidth = 0.2, 
    fixed.pi0 = FALSE) 
{
    require(methods)
    require(data.table)
    fp = fpi0(p, X, lambda, method, smooth.df, pi0.bandwidth, 
        check.assumptions)
    pi0 = fp$pi0
    qX = fp$qX
    dt = data.table(pvalue = p, X = X, qX = qX, pi0 = pi0)
    if (fixed.pi0) {
        dt$pi0 = min(dt$pi0)
    }
    dt.ordered = dt[order(pvalue), ]
    densities = cumulative.densities(dt.ordered, bandwidth = bandwidth)
    dt$cumulative.density = densities[rank(dt$pvalue)]
    dt[, `:=`(fqvalue, pmin(pi0 * pvalue/cumulative.density, 
        1))]
    dt
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
