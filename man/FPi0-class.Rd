\docType{class}
\name{FPi0-class}
\alias{FPi0-class}
\title{An estimation of functional pi0 on some variable z}
\description{
This is an estimation of pi0 depending on a variable z for
a particular combination of p-values and z. Aside from the
functional pi0 estimates, the class contains information on
the selection of the tuning parameter lambda, which
controls a bias/variance tradeoff.
}
\section{Slots}{
  \describe{ \item{\code{table}: a \code{data.table}
  containing the p-values, z, z0, and estimate fpi0}
  \item{\code{tableLambda}: a \code{data.table} containing
  the fPi0 estimates for each choice of tuning parameter
  lambda} \item{\code{MISE}: a \code{data.table} containing
  estimates of bias, variance and mean integrated squared
  error (MISE) for each choice of lambda}
  \item{\code{lambda}: choice of tuning parameter lambda}
  \item{\code{method}: character string indicating the
  method used to fit the shape} \item{\code{fits}: list of
  objects produced by fitting at each lambda (i.e. glmfit
  or gamfit objects)} }
}

