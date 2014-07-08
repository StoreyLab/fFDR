##################################
######## AllClasses.R
########
########    all classes in fFDR
##################################
####################################################################


#' An estimation of functional pi0 on some variable z
#' 
#' This is an estimation of pi0 depending on a variable z for a particular
#' combination of p-values and z. Aside from the functional pi0 estimates,
#' the class contains information on the selection of the tuning parameter
#' lambda, which controls a bias/variance tradeoff.
#' 
#'@section Slots:
#'  \describe{
#'    \item{\code{table}:}{a \link{data.table} containing the p-values, z, z0, and estimate fpi0}
#'    \item{\code{tableLambda}:}{\link{data.table} containing the fPi0 estimates for each
#' choice of tuning parameter lambda}
#'    \item{\code{MISE}:}{a \link{data.table} containing estimates of bias, variance and mean
#' integrated squared error (MISE) for each choice of lambda}
#'    \item{\code{lambda}:}{choice of tuning parameter lambda}
#'    \item{\code{method}:}{character string indicating the method used to fit the shape}
#' }
#' 
#' 
setClass("fPi0",
         representation(
             table = "data.table",
             tableLambda = "data.table",
             MISE = "data.table",
             lambda = "numeric",
             method = "character"
         )
)

setClassUnion("fPi0OrNULL", c("fPi0", "NULL"))

#' Functional q-values calculated for a set of p-values and z
#' 
#' This controls false discovery rate taking into consideration the
#' dependence of p-values on a factor z.
#' 
#'@section Slots:
#'  \describe{
#'    \item{\code{table}:}{a \link{data.table} containing the p-values, z, z0, estimated fpi0, estimated density, local FDR, and the computed f-qvalues.}
#'    \item{\code{fPi0}:}{An \linkS4class{fPi0} object, or NULL if functional pi0 estimation was not used}
#'    \item{\code{density}:}{A \link{data.table} containing the estimated density of each point in the p-value/z space, used to calculate the local FDR}
#' }
#' 
#' 
setClass("fqvalue",
         representation(
             table = "data.table",
             fPi0 = "fPi0OrNULL",
             density = "data.table"
         )
)


#' Factorial simulation of the fqvalue package
#' 
#' The results of the fqvalue package, including power and FDR, on simulated
#' t-tests, where the simulation and fqvalue parameters are set up in
#' factorial combinations
#' 
#'@section Slots:
#'  \describe{
#'    \item{\code{table}:}{a \link{data.table} containing the p-values, z, z0, estimated fpi0, estimated density, local FDR, and the computed f-qvalues.}
#'    \item{\code{fPi0}:}{An \linkS4class{fPi0} object, or NULL if functional pi0 estimation was not used}
#'    \item{\code{density}:}{A \link{data.table} containing the estimated density of each point in the p-value/z space, used to calculate the local FDR}
#' }
#' 
#' 
setClass("Simulation",
         representation(
             table = "data.table",
             parameters = "character"
         )
)
