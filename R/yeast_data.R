#' Yeast eQTL data
#'
#' Data from an eQTL experiment on baker's yeast, with gene expressions in each of the 109 genotyped strains under two conditions, glucose and ethanol.
#'
#' @usage data(yeast)
#'
#' @format A RData object.
#'
#' @keywords datasets
#'
#' @references Smith, Erin N and Kruglyak, Leonid. (2008). Gene-environment interaction in yeast gene expression. PLoS Biology 6(4), e83.
#' (\href{https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0060083}{PLoS Biology})
#'
#' @source \href{https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0060083}{PLoS Biology}
#' 
#' @examples
#' library(fFDR)
#' data(yeast)
#' 
#' yeast_summary()
#' 
#' @export
yeast_summary <- function () {
  
  list(exp.e = exp.e(),
       exp.g = exp.g(),
       exp.pos = exp.pos(),
       marker = marker(),
       marker.pos = marker.pos()) %>%
    utils::str()
}

#' @rdname yeast_summary
#' @name exp.e
exp.e <- function () {
  utils::data(yeast)
  exp.e
}

#' @rdname yeast_summary
#' @name exp.g
exp.g <- function () {
  utils::data(yeast)
  exp.g
}

#' @rdname yeast_summary
#' @name exp.pos
exp.pos <- function () {
  utils::data(yeast)
  exp.pos
}

#' @rdname yeast_summary
#' @name marker
marker <- function () {
  utils::data(yeast)
  marker
}

#' @rdname yeast_summary
#' @name marker.pos
marker.pos <- function () {
  utils::data(yeast)
  marker.pos
}
