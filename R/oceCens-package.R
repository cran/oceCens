#' oceCens: A package for ordered composite endpoints with censoring.
#'
#'  An ordered composite endpoint combines several time-to-event endpoints
#'  into one score. The package compares two groups with two parameters, the
#'  win ratio, P[Y1>Y0]/P[Y0>Y1], and the Mann-Whitney parameter,
#'  P[Y1>Y0]+(1/2)P[Y1=Y0], where Y1 and Y0 are the oce scores in the two groups.
#'  The main function is \code{\link{oceTest}}, which calls many of the other
#'  functions and has several different methods for estimation.
#'  Statistical details are in Follmann, et al 2020.
#'
#'
#' @references Follmann, D., Fay, M. P., Hamasaki, T., and Evans, S. (2020). Analysis of
#' ordered composite endpoints. Statistics in Medicine, 39(5), 602-616.
#'
#' @docType package
#' @name oceCens
#' @aliases oceCens-package
NULL
