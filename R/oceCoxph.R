#' Estimate win ratio or Mann-Whitney parameter using Cox Proportional Hazards
#'
#' Usually called from within \code{\link{oceTest}}, but useful for
#' getting \code{coxph} output details or customizing graphics. Estimation done using
#' coxph (partial likelihood methods).
#'
#' @param oceData output from \code{\link{oceFormat}}
#'
#' @return An \code{oceCoxph} object, which is a list with the following
#' elements (where Yg=ordered composite endpoint score for group=g):
#' \describe{
#'  \item{oceNames}{long names for each of the ordered endpoints}
#'  \item{TAU}{maximum of all the time-to-event variables (even censored ones)}
#'  \item{coxOutput}{output from coxph function}
#'  \item{int01}{estimate of P[Y0>Y1]}
#'  \item{int10}{estimate of P[Y1>Y0]}
#'  \item{WR}{win ratio, estimate of P[Y1>Y0]/P[Y0>Y1]}
#'  \item{MW}{desirability of outcome ranking,
#'     estimate of P[Y1>Y0]+(1/2)P[Y1=Y0]}
#'  }
#'
#' @importFrom survival coxph survfit
#' @importFrom stats coef
#' @seealso For an example using plotting see \code{\link{plot.oceCoxph}}.
#'  For Cox regression with other covariates, see
#' \code{vignette("Using oceCens",package="oceCens")}.
#'
#' @export
oceCoxph <- function(oceData) {
  IE <- oceData$status
  Z <- oceData$group
  START <- oceData$START
  STOP <- oceData$STOP
  regord <- coxph(Surv(START, STOP, IE) ~ Z)
  datfra <- data.frame("Z" = 0)
  regordfit <- survfit(regord, newdata = datfra)

  out <- oceSurv2WRMW(
    time0 = regordfit$time,
    surv0 = regordfit$surv,
    time1 = regordfit$time,
    surv1 = regordfit$surv ^ exp(regord$coef)
  )


  out <-
    list(
      oceNames = oceData$oceNames,
      TAU = oceData$TAU,
      coxphOutput = regord,
      int10 = out$int10,
      int01 = out$int01,
      WR = unname(exp(-coef(regord)) ),
      MW = out$MW
    )
  class(out) <- "oceCoxph"
  out

}
