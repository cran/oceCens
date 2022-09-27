#' Take coxph object and translate results to win ratios.
#'
#' Let \code{cout} a coxph object, then Using normal approximations and the output from the
#' \code{cout$coefficients} and \code{cout$var}. If the cluster argument is used in the coxph
#' call, then \code{cout$var} is the robust variance (see \code{\link[survival]{coxph}}.
#'
#' @param coutput a coxph object created by \code{\link[survival]{coxph}}.
#' @param conf.level confidence level.
#'
#' @details
#' The function takes a beta coefficient and returns the win ratio version: exp(-beta).
#' Confidence intervals are calculated by
#' exp(-beta -/+ qnorm(1-(1-conf.level)/2)*sqrt(coutput$var)).
#' P-values are two-sided.
#'
#' @return A vector or matrix with 4 elements (or columns) giving the win ratio,
#' the lower and upper confidence limits, and the two-sided p-value.
#'
#' @export
#'
#' @importFrom stats pnorm qnorm
#'
#' @references
#' Follmann, D., Fay, M. P., Hamasaki, T., and Evans, S. (2020). Analysis of
#' ordered composite endpoints. Statistics in Medicine, 39(5), 602-616.
#'
#' @examples
#' data(simScenario5)
#' xform<-oceFormat(data=simScenario5,oceTime=c("T1","T2","T3"),
#'    oceStatus=c("I1","I2","I3"),
#'    group="Z",outputDataFrame=TRUE)
#' # perform cox regression using time varying treatment efects, IZ1,IZ2, IZ3
#' # associated with the 3 prioritized endpoints
#' cout<- coxph(Surv(START, STOP, status) ~ IZ1+IZ2+IZ3, data=xform$data)
#' coxph2WR(cout)
coxph2WR<-function(coutput,conf.level=0.95){
  # outputs parameter estimates as exp( - coef)
  # to give 1/HR
  WR<- exp( - coutput$coefficients)
  sd<- sqrt(diag(coutput$var))
  wald<- - coutput$coefficients/sd
  # calculate two-sided p-value based on the wald statistic
  p.value<-  2*(1-pnorm(abs(wald)))
  a<- 1-conf.level
  upperCL<- exp( - coutput$coefficients + qnorm(1-a/2)*sd)
  lowerCL<- exp( - coutput$coefficients - qnorm(1-a/2)*sd)
  output<- cbind(WR,lowerCL,upperCL,p.value)
  attr(output,"conf.level")<- conf.level
  output
}
