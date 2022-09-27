#' Estimate win ratio or Mann-Whitney parameter using NPMLE
#'
#' Estimation done using NPMLE (nonparametric maximum likelihood estimators of
#' survival).
#'
#' @param oceData output from \code{\link{oceFormat}}
#'
#' @return An object of class 'oceNPMLE', which is a list with the following
#' elements (where Yg=ordered composite endpoint score for group=g):
#' \describe{
#'  \item{oceNames}{long names for each of the ordered endpoints}
#'  \item{TAU}{maximum of all the time-to-event variables (even censored ones)}
#'  \item{KM0}{survfit output for group=0 subset}
#'  \item{KM1}{survfit output for group=1 subset}
#'  \item{WR}{win ratio, estimate of P[Y1>Y0]/P[Y0>Y1]}
#'  \item{MW}{desirability of outcome ranking, estimate of P[Y1>Y0]+(1/2)P[Y1=Y0]}
#'  }
#'
#' @export
#'
#' @seealso See \code{\link{plot.oceNPMLE}} for an example with plotting.
oceNPMLE<-function(oceData){

  IE<- oceData$status
  ZZ<-oceData$group
  START<- oceData$START
  STOP<- oceData$STOP
  # Calculate NPMLE for Group Z=0
  KM0<-survfit( Surv(START[ZZ==0],STOP[ZZ==0],IE[ZZ==0])~1)
  # Calculate NPMLE for Group Z=1
  KM1<-survfit( Surv(START[ZZ==1],STOP[ZZ==1],IE[ZZ==1])~1)

  summary0<- summary(KM0)
  summary1<- summary(KM1)

  out<- oceSurv2WRMW(time0=summary0$time,
                     surv0=summary0$surv,
                     time1=summary1$time,
                     surv1=summary1$surv)



  out<-list(oceNames=oceData$oceNames,TAU=oceData$TAU,KM0=KM0,KM1=KM1,WR=out$WR,MW=out$MW)
  class(out)<-"oceNPMLE"
  out
}
