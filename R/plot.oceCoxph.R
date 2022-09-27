#' Plot \code{oceNPMLE} object.
#'
#' @param x oceCoxph object (see \code{\link{oceCoxph}}).
#' @param linesonly logical, add lines to an existing plot?
#' @inheritParams plot.oceNPMLE
#' @param ... Extra arguments (e.g., lwd=3) added to both \code{lines} functions.
#'
#' @return The function invisibly (see \code{\link[base]{invisible}}) returns
#'   a list with 4 elements: (time0, surv0, time1, and surv1)
#' @export
#'
#' @seealso Example in \code{\link{plot.oceNPMLE}} shows adding lines from
#'   the coxph output to an existing plot.
#'
#' @importFrom graphics lines text
#' @examples
#' # need to first run oceFormat and oceCoxph
#' data(simScenario5)
#' dataFormt<-oceFormat(data=simScenario5, oceTime=c("T1","T2","T3"),
#'    oceStatus=c("I1","I2","I3"), group=c("Z"),
#'    oceNames = c("Death","Stroke/MI","Bleed"))
#' coxOutput<- oceCoxph(dataFormt)
#' plot(coxOutput, xlab="Custom x label")
#'
plot.oceCoxph<-function(x, linesonly=FALSE,
                        xlab="Ordering Score",
                        ylab="Proportion with a larger ordering score",
                        col=c("red","blue"),...){

  datfra <- data.frame("Z" = 0)
  cout <- survfit(x$coxphOutput, newdata = datfra)
  time0 <- time1<- cout$time
  surv0 <- cout$surv
  surv1 <- cout$surv^exp(x$coxphOutput$coef)
  TAU<-x$TAU
  oceNames<- x$oceNames
  k<- length(oceNames)
  # set lots of plotting parms to NULL, so we can give them



  # Do not plot the confidence intervals, they will not be correct, because we have ignored the correlation between endpoints

  if (length(col)==1){ col<-rep(col,2) }

  if (!linesonly) plot(c(0,k*TAU),c(0,1),type="n",ylim=c(0,1),las=1,ylab=ylab,xlab=xlab)

  # plot lines in steps (like Kaplan-Meier or cdf) instead of connecting
  # points
  steplines<-function(time,surv,COL=col[1]){
    m<- length(surv)
    if (length(time)!=m) stop("time length does not match surv length")
    surv<- c(1,1,rep(surv[-m], each=2),surv[m])
    time<- c(0,rep(time[-m], each=2),time[m],Inf)
    lines(time,surv,col=COL,...)
  }

  steplines(time0,surv0,COL=col[1])
  steplines(time1,surv1,COL=col[2])

  # add oceNames to each area
  if (!linesonly){
    for (j in 1:k){
      lines(c(j*TAU,j*TAU),c(0,1))
      text((j-0.5)*TAU,.9,oceNames[j])
    }
  }
  invisible(list(time0=time0,surv0=surv0,time1=time1,surv1=surv1))
}
