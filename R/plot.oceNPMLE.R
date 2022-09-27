#' Plot \code{oceNPMLE} object.
#'
#' @param x oceNPMLE object (see \code{\link{oceNPMLE}}).
#' @param xlab x label
#' @param ylab y label
#' @param ylim limits for the y axis, defaults to c(0,1)
#' @param col color vector, col[1] for group=0 and col[2] for group=1.
#' @param mark.time logical, should censored values be plotted?
#' @param ... Extra arguments (e.g., lwd=2) added to  \code{lines} functions.
#'
#' @export
#'
#' @return No return value, called for side effects.
#'
#' @examples
#' data(simScenario5)
#' dataFormt<-oceFormat(data=simScenario5, oceTime=c("T1","T2","T3"),
#'    oceStatus=c("I1","I2","I3"), group=c("Z"),
#'    oceNames = c("Death","Stroke/MI","Bleed"))
#' npmleOutput<- oceNPMLE(dataFormt)
#' plot(npmleOutput, xlab="Custom x label", mark.time=FALSE, lwd=2)
#' # can add lines from coxph output
#' coxOutput<- oceCoxph(dataFormt)
#' plot(coxOutput,linesonly=TRUE, col=c("orange","purple"),lwd=2)
#' legend("bottomleft",
#'    legend=c("grp=0, NPMLE","grp=1, NPMLE","grp=0, coxph","grp=1, coxph"),
#'    col=c("red","blue","orange","purple"),lty=c(1,1,1,1),lwd=2)
#'
plot.oceNPMLE<-function(x,
                        xlab="Ordering Score",
                        ylab="Proportion with a larger ordering score",
                        ylim=c(0,1),
                        col=c("red","blue"),mark.time=TRUE,...){

  KM0<-x$KM0
  KM1<-x$KM1
  TAU<-x$TAU
  oceNames<- x$oceNames
  k<-length(oceNames)
  # set lots of plotting parms to NULL, so we can give them



  # Do not plot the confidence intervals, they will not be correct, because we have ignored the correlation between endpoints

  if (length(col)==1){ col<-rep(col,2) }

  plot(c(0,k*TAU),c(0,1),type="n",las=1,ylab=ylab,xlab=xlab,ylim=ylim)
  lines(KM0,conf.int=F,col=col[1],mark.time=mark.time,...)
  lines(KM1,conf.int=F,col=col[2],mark.time=mark.time,...)

  k<-length(x$oceNames)
  TAU<- x$TAU
  for (j in 1:k){
    lines(c(j*TAU,j*TAU),c(0,1))
    text((j-0.5)*TAU,.9,oceNames[j])
  }

}
