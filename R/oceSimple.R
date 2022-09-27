#' Estimate win ratio or Mann-Whitney parameter using Simple Method
#'
#' Usually called from within \code{\link{oceTest}}.
#' Estimation done using simple method and output from \code{\link{oceCoxph}}.
#'
#' @param oceData output from \code{\link{oceFormat}}.
#' @param oceCoxOutput output from \code{\link{oceCoxph}}, if NULL recalculates
#'   using \code{oceData} and \code{\link{oceCoxph}}.
#'
#' @return A list with the following elements (where Yg=ordered composite
#'   endpoint score for group=g):
#' \describe{
#'  \item{int01}{estimate of P[Y0>Y1] (calculated from \code{\link{oceCoxph}})}
#'  \item{int10}{estimate of P[Y1>Y0] (calculated from \code{\link{oceCoxph}})}
#'  \item{WR}{win ratio, estimate of P[Y1>Y0]/P[Y0>Y1]}
#'  \item{MW}{desirability of outcome ranking,
#'    estimate of P[Y1>Y0]+(1/2)P[Y1=Y0]}
#'  }
#'
oceSimple<-function(oceData,oceCoxOutput=NULL){
  time<- oceData$timeMatrix
  status<- oceData$statusMatrix
  Z<- oceData$Z
  n<- nrow(time)
  k<- ncol(time)


  WIN<-0; LOSE<-0;
  index.i<- c(1:n)[Z==1]
  index.j<- c(1:n)[Z==0]

  for(i in index.i){
    for(j in index.j){
      atRisk<- TRUE
      WINij<-0
      LOSEij<-0
      for (h in 1:k){
        WINij<- WINij + (atRisk)*(time[i,h]>time[j,h])*(status[j,h]==1)
        LOSEij<- LOSEij + (atRisk)*(time[i,h]<time[j,h])*(status[i,h]==1)
        if (WINij>0 | LOSEij>0) atRisk<- FALSE
      }
      WIN  <- WIN  + WINij
      LOSE <- LOSE +  LOSEij
    }
  }

  if (is.null(oceCoxOutput)){
    oceCoxOutput<-oceCoxph(oceData)
  }
  int10<- oceCoxOutput$int10
  int01<- oceCoxOutput$int01
  MW  <- (WIN/(WIN+LOSE))*(int10+ int01) + (1/2)*(1-int10-int01)
  WR<-WIN/LOSE

  out<- list(MW=MW,WR=WR,int10=int10,int01=int01)
  class(out)<- "oceSimple"
  out
}
