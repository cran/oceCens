#' Converts Survival Fits to Win Ratio and Mann-Whitney Estimates
#'
#' @param time0 vector of times for group=0 subset
#' @param surv0 vector of survival distribution values for group=0
#' @param time1 vector of times for group=1 subset
#' @param surv1 vector of survival distribution values for group=1
#'
#' @return A list with the following elements (where Yg=ordered composite endpoint score for group=g):
#' \describe{
#'  \item{int01}{estimate of P[Y0>Y1]}
#'  \item{int10}{estimate of P[Y1>Y0]}
#'  \item{WR}{win ratio, int10/int01}
#'  \item{MW}{estimate of P[Y1>Y0]+(1/2)P[Y1=Y0]}
#'  }
#'
#' WR=WR,MW=MW,int10=int10,int01=int01
#' @export
#'
#' @details
#' This is an interval function called by \code{\link{oceNPMLE}} or
#' \code{\link{oceCoxph}}.
#'
#' @keywords internal
oceSurv2WRMW<-function(time0,surv0,time1,surv1){
  II0   <-seq(1:length(time0))

  #  Survival function is  evaluated at a specific time [0,oo),
  S0 <- function(s){
    ifelse(s<min(time0),1, surv0[max(II0[time0<=s])])
  }

  # Jump is evaluated at an index corresponding to the ith event time in group 0.  Returns NA if arguemnt is not in (1,2,...length(surv0))
  dF0 <- function(i){
    ifelse(i==1,1-surv0[i],surv0[i-1] -surv0[i])
  }


  #  Get the Survival functions and jump functions for group 1
  II1   <-seq(1:length(time1))

  # Survival function is evaluated at a specific time (0,oo)
  S1 <- function(s){
    ifelse(s<min(time1),1, surv1[max(II1[time1<=s])])
  }

  # Jump evaluated at an index corresponing to the ith event time in group 1.  Returns NA if arguemnt is not in (1,2,...length(surv1))
  dF1 <- function(i){
    ifelse(i==1,1-surv1[i],surv1[i-1] -surv1[i])
  }

  # do numeric integrations
  J0<-length(time0)
  int10<-0
  for(i in 1:J0){
    int10 <- int10 + S1(time0[i])*dF0(i)
  }
  J1<-length(time1)
  int01<-0
  for(i in 1:J1){
    int01 <- int01 + S0(time1[i])*dF1(i)
  }
  WR<- int10/int01
  MW<-int10 + 1/2*(1-int10 - int01)
  list(WR=WR,MW=MW,int10=int10,int01=int01)
}
