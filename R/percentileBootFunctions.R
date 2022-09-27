#' Percentile Bootstrap Two-sided Confidence Intervals and p-values
#'
#' Input vector of bootstrap replicates and get either the two-sided percentile
#' confidence interval or the compatible two-sided p-value.
#'
#' @param Ti A numeric vector of bootstrap replicates of an estimate.
#' @param conf.level Confidence level.
#' @param theta0 Null hypothesis value of estimand.
#'
#' @return \code{percci} returns only a two-sided confidence interval and
#'  \code{percpval} returns only a two-sided p-value.
#'
#' @importFrom stats median
#' @export
#'
#' @details
#' Simple functions, where \code{percci} gives two-sided confidence intevals
#' and \code{percpval} gives two-sided p-values.
#'
#' We get a two-sided p-value by inverting the percentile Bootstrap
#' confidence interval. This is not straightforward if there are not enough
#' bootstrap samples and/or if the minimum and maximum of the replicates do not
#' cover the null value. If there are B bootstrap resamples, then the interval
#' from the minimum to the maximum has confidence level =1- 2/(B+1).
#' We can see this because the percentile interval
#' (see Efron and Tibshirani, 1993, p. 160 bottom) is
#' T[k], T[B+1-k]   where k=floor( (B+1)*(1-conf.level)/2),
#' where T is an ordered vector of B test statistics calculated from B
#' bootstrap replicates (T=Ti[order(Ti)]).
#' Therefore, if conf.level > 1 - 2/(B+1) then we cannot get a percentile
#' interval, so if the min and max of T do not surround theta0, then
#' a two-sided p-value can be stated to be p<= 2/(B+1). If the p-value
#' is 2/(B+1), then it is the lowest possible for that B, and increasing
#' B may produce a lower p-value.
#'
#' @references
#'  Efron, B and Tibshirani, RJ (1993) An Introduction to the Bootstrap.
#'  Chapman and Hall.
#'
#' @examples
#' \donttest{
#'   set.seed(123)
#'   y<- rnorm(100)+0.1
#'   nB<- 1e5
#'   Tstat<- rep(NA,nB)
#'   for (i in 1:nB){
#'     Tstat[i]<-mean( sample(y,replace=TRUE) )
#'    }
#'    # two-sided bootstrap percentile p-value
#'    # that mean is different from 0
#'    percpval(Tstat,theta0=0)
#'    # 95% percentile interval
#'    percci(Tstat)
#'    # compare to t-test
#'    t.test(y)
#'
#'    # to show that the functions are close to compatiable
#'    # set confidence level to 1-pvalue
#'    pval<-percpval(Tstat,theta0=0)
#'    confLevel<- 1-pval
#'    pval
#'    # then lower limit should be close to 0
#'    percci(Tstat, conf.level=confLevel)
#'   }
percci<-function(Ti,conf.level=.95){
  ### see Efron and Tibshirani, p. 160 bottom

  alpha<- (1-conf.level)/2
  B<-length(Ti)
  k<-floor((B+1)*alpha)
  if (k==0){
    warning("increase number of bootstrap samples")
    ci<-c(-Inf,Inf)
  } else {
    oTi<-Ti[order(Ti)]
    ci<-oTi[c(k,B+1-k)]
  }
  attr(ci,"conf.level")<- conf.level
  ci
}
#' @describeIn percci Bootstrap percentile p-values
#' @export
#' @md
percpval<-function(Ti,theta0=0){
  ### obtain a two-sided p-value by inverting the percentile Bootstrap
  ### ci. If there are B bootstrap resamples, then the interval from the
  ### minimum to the maximum has confidence level =1- 2/(B+1)
  ### We can see this because the percentile interval
  ### (see Efron and Tib, p. 160 bottom)
  ### is
  ### T(k), T(B+1-k)   where k=floor( (B+1)*(1-conf.level)/2)
  ### so if conf.level > 1 - 2/(B+1) then we cannot get a percentile interval
  ### so if the min and max of Ti do not surround 0, then
  ### a two-sided p-value can be stated to be p< 2/(B+1)
  B<-length(Ti)
  minp<- 2/(B+1)
  if ( min(Ti) >=theta0 | max(Ti) <= theta0 ){
    p<- minp
  } else {
    ### define k such that
    ### T(k-1) < theta0 <= T(k)
    ### if theta0<=median
    ### T(B+1-k) > theta0 > T(B+2-k)
    ### if theta0> median
    if (theta0<=median(Ti)){
      k<-length(Ti[Ti<theta0]) + 1
    } else {
      k<-length(Ti[Ti>theta0])+1
    }
    p<- min(1,2*k/(B+1))
  }
  p
}
