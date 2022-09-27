#' Tests for ordered composite endpoints with censoring.
#'
#' An ordered composite endpoint (oce) is a way of ranking responses by
#' ordering several types of responses by order of importance. Rank by
#' the most important response, then break ties with the next most important,
#' and so on. The tests here are based on two sample tests. Let Y0 and Y1
#' be the oce score in the control arm and treatment arm, respectively. Then
#' here we estimate both the win ratio (WR), P[Y1>Y0]/P[Y0>Y1], or the
#' Mann-Whitney parameter, P[Y1>Y0] + (1/2) Pr[Y1=Y0]. Different methods are
#' used to estimate those parameters, and inferences are done by bootstrap
#' percentile methods.
#'
#' @inheritParams oceFormat
#' @param method Estimation method, one of 'all', 'npmle', 'coxph' or 'simple'.
#'   Default is 'all' which calculates all of the three methods. See details.
#' @param ciMethod confidence interval method, default is 'bootstrap'
#' @param conf.int Logical, should confidence intervals be calculated.
#' @param conf.level confidence level.
#' @param nBoot number of bootstrap replicates (ignored if conf.int=FALSE).
#' @param plot logical, plot oce score by group as survival functions
#'    (NPMLE version, except if method='coxph'). For more control over those
#'    plots see either \code{\link{plot.oceNPMLE}} or
#'    \code{\link{plot.oceCoxph}}.
#' @param ... holder space for future arguments.
#'
#' @details
#' This idea is to stack the time to first event for the k different types of
#' events. So if TAU is the maximum time that any individual is in the
#' study, then the primary type of event has scores that fall into (0,TAU],
#' the secondary type has scores that fall into (TAU,2*TAU], and so on.
#' Then we rank by the primary type (e.g., death), but if there are many ties
#' in the primary type (e.g., many that did not die during the study), then
#' we break ties by the secondary type of event, and so on.
#'
#' The difficulty is when there is censoring in time, because that imposes
#' interval censoring on the score scale. This can be handled with interval
#' censoring methods (although in a non-standard way). The 'npmle' method
#' calculates a nonparametric maximum likelihood estimate of the 'survival'
#' distribution of the ordering score for each arm, then gets the estimates
#' by numeric integration. The 'coxph' method uses an interval censored
#' proportional hazards model treating the oce scores as time
#' using \code{\link[survival]{coxph}} from the
#' survival R package. The 'simple' method uses part of the 'coxph'
#' method together with a more simple estimator. Each method produces
#' a win ratio (P[Y1>Y0]/P[Y0>Y1]) and a Mann-Whitney
#' (P[Y1>Y0] + (1/2) Pr[Y1=Y0]) estimate. Details are given in
#' Follmann, et al (2020).
#'
#' When \code{ciMethod="bootstrap"} inferences are done by nonparametric
#' bootstrap percentile method (see
#' \code{\link{percci}}) in order to account for the correlation among the
#' different types of responses. When \code{ciMethod="WLW"} and
#' \code{method="coxph"}, then the win ratio is calculated by the Cox model
#' with the standard errors of the log(HR) or log(WR) calculated by the robust
#' sandwich method suggested by Wei, Lin, and Weissfeld (1989).
#' P-values are all two-sided and test the
#' null hypothesis of no difference between the arms (for the win ratio, the
#' null value is 1, while for the MW the null value is 0).
#'
#' For access to the \code{coxph} output see \code{\link{oceCoxph}}, or for the
#' NPMLE output see \code{\link{oceNPMLE}}.
#'
#' @return If \code{conf.int=FALSE} then a vector of estimates determined
#'  by \code{method} results. If \code{conf.int=TRUE} then a matrix is returned
#'  with a row for each estimate, and 4 columns for the Estimate, lower
#'  confidence limit, upper confidence limit, and two-sided p-value.
#'
#' @export
#'
#' @references
#' Follmann, D., Fay, M. P., Hamasaki, T., and Evans, S. (2020). Analysis of
#' ordered composite endpoints. Statistics in Medicine, 39(5), 602-616.
#'
#' Wei, L. J., Lin, D. Y., & Weissfeld, L. (1989). Regression analysis of
#' multivariate incomplete failure time data by modeling marginal distributions.
#' Journal of the American statistical association, 84(408), 1065-1073.
#'
#' @examples
#' data(simScenario5)
#' oceTest(data=simScenario5, oceTime=c("T1","T2","T3"),
#'  oceStatus=c("I1","I2","I3"), group=c("Z"), id = "PATID",
#'  oceNames = c("Death","Stroke/MI","Bleed"), method=c("all"))
oceTest<-function(data,
                 oceTime,
                 oceStatus,
                 group,
                 id = NULL,
                 oceNames = NULL,
                 method=c("all","npmle","coxph","simple"),
                 ciMethod=c("WLW","bootstrap"),
                 conf.int=FALSE,
                 conf.level=0.95,
                 nBoot=2000,plot=FALSE,...){
  # check that group variable is 0 or 1
  if (!all(sort(unique(data[,group]))==c(0,1))) stop("group variable must have only values of 0 and 1")

  if (conf.int){
    # get starting time to estimate how long the bootstrap will take
    comp.time0<- proc.time()
  }
  xform<- oceFormat(data=data,
                    oceTime=oceTime,
                    oceStatus=oceStatus,
                    group=group,
                    id = id,
                    oceNames = oceNames)

  method<- match.arg(method)
  ciMethod<- match.arg(ciMethod)

  if (method=="all"){
   outNPMLE<- oceNPMLE(xform)
   outCox<- oceCoxph(xform)
   outSimple<- oceSimple(xform,outCox)
   output<- c(WR.npmle=outNPMLE$WR,
              WR.coxph=outCox$WR,
              WR.simple=outSimple$WR,
              MW.npmle=outNPMLE$MW,
              MW.coxph=outCox$MW,
              MW.simple=outSimple$MW)
  } else if (method=="npmle"){
    outNPMLE<- oceNPMLE(xform)
    output<- c(WR.npmle=outNPMLE$WR,
               MW.npmle=outNPMLE$MW)
  } else if (method=="coxph"){
    if (conf.int & ciMethod=="WLW"){
      # run analyses later
    } else {
      outCox<- oceCoxph(xform)
      output<- c(WR.coxph=outCox$WR,
                 MW.coxph=outCox$MW)
    }

  } else if (method=="simple"){
    # although we need to run the oceCoxph
    # to get the results for oceSimple
    # do not output results from coxph
    # to avoid appearance/temptation of
    # running both and picking the
    # largest effect
    outCox<- oceCoxph(xform)
    outSimple<- oceSimple(xform,outCox)
    output<- c(WR.simple=outSimple$WR,
               MW.simple=outSimple$MW)
  }

  if (plot){
    if (method !="coxph"){
      plot(outNPMLE)
    } else {
      plot(outCox)
    }
  }
  if (conf.int & ciMethod=="bootstrap"){

     B<- nBoot
     # use recursive programming, but make sure to set cont.int=FALSE
     T0<- output
     Tnames<- names(T0)
     m<- length(T0)
     Tmatrix<- matrix(NA,B, m,dimnames=list(NULL,Tnames))
     I<- 1:nrow(data)
     # to get a good estimate of the computation time
     # run the first 9 reps
     if (B<10) stop("nBoot must be at least 10")
     for (i in 1:9){
       data.i<- data[sample(I,replace=TRUE), ]
       Tmatrix[i,]<- oceTest(data=data.i,
                             oceTime=oceTime,
                             oceStatus=oceStatus,
                             group=group,
                             id = id,
                             oceNames = oceNames,
                             method=method,
                             conf.int=FALSE)
     }
     # estimate how much time the bootstrapping will take,
     # this is a rough estimate
     comp.time1<- proc.time()
     time.for.10reps<- comp.time1["elapsed"] - comp.time0["elapsed"]
     run.time.seconds<- round(time.for.10reps*B/10,1)
     if (run.time.seconds<61){ runTime<- paste0(run.time.seconds, " seconds....")
     } else {
       runTime<- paste0(round(run.time.seconds/60,1), " minutes....")
     }
     MESSAGE<- paste0("Running ",B," bootstrap replicates. Rough estimate of time until completion: ",runTime)
     message(MESSAGE)
     for (i in 10:B){
       data.i<- data[sample(I,replace=TRUE), ]
       Tmatrix[i,]<- oceTest(data=data.i,
                             oceTime=oceTime,
                             oceStatus=oceStatus,
                             group=group,
                             id = id,
                             oceNames = oceNames,
                             method=method,
                             conf.int=FALSE)
     }
     lo<- paste0("lower ",round(100*conf.level,2),"% CL")
     hi<- paste0("upper ",round(100*conf.level,2),"% CL")
     OUT<- matrix(NA,m,4,dimnames=list(Tnames,c("Estimate",lo,hi,"two-sided p-value")))
     OUT[,1]<-T0



     for (j in 1:m){
       ci<- percci(Tmatrix[,j],conf.level=conf.level)
       OUT[j,2:3]<-ci
       # Null hypothesis value for win ratio=1
       if (grepl("WR",Tnames[j])){
         THETA0<-1
       } else {
         # null hypothesis value for MW is 1/2
         THETA0<- 0.5
       }
       OUT[j,4]<- percpval(Tmatrix[,j], theta0=THETA0)
     }
     output<- OUT
  } else if (conf.int & ciMethod=="WLW"){
    if (method !="coxph") stop("ciMethod='WLW' only available when method='coxph' ")
    # code is a little inefficient, because we run coxph twice.
    IE <- xform$status
    Z <- xform$group
    START <- xform$START
    STOP <- xform$STOP
    ID<- xform$id
    regord <- coxph(Surv(START, STOP, IE) ~ Z +cluster(ID))
    output<-  coxph2WR(regord,conf.level=conf.level)
  }

  output
}
