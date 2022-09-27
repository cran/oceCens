#' Format ordered composite endpoint.
#'
#' Usually called from within \code{\link{oceTest}}. Input data.frame
#' with one row for each individual and columns for
#' k time-to-event outcomes, k status variables, and a group variable.
#' Format output so that each individual has several rows representing
#' different intervals at risk. Returns a list with elements used for later
#' calculations.
#'
#'
#'
#' @param data data.frame name, must have variables with names listed in
#'  \code{oceTime}, \code{oceStatus}, \code{group}
#' @param oceTime character vector with ordered (primary is first) names of
#'  different time-to-event variables.
#' @param oceStatus character vector with ordered names of status
#'  (0=censored, 1=event) variables.
#' @param group name of group variable.
#' @param id name of ID variable, NULL creates integer IDs.
#' @param oceNames long names of ordered endpoints, NULL uses \code{oceTime}.
#' @param outputDataFrame logical, output a data.frame in the list, defaults to
#'   FALSE for speed in the bootstrap.
#'
#'
#' @export
#'
#' @return A list with the following elements:
#' \describe{
#'  \item{timeMatrix}{n by k matrix with input values for k time-to-event values
#'   for each of n individuals}
#'   \item{statusMatrix}{n by k matrix of status values}
#'   \item{Z}{n vector of group variable with elements either 0 or 1}
#'   \item{oceNames}{k vector of long oceNames (for plotting labels)}
#'   \item{id}{m vector of individual ids, one element for each
#'      interval, so m>n}
#'   \item{group}{m vector of group values, either 0 or 1}
#'  \item{status}{m vector of status for each interval}
#'  \item{START}{m vector, START of interval}
#'  \item{STOP}{m vector, end of interval}
#'  \item{TAU}{maximum of the time-to-event outcomes}
#'  \item{IZMatrix}{m by k matrix, with jth column an indicator of representing
#'     ordering score 'time' for the jth endpoint}
#'  \item{data}{a data.frame output if outputDataFrame=TRUE, with variables:
#'    id, group, status, START, STOP, IZ1,...,IZk (columns of IZMatrix)}
#'  }
#' @examples
#'  d.temp<-data.frame(T1=c(1,4,3,6),s1=c(0,0,1,0),T2=c(4,1,5,3),
#'    s2=c(1,0,0,1),z=c(0,0,1,1))
#'  d.temp
#'  x<-oceFormat(data=d.temp,oceTime=c("T1","T2"),oceStatus=c("s1","s2"),
#'    group="z",outputDataFrame=TRUE)
#'  # put time to second event starting at TAU
#'  x$TAU
#'  x$data
#'
oceFormat <-
  function(data,
           oceTime,
           oceStatus,
           group,
           id = NULL,
           oceNames = NULL,
           outputDataFrame=FALSE) {

    # k=number of outcome types, in order
    k <- length(oceTime)
    if (is.null(oceNames)) {
      oceNames <- oceTime
    }
    n <- nrow(data)
    if (is.null(id)) {
      PATID <- 1:n
    } else {
      PATID <- data[, id]
    }
    Z <- data[, group]

    time <- status <- matrix(NA, nrow(data), k)
    for (j in 1:k) {
      time[, j] <- data[, oceTime[j]]
      status[, j] <- data[, oceStatus[j]]
    }
    TAU <- max(time)
    SUBID <- NULL
    START <- NULL
    STOP <- NULL
    IE <- NULL
    ZZ <- NULL

    # initialize atRisk vector, number at risk, and start of interval
    # for kth outcome
    atRisk <- rep(TRUE, n)
    nAtRisk <- n
    beginInterval <- 0
    # see Follmann, et al, 2020 Stat in Med, p. 615 for definition of IZj
    # it is Z times an indicator that the time period from START to STOP
    # is for endpoint j
    IZ<- NULL
    for (j in 1:k) {

      START <- c(START, rep(beginInterval, nAtRisk))
      STOP <- c(STOP, beginInterval + time[atRisk, j])
      IE <- c(IE, status[atRisk, j])
      ZZ <- c(ZZ, Z[atRisk])
      SUBID <- c(SUBID, PATID[atRisk])

      IZ.iter<- matrix(0,nAtRisk,k)
      IZ.iter[,j]<- Z[atRisk]
      IZ<- rbind(IZ,IZ.iter)


      # update at risk indicator for next iteration
      atRisk <- atRisk & status[, j] == 0
      nAtRisk <- sum(atRisk)
      # update beginInterval for next iteration
      beginInterval <- beginInterval + TAU
    }


    outputList<- list(
      timeMatrix = time,
      statusMatrix = status,
      Z = Z,
      oceNames = oceNames,
      id = SUBID,
      group = ZZ,
      status = IE,
      START = START,
      STOP = STOP,
      TAU = TAU,
      IZMatrix=IZ
    )

    if (outputDataFrame){
      dimnames(IZ)[[2]]<- paste0("IZ",1:k)
      data<- data.frame(group=ZZ,id=SUBID, status=IE, START=START, STOP=STOP,IZ)
      outputList<- c(outputList,list(data=data))
    }

    return( outputList)
  }

