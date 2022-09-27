## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=TRUE---------------------------------------------------------------
library(oceCens)
data(simScenario5)
head(simScenario5)

## ---- eval=TRUE---------------------------------------------------------------
oceTest(data=simScenario5, oceTime=c("T1","T2","T3"),
 oceStatus=c("I1","I2","I3"), group=c("Z"), id = "PATID",
 oceNames = c("Death","Stroke/MI","Bleed"))

## ---- eval=TRUE---------------------------------------------------------------
oceTest(data=simScenario5, oceTime=c("T1","T2","T3"),
 oceStatus=c("I1","I2","I3"), group=c("Z"), id = "PATID",
 oceNames = c("Death","Stroke/MI","Bleed"), conf.int=TRUE, method=c("coxph"),ciMethod="WLW")

## -----------------------------------------------------------------------------
xform<-oceFormat(data=simScenario5,oceTime=c("T1","T2","T3"),
    oceStatus=c("I1","I2","I3"),
    group="Z",outputDataFrame=TRUE)

## -----------------------------------------------------------------------------
# perform cox regression using time varying treatment efects, IZ1,IZ2, IZ3
 # associated with the 3 prioritized endpoints
cout<- coxph(Surv(START, STOP, status) ~ IZ1+IZ2+IZ3+cluster(id), data=xform$data)

## -----------------------------------------------------------------------------
cout

## -----------------------------------------------------------------------------
coxph2WR(cout)

## -----------------------------------------------------------------------------
# create some made up data first...
# demographics has 400 rows, one row for each id
demographics<- data.frame(id=1:400,
                     age=sample(55:90,400,replace=TRUE),                  
                     sex=sample(c("M","F"),400,replace=TRUE))
# the formatted data has several rows for each id
dim(xform$data)
# when you merge, use all.x=TRUE to fill in the
# demographics for each  individual in the proper place.
simScenario5.wDemo<- merge(xform$data,demographics,by="id",all.x=TRUE)
dim(simScenario5.wDemo)
head(simScenario5.wDemo)

## -----------------------------------------------------------------------------
cout<- coxph(Surv(START, STOP, status) ~         IZ1+IZ2+IZ3+age+sex+cluster(id), 
             data=simScenario5.wDemo)
cout

## -----------------------------------------------------------------------------
coxph2WR(cout)

## ----fig1, fig.width=7, fig.height=7------------------------------------------
dataFormt<-oceFormat(data=simScenario5,
   oceTime=c("T1","T2","T3"),
   oceStatus=c("I1","I2","I3"), group=c("Z"),
   oceNames = c("Death","Stroke/MI","Bleed"))
npmleOutput<- oceNPMLE(dataFormt)
plot(npmleOutput, xlab="Custom x label", mark.time=FALSE, lwd=2)
#  add lines from coxph output
coxOutput<- oceCoxph(dataFormt)
plot(coxOutput,linesonly=TRUE, col=c("orange","purple"),lwd=2)
legend("bottomleft",
   legend=c("grp=0, NPMLE","grp=1, NPMLE",
            "grp=0, coxph","grp=1, coxph"),
   col=c("red","blue","orange","purple"),lty=c(1,1,1,1),lwd=2)

