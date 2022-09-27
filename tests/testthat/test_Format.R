context("simple cases: oceFormat")

d<-data.frame(T1=c(1,4,3,6),s1=c(0,0,1,0),T2=c(4,1,5,3),s2=c(1,0,0,1),z=c(0,0,1,1))

test_that("Simple case: 2 times, 4 individuals",{

  expect_equal(
    oceFormat(oceTime=c("T1","T2"),
              oceStatus=c("s1","s2"),
              group="z",
              id = NULL,
              oceNames = NULL,
              data=d)$STOP,
    c(1,  4,  3,  6, 10,  7,  9)
  )
})
