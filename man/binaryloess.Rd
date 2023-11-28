\name{binaryloess}
\alias{binaryloess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Using loess to Check Functional Form for Logistic Regression
}
\description{
This function plots a smoothed line of how the average value of Y changes with X in order to check functional form for logistic regression. 
}
\usage{
binaryloess(x, y, scale = c("logit", "linear"), span = 0.7, weights = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
tbdx
}
  \item{y}{
tbdy
}
  \item{scale}{
tbdscale
}
  \item{span}{
smoothing parameter, passed to loess. If less than 1, the neighbourhood includes proportion a of the points. If greater than 1, all points are used, with the maximum distance assumed to be a^(1/p) times the actual maximum distance for p explanatory variables.
\cr
Missing records are removed first.
}
  \item{weights}{
sampling weights, passed to loess
}
  \item{\dots}{
passed to plotting function
}
}
\details{
This function comes from Jonathan Bartlett ()https://thestatsgeek.com/2014/09/13/checking-functional-form-in-logistic-regression-using-loess/).
}
\examples{

set.seed(1234)
n <- 1000
x <- rnorm(n)
xb <- -2+x
pr <- exp(xb)/(1+exp(xb))
y=rbern(n, pr)

par(mfrow=c(1,2))
binaryloess(x, y, scale = "logit", span = 0.7, weights = NULL, ylab="logit(p)")
binaryloess(x, y, scale = "linear", span = 0.7, weights = NULL, ylab="prob")

}
