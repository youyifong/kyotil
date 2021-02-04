\name{VEplot}
\alias{VEplot}
\alias{VEplot.cox.zph}
\alias{VEplot.glm}
\alias{myplot.cox.zph}
\title{
Vaccine Efficacy Plots}
\description{
Vaccine efficacy plots. 
}
\usage{

VEplot (object, ...) 
 
\method{VEplot}{cox.zph} (object, resid = TRUE, se = TRUE, df = 4, nsmo = 40, 
    var, ylab="VE", xlab="Time", xaxt="s", cex.axis=1, ...) 

\method{VEplot}{glm}(object, X1, X2, x, ...)

\method{myplot}{cox.zph} (object, resid = TRUE, se = TRUE, df = 4, nsmo = 40, var, 
    coef.transform=NULL, 
    ylab=NULL, 
    xlab="Time", xaxt="s", cex.axis=1, 
    ...) 

}

\arguments{

  \item{object}{An object}
  \item{resid}{Boolean, whether to plot residuals}
  \item{se}{Boolean, whether to plot confidence band}
  \item{df}{degrees of freedom}
  \item{nsmo}{number of points used to plot the fitted spline}
  \item{var}{estimated variance matrix from the Cox model fit}
  \item{xlab}{x label}
  \item{xaxt}{x axis}
  \item{cex.axis}{cex for axis}
  \item{ylab}{y label}
  \item{coef.transform}{a function to transform Cox hazard ratio estimate}
  \item{X1}{a matrix of dimension k by p, where k is the length of x (see below) and p is the length of coef(object)}
  \item{X2}{a matrix of dimension k by p, where k is the length of x (see below) and p is the length of coef(object)}
  \item{x}{a vector of length k that represents the x coordinate of the VE plot}
  \item{...}{additional parameters}
  
}

\details{
VEplot and myplot.cox.zph are extensions of survival::plot.cox.zph to plot VE curve and other transformations.

myplot.cox.zph adds the following parameters to the original list of parameters in plot.cox.zph:
    coef.transform: a function to transform the coefficients
    ylab: y axis label 
    xlab: x axis label
    
VEplot.glm computes a series of k VEs: for i in 1...k, VE[i] = P(Y=1|X1[i,])/P(Y=1|X2[i,]). It returns a 3 by k matrix, whose first row contains VE estimates and the second and third rows contain lower and upper bounds, respectively.
    
}
\author{
Youyi Fong, Dennis Chao
}
\references{
Durham, Longini, Halloran, Clemens, Azhar and Rao (1998) "Estimation of vaccine efficacy in the presence of waning: application to cholera vaccines." American Journal of Epidemiology 147(10): 948-959.
}
\examples{

library(survival)
vfit <- coxph(Surv(time,status) ~ trt + factor(celltype) + 
              karno + age, data=veteran, x=TRUE) 
temp <- cox.zph(vfit) 

par(mfrow=c(2,2))
for (v in c("trt","age")) {
    VEplot(temp, var=v, resid=FALSE, main=v, ylab="VE", cex.axis=1.5)
    plot(temp, var=v, resid=FALSE, main=v)
}

library(survival)
fit <- glm(status ~ trt + trt*age, data=veteran) 
summary(fit)
age=seq(min(veteran$age),max(veteran$age),length=10)
out = VEplot(fit, X1=cbind(1,1,age,1*age), X2=cbind(1,0,age,0*age), x=age)
out

}
