\name{cox.zph.2}
\alias{cox.zph.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test the Proportional Hazards Assumption of a Cox Regression (a slightly modified version)
}
\description{
A slightly modified test of the proportional hazards assumption for a Cox regression model fit (coxph). This version corrects some conservativeness of the test.
}
\usage{
cox.zph.2(fit, transform = "km", global = TRUE, exact=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
fit
}
  \item{transform}{
transform
}
  \item{global}{
global
}
  \item{exact}{
Boolean. If FALSE, this function is an identical copy of cox.zph. If TRUE, it computes the variance of the test statistic exactly, instead of approximately.
}
}
\details{
When the model uses time-dependent covariates, the approximation used in Grambsch and Therneau resulted in conservativeness of the test. This is "fixed" here at a cost of up to 2.5 times longer execution time. 
}
\references{
Fong, Y. and Halloran, M Elizabeth and Gilbert, P.  Using Time-Dependent Age Group in Cox Regression Analysis of Vaccine Efficacy Trials, Just Another Epi Journal, in prep.
}
\seealso{
\code{\link{cox.zph}}
}
\examples{
library(survival)
fit <- coxph(Surv(futime, fustat) ~ age + ecog.ps,  
             data=ovarian) 
temp <- cox.zph(fit) 
print(temp)        
temp.2 <- cox.zph.2(fit) 
print(temp.2)        

}
