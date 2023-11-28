\name{Deming}
\alias{Deming}
\alias{coef.Deming}
\alias{vcov.Deming}
\alias{getFixedEf.Deming}
\alias{predict.Deming}
\alias{summary.Deming}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit Deming regression.
}
\description{
Deming regression fit. Assume x and y variances are the same. Slightly modified from MethComp R package.
}
\usage{
Deming(x, y, vr = sdr^2, sdr = sqrt(vr), boot = TRUE, keep.boot = FALSE, 
    alpha = 0.05) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
tbd
}
  \item{y}{
tbdy
}
  \item{vr}{
tbdvr
}
  \item{sdr}{
tbdsdr
}
  \item{boot}{
tbdboot
}
  \item{keep.boot}{
tbdkeep.boot
}
  \item{alpha}{
tbdalpha
}
}

\examples{

\dontrun{
set.seed(1)
x=rnorm(100,0,1)
y=x+rnorm(100,0,.5)
x=x+rnorm(100,0,.5)
fit=Deming(x,y, boot=TRUE)
summary(fit)
plot(x,y)
abline(fit)
# compare with lm fit
fit.1=lm(y~x, data.frame(x,y))
summary(fit.1)
abline(fit.1, col=2)
}
}
