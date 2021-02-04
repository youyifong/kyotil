\name{marginal.risk}
\alias{marginal.risk}
\alias{marginal.risk.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute maringalized risk
}
\description{
Computes risk of disease as a function of marker s by marginalizing over a covariate vector Z.
}
\usage{
marginal.risk(fit.risk, fit.s, data, ss=NULL, weights=rep(1, nrow(data)), t=NULL)
marginal.risk.boot(formula, marker.name, data, B, ci.type="quantile", numCores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.risk}{
A regression object, e.g. with formula y~Z+marker. Need to support predict(fit.risk)
}
  \item{fit.s}{
A linear regression object, e.g. with formula marker~Z. Need to support predict(fit.risk) and summary(fit.s)$sigma. Make sure fit.risk and fit.s are fit using the same data (e.g. by taking care of missing data)
}
  \item{data}{
A data frame
}
  \item{ss}{
A vector of marker values
}
  \item{weights}{
Inverse prob sampling weight, optional
}
  \item{t}{
A number. The time at which distribution function will be assessed
}

  \item{formula}{
A formula for the risk model san the marker
}
  \item{marker.name}{
Character
}
  \item{B}{
Number of bootstrap iterations
}
  \item{ci.type}{
Bootstrap confidence interval type
}
  \item{numCores}{
Number of cores for parallel processing
}
}
\value{
A matrix of two columns named marker and prob.
}
\details{
See the vignette file for more details.
}
\references{
Gilbert, P. B., Gabriel, E. E., Miao, X., Li, X., Su, S. C., Parrino, J., Chan, I. S. (2014). Fold rise in Antibody Titers by gpELISA is an Excellent Correlate of Protection for a Herpes Zoster Vaccine, Demonstrated via the Vaccine Efficacy Curve. J. Infect. Dis, 210, 1573-1581.
}
\examples{


\dontrun{

marker.name="EIA.log10d14overd0"
data=subset(dat.wide, trt==1)
formula=y2~riskScore.mat.endpoint1+vacc2birth

dat.tmp=data[!is.na(data[[marker.name]]), ] # make sure there is no missing data
fit.risk=glm(update(formula, as.formula(paste0("~.+",marker.name))), dat.tmp, family=binomial)
fit.s=lm(update(formula, as.formula(paste0(marker.name,"~."))), dat.tmp) 
res=marginal.risk(fit.risk, fit.s, dat.tmp)

plot(prob~marker, res, type="l", xlab=marker.name, ylab="probability of endpoint 2")
abline(h=mean(data$y2), col="gray", lty=2)

}

}
