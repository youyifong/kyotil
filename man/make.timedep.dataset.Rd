\name{make.timedep.dataset}
\alias{make.timedep.dataset}
\title{
Create Dataset for Time-dependent Covariate Proportional Hazard Model Analaysi
}
\description{
Returns a data frame that is suitable for time-dependent covariate Cox model fit. 
}
\usage{
make.timedep.dataset(dat, X, d, baseline.ageyrs, t.1, t.2 = NULL)
}
\arguments{
  \item{dat}{
data frame
}
  \item{X}{
string. Name of the followup time column in dat. Unit needs to be years.
}
  \item{d}{
string. Name of the followup time column in dat.
}
  \item{baseline.ageyrs}{
string. Name of the followup time column in dat.
}
  \item{t.1}{
numerical. Cutoff for age group
}
  \item{t.2}{
numerical. Second cutoff for age group 
}
}
\details{
The function assumes that the followup length is such that only one change of age group is possible.

}
\value{
Returns a data frame with the following columns added:  tstart, tstop, .timedep.agegrp, .baseline.agegrp
  \item{tstart}{left bound of time interval}
  \item{tstop}{right bound of time interval}
  \item{.timedep.agegrp}{time-dependent age group}
  \item{.baseline.agegrp}{baseline age group}
}
\references{
Therneau, T. and Crowson, C. Using Time Dependent Covariates and Time Dependent Coefficients in the Cox Model. A vignette from the R package surival.
}
\author{
Youyi Fong
}

\examples{

library(survival)

n=3000; followup.length=5; incidence.density=0.015; age.sim="continuous"

dat.0=sim.dat.tvarying.two(n, followup.length, incidence.density, age.sim, seed=1)
dat=subset(dat.0, for.non.tvarying.ana, select=c(ptid, X, d, baseline.age, trt))
dat.timedep = make.timedep.dataset (dat, "X", "d", "baseline.age", 6)
coxph(Surv(tstart,tstop,d) ~ trt*.timedep.agegrp, dat.timedep)


}
\keyword{ time varying }
