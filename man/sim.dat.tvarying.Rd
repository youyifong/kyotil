\name{sim.dat.tvarying.two}
\alias{sim.dat.tvarying.two}
\alias{sim.dat.tvarying.three}
\title{
Simulation Functions for Time-dependent Proportional Hazard Model
}
\description{
sim.dat.tvarying.three simulates from a model with time varing age group variale of three levels, sim.dat.tvarying.two two.
}
\usage{
sim.dat.tvarying.three(n, followup.length, incidence.density, 
    age.sim = c("tvaryinggroup", "baselinegroup", "continuous","bt"),
    random.censoring.rate = 0.05, seed)
    
sim.dat.tvarying.two(n, followup.length, incidence.density, 
    age.sim = c("tvaryinggroup", "baselinegroup", "continuous","bt"), 
    random.censoring.rate = 0.05, seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
integer. Sample size.}
  \item{followup.length}{
numeric. Length of followup, in years.
}
  \item{incidence.density}{
numeric. Incidence rate per year.
}
  \item{age.sim}{
string. Choose between one of three possibilities. tvaryinggroup: age group is time-varying covariate; baselinegroup: age group is a baseline covariate; continuous: age is a continuous covariate; bt: age group by treatment interaction uses baseline age group, while age group main effect uses time-dependent age group
}
  \item{random.censoring.rate}{
numeric. Amount of random censoring.
}
  \item{seed}{
integer. Random number generator seed.
}
}
\details{

In sim.dat.tvarying.three, baseline age is uniformly distributed between 2.0 and 16.0, and divivded into three groups at 6 and 12.
In sim.dat.tvarying.two, baseline age is uniformly distributed between 2.0 and 12.0, and divivded into two groups at 6.

}
\value{
Return a data frame with the following columns:
    \item{\code{ptid}}{subject identifier}
    \item{\code{trt}}{treatment indicator 0/1}
    \item{\code{for.non.tvarying.ana}}{Boolean, used to subset dataset for non-time dependent analysis}
    \item{\code{C}}{censoring time}
    \item{\code{baseline.age}}{age years at baseline}
    \item{\code{agegrp}}{a factor with levels \code{[0,6)} \code{[6,12)} \code{[12,100)}}
    \item{\code{baseline.agegrp}}{a factor with levels \code{[0,6)} \code{[6,12)} \code{[12,100)}}
    \item{\code{tstart}}{left bound of time interval}
    \item{\code{tstop}}{right bound of time interval}
    \item{\code{d}}{event indicator}
    \item{\code{X}}{followup time, in years}
}
\author{
Youyi Fong
}

\seealso{
\code{\link{make.timedep.dataset}}
}
\examples{

library(survival)

dat=sim.dat.tvarying.three(n=6000,followup.length=3, incidence.density=0.05, 
    age.sim="tvaryinggroup", seed=1)
f.tvarying = Surv(tstart,tstop,d) ~ trt*agegrp 
f =          Surv(X,d)            ~ trt*baseline.agegrp 
fits=list()
fits[["tvarying"]]=coxph(f.tvarying, dat)
fits[["baseline"]]=coxph(f, subset(dat, for.non.tvarying.ana))
fits







}
