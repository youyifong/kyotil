\name{iorw}
\alias{iorw}
\alias{print.iorw}
\title{
Causal Mediation Analysis of Cowling et al.
}
\description{
Estimate the total, direct, and indirect effects using IORW method (inverse odds ratio weighting) and compute 95% percentile bootstrap confidence intervals. Modified from code written by Ben Cowling and Wey Wen Lim.
}
\usage{
iorw(formula.effect, formula.mediators, data, family =
 NULL, nboot = 10000, numCores = 1, save.steps = FALSE,
 verbose = FALSE)

\method{print}{iorw}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula.effect}{
a formula object for the total and direct effect regression. The first term on the right is assumed to be the binary treatment/exposure variable.
}
  \item{formula.mediators}{
a formula object for logistic regression. It should be of the form: ~ mediation marker1 + mediation marker2.
}
  \item{data}{
a data frame.
}
  \item{family}{
if Cox regression, leave as NULL; otherwise, it will be passed to glm().
}
  \item{nboot}{
an integer. Number of bootstrap replicates.
}
  \item{numCores}{
an interger. Number of cores to use for parallel procesing.
}
  \item{save.steps}{
boolean. Whether or not to save the fits from the three steps and the weights.
}
  \item{x}{
Object of type iorw
}
  \item{verbose}{
boolean. 
}
  \item{...}{
Additional arguments passed to the print function. 
}
}
\details{
Code by Cowling and Lim was downloaded from https://datadryad.org/stash/dataset/doi:10.5061/dryad.cv37539
\cr
If a bootstrap replicate generates warnings during regression, NA will be returned for that replicate. The number of such occurrences is recorded in an attribute of boot.perc in the return value.
\cr
It does not handle sampling weights yet.

}
\value{
Point estimates and percentile bootstrap confidence intervals.
}
\references{
Cowling, B. J., Lim, W. W., Perera, R. A., Fang, V. J., Leung, G. M., Peiris, J. M., & Tchetgen Tchetgen, E. J. (2019). Influenza hemagglutination-inhibition antibody titer as a mediator of vaccine-induced protection for influenza B. Clinical Infectious Diseases, 68(10), 1713-1717.
\cr
Nguyen, Q. C., Osypuk, T. L., Schmidt, N. M., Glymour, M. M., & Tchetgen Tchetgen, E. J. (2015). Practical guidance for conducting mediation analysis with multiple mediators using inverse odds ratio weighting. American journal of epidemiology, 181(5), 349-356.
\cr
Tchetgen Tchetgen, E. J. (2013). Inverse odds ratio-weighted estimation for causal mediation analysis. Statistics in medicine, 32(26), 4567-4580.
\cr
Imai, K., Keele, L., & Tingley, D. (2010). A general approach to causal mediation analysis. Psychological methods, 15(4), 309.
}
\author{
Youyi Fong, based on code by Cowling and Lim
}
\examples{


#### Cox regression

# without adjusting for baseline markers
library(survival)
formula.effect=Surv(surv_time, flu)~vaccine+age
formula.mediators=~log2(postvax.B.Brisbane/5)
res.1=iorw(formula.effect, formula.mediators, kid, nboot=10, numCores=1); res.1
stopifnot(max(abs(res.1$boot[1,] - c(0.2029779,0.6070105,0.3039110,0.4283389,0.2124268)))<1e-6)

# adjust for baseline markers
formula.effect=Surv(surv_time, flu)~vaccine+log2(prevax.B.Brisbane)+age
formula.mediators=~log2(postvax.B.Brisbane/5)
res.2=iorw(formula.effect, formula.mediators, kid, nboot=10, numCores=1); res.2


#### Logistic regression

# without adjusting for baseline markers
formula.effect=flu~vaccine+age
formula.mediators=~log2(postvax.B.Brisbane/5)
res.3=iorw(formula.effect, formula.mediators, kid, family=binomial(), nboot=10, numCores=1); res.3
stopifnot(max(abs(res.3$boot[1,] - c(0.1960024,0.6154349,0.2937164,0.4145470,0.2168644)))<1e-6)

# adjust for baseline markers
formula.effect=flu~vaccine+log2(prevax.B.Brisbane)+age
formula.mediators=~log2(postvax.B.Brisbane/5)
res.4=iorw(formula.effect, formula.mediators, kid, family=binomial(), nboot=10, numCores=1); res.4






}
