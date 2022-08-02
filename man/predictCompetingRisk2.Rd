\name{pcr2}
\alias{pcr2}
\alias{predictCompetingRisk2}
\title{
Cumulative Incidence Function (CIF) Under Competing Risk - Approach 2
}
\description{
Fits cause-specific Cox models to each cause to compute cumulative incidence function for the cause of interest under competing risk. Weights are allowed in the optional arguments.

}
\usage{
predictCompetingRisk2(formula, formula.all, data, t0, newdata = data, stype
         = 2, ctype = 2, ...)
pcr2(formula.list, data, t0, newdata=data, ...)
}
\arguments{
  \item{formula.list}{
list of formulae for cause-specific failures. Assume the first cause to be the cause of interest
}
  \item{data}{
data frame
}
  \item{t0}{
the time till which cumulative incidence function is computed
}
  \item{newdata}{
new data for making prediction, default to the data for fitting the models
}
  \item{\dots}{
optional arguments that are passed to coxph, the most import of which is weights
}
}
\details{
The difference between predict.competing.risk and predict.competing.risk.2 is that instead of fitting a model to the overall failure, a model is fit for each cause, including the cause of interest. 
The overall survival is computed by adding together the cumulative hazard from individual causes.

See ?predict.competing.risk for additional details. 



}
\value{
A vector of real numbers as the risk till t0 for each subject in newdata
}
\references{

Competing Risk Analysis
Columbia Public Health
https://www.publichealth.columbia.edu/research/population-health-methods/competing-risk-analysis

Introduction to the Analysis of Survival Data in the Presence of Competing Risks
Peter C Austin, Douglas S Lee, Jason P Fine
https://www.ahajournals.org/doi/full/10.1161/CIRCULATIONAHA.115.017719

riskRegression: Predicting the Risk of an Event using Cox Regression Models by Brice Ozenne, Anne Lyngholm Sørensen, Thomas Scheike, Christian Torp-Pedersen, Thomas Alexander Gerds
https://journal.r-project.org/archive/2017/RJ-2017-062/RJ-2017-062.pdf

}
\seealso{
  \code{\link{predict.competing.risk}}.
}
\examples{


library(survival) 

# prepare a dataset with competing risk
lung1=lung[order(lung$time),]
lung1$status=lung1$status-1
lung1$status[1:50]=2
with(lung1, table(status))
lung1$status.1=ifelse(lung1$status==1,1,0)
lung1$status.2=ifelse(lung1$status==2,1,0)
lung1$status.a=ifelse(lung1$status==0,0,1)
lung1$wt=rep(1, nrow(lung1))

t0=1000
formula.list=list(
    Surv(time, status.1) ~ sex,
    Surv(time, status.2) ~ age
)

cif.2=pcr2(formula.list, lung1, t0)

## results match those from riskRegression_2022.07.13 when weights are all 1
#library(riskRegression)
#fit2 <- CSC(list(Hist(time,status)~sex, Hist(time,status)~age), data=lung1, surv.type="hazard",cause=1)
#r2=c(predictRisk(fit2,cause=1,product.limit=F,newdata=lung1,times=t0))
#head(cbind(cif, r2))

## pcr can produce cif that is greater than 1, especially when the covariates for the overall failure model differs from the covariates for the cause-specific model
## because the survival prob for the overall failure for a subject can be arbitrarily high
#cif.1=pcr(Surv(time, status.1) ~ sex, Surv(time, status.a) ~ age, lung1, t0)
##
#fit3 <- CSC(list(Hist(time,status)~sex, Hist(time,status)~age), data=lung1, surv.type="surv",cause=1)
#r3=c(predictRisk(fit3,cause=1,product.limit=F,newdata=lung1,times=t0))


# dealing with weights
lung1$wt=c(rep(2,50), rep(1, nrow(lung1)-50))
cif=predict.competing.risk.2(formula.list, lung1, t0, weights=lung1$wt)


}
