\name{pcr}
\alias{pcr}
\alias{predict.competing.risk}
\title{
Cumulative Incidence Function (CIF) Under Competing Risk
}
\description{
Computes cumulative incidence function under competing risk.
}
\usage{
pcr(formula, formula.all, data, t0, newdata, ...)
}
\arguments{
  \item{formula}{
formula for the cause-specific failure
}
  \item{formula.all}{
formula for all-cause failure
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
  \item{stype}{
computation of the survival curve, 1=direct, 2= exponenial of the cumulative hazard. Default 2, which is the default of basehaz and predict.coxph
}
  \item{ctype}{
whether the cumulative hazard computation should have a correction for ties, 1=no, 2=yes. Default 2, which is the default of basehaz and predict.coxph
}
  \item{\dots}{
optional arguments that are passed to coxph, the most import of which is weights
}
}
\details{
Fits a cause-specific Cox model and a all-cause Cox model. Weights are allowed in the optional arguments.

When there is only one cause, CIF is conceptually 1 - surival prob.

}
\value{
A vector of real numbers
}
\references{
https://www.publichealth.columbia.edu/research/population-health-methods/competing-risk-analysis
}
\examples{


library(survival) 

lung1=lung[order(lung$time),]
lung1$status=lung1$status-1
with(lung1, table(status))
lung1$status.1=ifelse(lung1$status==1,1,0)
lung1$status.a=ifelse(lung1$status==0,0,1)
lung1$wt=rep(1, nrow(lung1))


# when there is no covariate and one cause, CIF = 1 - KM estimate of survival prob
# stype=2 is surv=prod limit
fitKM <- survfit(Surv(time, status.1) ~ 1, data=lung1, stype=1, ctype=2)
cbind(summary(fitKM)$cumhaz, exp(-summary(fitKM)$cumhaz), summary(fitKM)$surv)[1:2,]
#[1,] 0.004385965 0.9956236 0.9956140
#[2,] 0.017660474 0.9824946 0.9824561
cif=predict.competing.risk(Surv(time, status.1) ~ 1, Surv(time, status.1) ~ 1, lung1, t0=11, newdata=lung1[1,,drop=F], weights=lung1$wt, stype=1, ctype=2)
cif # 0.01754386
1-cif # 0.9824561 = summary(fitKM)$surv at t=11


# when there are covariates and one cause, CIF and 1-exp(-H) are close to each other when H is small but not close when H is large
form  =Surv(time, status.1) ~ age
form.a=Surv(time, status.a) ~ age
par(mfrow=c(1,2))
for (t0 in c(12,1000)) {
    fit=coxph(form, lung1, weights=lung1$wt)
    lung2=lung1; lung2$time=t0
    r=predict(fit, type="expected", newdata=lung2)
    print(head(basehaz(fit, centered=T)))
    cif=predict.competing.risk(form, form.a, lung1, t0, newdata=lung1, weights=lung1$wt, stype=2, ctype=2)
    plot(cif, 1-exp(-r), xlab="Cumulative incidence function", ylab="Expected number of events from predict.coxph", main="t0: a later time point"); abline(0,1)
    print(head(cbind(cif, 1-exp(-r))))
}



## results match those from riskRegression_2022.07.13 when weights are all 1
#lung1=lung[order(lung$time),]
#lung1$status=lung1$status-1
#lung1$status[1:50]=2
#with(lung1, table(status))
#lung1$status.1=ifelse(lung1$status==1,1,0)
#lung1$status.a=ifelse(lung1$status==0,0,1)
#lung1$wt=rep(1, nrow(lung1))
#t0=1000
#library(riskRegression)
#fit2 <- CSC(Hist(time,status)~age, data=lung1, surv.type="survival",cause=1)
#r2=c(predictRisk(fit2,cause=1,product.limit=F,newdata=lung1,times=t0))
##r3=c(predictRisk(fit2,cause=1,product.limit=T,newdata=lung1,times=t0))# slightly different from r2
#cif=predict.competing.risk(form, form.a, lung1, t0, newdata=lung1, weights=lung1$wt, stype=2, ctype=2)
#head(cbind(cif, r2))


}
