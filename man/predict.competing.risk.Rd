\name{pcr}
\alias{pcr}
\alias{predict.competing.risk}
\title{
Cumulative Incidence Function Under Competing Risk
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

}
\value{
A vector of real numbers
}
\references{
https://www.publichealth.columbia.edu/research/population-health-methods/competing-risk-analysis
}
\examples{


library(survival) 

# make a second status
with(lung, table(status))

lung1=lung[order(lung$time),]

lung1$status=lung1$status-1
lung1$status[1:50]=2
with(lung1, table(status))
lung1$status.1=ifelse(lung1$status==1,1,0)
lung1$status.a=ifelse(lung1$status==0,0,1)
lung1$wt=rep(1, nrow(lung1))
#lung1$wt=c(rep(1,100), rep(2,nrow(lung1)-100))


form  =Surv(time, status.1) ~ age
form.a=Surv(time, status.a) ~ age

par(mfrow=c(1,2))
fitKM <- survfit(form, data=lung1, stype=1, ctype=1)
plot(exp(-summary(fitKM)$cumhaz), summary(fitKM)$surv); abline(0,1)
summary(fitKM)$cumhaz[1:10]
summary(fitKM)$surv[1:10]
fitKM <- survfit(form, data=lung1, stype=2, ctype=1)
plot(exp(-summary(fitKM)$cumhaz), summary(fitKM)$surv); abline(0,1)


## results match those from riskRegression_2022.07.13 when there are no weights
#t0=12
#library(riskRegression)
#fit2 <- CSC(Hist(time,status)~age, data=lung1, surv.type="survival",cause=1)
#r2=c(predictRisk(fit2,cause=1,product.limit=F,newdata=lung1,times=t0))
##r3=c(predictRisk(fit2,cause=1,product.limit=T,newdata=lung1,times=t0))# slightly different from r2
#cif=predict.competing.risk(form, form.a, lung1, t0, newdata=lung1, weights=lung1$wt, stype=2, ctype=2)
#cbind(cif, r2)


t0=12


# cif does not match r when there is only 1 cause, there is 1% shift
fit=coxph(form, lung1, weights=lung1$wt)
lung2=lung1; lung2$time=t0
r=predict(fit, type="expected", newdata=lung2, stype=1, ctype=1)

#exp(coef(fit)*lung1[1:2,"age"]) * basehaz(fit)[1,1]

cif=predict.competing.risk(form, form.a, lung1, t0, newdata=lung1, weights=lung1$wt, stype=2, ctype=2)#; 1-cif[1:3]
head(cbind(cif, 1-exp(-r)))
plot(cif, 1-exp(-r), xlab="Cumulative incidence function", ylab="1-Survival from predict.coxph"); abline(0,1)

#formula=form; formula.all=form.a; data=lung1; newdata=data[1,,F]


}
