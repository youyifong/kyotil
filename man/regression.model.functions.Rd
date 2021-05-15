\name{regression.model.functions}
\alias{coef.geese}
\alias{coef.tps}
\alias{getFixedEf}
\alias{getFixedEf.coxph}
\alias{getFixedEf.gam}
\alias{getFixedEf.tps}
\alias{getFixedEf.gee}
\alias{getFixedEf.geese}
\alias{getFixedEf.glm}
\alias{getFixedEf.svyglm}
\alias{getFixedEf.svy_vglm}
\alias{getFixedEf.svycoxph}
\alias{getFixedEf.inla}
\alias{getFixedEf.lm}
\alias{getFixedEf.glmerMod}
\alias{getFixedEf.merMod}
\alias{getVarComponent.glmerMod}
\alias{getVarComponent.mer}
\alias{getVarComponent.merMod}
\alias{getFixedEf.lmerMod}
\alias{getFixedEf.lme}
\alias{getFixedEf.logistf}
\alias{getFixedEf.matrix}
\alias{getFixedEf.mer}
\alias{getFixedEf.MIresult}
\alias{getFixedEf2}
\alias{getFormattedSummary}
\alias{getVarComponent}
\alias{getVarComponent.hyperpar.inla}
\alias{getVarComponent.matrix}
\alias{getVarComponent.lmerMod}
\alias{interaction.table}
\alias{predict.geese}
\alias{predict.tps}
\alias{residuals.geese}
\alias{vcov.geese}
\alias{vcov.tps}
\alias{vcov.logistf}
\alias{risk.cal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regression Model Functions
}
\description{
getFormattedSummary prints a table of regression coefficient estimates and standard errors.
}
\usage{

getFormattedSummary(fits, type = 12, est.digits = 2, se.digits = 2,
 robust, random = FALSE, VE = FALSE, to.trim = FALSE,
 rows = NULL, coef.direct = FALSE, trunc.large.est =
 TRUE, scale.factor = 1, p.digits = 3, remove.leading0
 = FALSE, p.adj.method = "fdr", ...)

getVarComponent(object, ...)

getFixedEf(object, ...)

risk.cal(risk, binary.outcome, weights = NULL, ngroups = NULL, 
    cuts = NULL, main = "", add = FALSE, show.emp.risk = TRUE, 
    lcol = 2, ylim = NULL, scale = c("logit", "risk"))
interaction.table(fit, v1, v2, v1.type = "continuous", v2.type = "continuous",
logistic.regression = TRUE)


\method{getFixedEf}{coxph}(object, exp=FALSE,robust=FALSE, ...)

\method{getFixedEf}{gam}(object, ...)

\method{getFixedEf}{gee}(object, exp = FALSE, ...)

\method{getFixedEf}{geese}(object, robust = TRUE, ...)
\method{getFixedEf}{tps}(object, exp=FALSE, robust=TRUE, ...)

\method{getFixedEf}{glm} (object, exp = FALSE, robust = TRUE, ret.robcov = FALSE, 
    ...)  

\method{getFixedEf}{svyglm} (object, exp = FALSE, robust = TRUE, ...)  
\method{getFixedEf}{svy_vglm} (object, exp = FALSE, robust = TRUE, ...)  

\method{getFixedEf}{svycoxph} (object, exp = FALSE, robust = TRUE, ...)  

\method{getFixedEf}{inla}(object, ...)

\method{getFixedEf}{lm}(object, ...)

\method{getFixedEf}{lme}(object, ...)

\method{getFixedEf}{logistf}(object, exp = FALSE, ...)

\method{getFixedEf}{matrix}(object, ...)

\method{getFixedEf}{MIresult}(object, ...)

\method{getVarComponent}{hyperpar.inla}(object, transformation = NULL, ...)

\method{getVarComponent}{matrix}(object, ...)

\method{coef}{geese}(object, ...)
\method{coef}{tps}(object, ...)

\method{predict}{geese}(object, x, ...)
\method{predict}{tps}(object, newdata = NULL, type = c("link", "response"), ...)

\method{residuals}{geese}(object, y, x,...)

\method{vcov}{geese}(object, ...)
\method{vcov}{tps}(object, robust, ...)

\method{vcov}{logistf}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}

  \item{object}{
%%     ~~Describe \code{object} here~~
}

  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{coef.direct}{
%%     ~~Describe \code{fit} here~~
}
  \item{robust}{Boolean, whether to return robust variance estimate
}




  \item{exp}{
%%     ~~Describe \code{exp} here~~
}
  \item{remove.leading0}{
%%     ~~Describe \code{exp} here~~
}
  \item{p.adj.method}{
%%     ~~Describe \code{exp} here~~
}




  \item{cuts}{
%%     ~~Describe \code{fits} here~~
}
  \item{ret.robcov}{
%%     ~~Describe \code{fits} here~~
}

  \item{fits}{
%%     ~~Describe \code{fits} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{est.digits}{
%%     ~~Describe \code{est.digits} here~~
}
  \item{se.digits}{
%%     ~~Describe \code{se.digits} here~~
}
  \item{p.digits}{
%%     ~~Describe \code{se.digits} here~~
}
  \item{random}{
%%     ~~Describe \code{random} here~~
}
  \item{VE}{
%%     ~~Describe \code{random} here~~
}


  \item{transformation}{
%%     ~~Describe \code{transformation} here~~
}


  \item{weights}{
%%     ~~Describe \code{v1} here~~
}
  \item{v1}{
%%     ~~Describe \code{v1} here~~
}
  \item{v2}{
%%     ~~Describe \code{v2} here~~
}
  \item{v1.type}{
%%     ~~Describe \code{v1.type} here~~
}
  \item{v2.type}{
%%     ~~Describe \code{v2.type} here~~
}
  \item{logistic.regression}{
%%     ~~Describe \code{logistic.regression} here~~
}

  \item{newdata}{
%%     ~~Describe \code{x} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}

  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{to.trim}{
%%     ~~Describe \code{y} here~~
}
  \item{rows}{
%%     ~~Describe \code{y} here~~
}
  \item{risk}{
%%     ~~Describe \code{fit} here~~
}
  \item{binary.outcome}{
%%     ~~Describe \code{fit} here~~
}
  \item{ngroups}{
%%     ~~Describe \code{fit} here~~
}
  \item{main}{
%%     ~~Describe \code{fit} here~~
}
  \item{add}{
%%     ~~Describe \code{fit} here~~
}
  \item{show.emp.risk}{
%%     ~~Describe \code{fit} here~~
}
  \item{lcol}{
%%     ~~Describe \code{fit} here~~
}
  \item{ylim}{
%%     ~~Describe \code{fit} here~~
}
  \item{scale}{
%%     ~~Describe \code{fit} here~~
}

  \item{trunc.large.est}{
%%     ~~Describe \code{fit} here~~
}
  \item{scale.factor}{
%%     ~~Describe \code{fit} here~~
}


}
\details{
getFormattedSummary: from a list of fits, say lmer, inla fits, return formatted summary controlled by "type". For a matrix, return Monte Carlo variance
 random=TRUE returns variance components
 type=1: est
 type=2: est (se)
 type=3: est (2.5 percent, 97.5 percent)
 type=4: est   se

getFixedEf returns a matrix, first column coef, second column se, 

getFixedEf.matrix used to get mean and sd from a jags or winbugs sample, getVarComponent.matrix and getFixedEf.matrix do the same thing. Each column of samples is a variable

interaction.table expects coef and vcov to work with fit.
}
\examples{

## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
glm.D9 <- glm(weight ~ group)
getFormattedSummary (list(lm.D9, glm.D9), robust=FALSE)

}
