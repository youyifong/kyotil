\name{roc}
\alias{computeRoc}
\alias{fastauc}
\alias{plotRoc}
\alias{addRoc}
\alias{classification.error}
\title{ROC and AUC}
\description{
ROC/AUC methods.
\code{fastauc} calculates the AUC using a sort operation, instead of summing over pairwise differences in R. \cr
\code{computeRoc} computes an ROC curve. \cr
\code{plotRoc} plots an ROC curve. \cr
\code{addRoc} adds an ROC curve to a plot. \cr
\code{classification.error} computes classification error \cr
}
\usage{
fastauc (score, outcome, t0 = 0, t1 = 1, reverse.sign.if.nece = TRUE, quiet = FALSE)
computeRoc (score, outcome, reverse.sign.if.nece = TRUE, cutpoints
         = NULL)
plotRoc(x, add = FALSE, type = "l", diag.line=TRUE,...)
addRoc (x,...)
classification.error(score, outcome, threshold=NULL, verbose=FALSE)
}
\arguments{
  \item{score}{a vector. Linear combination or score.}
  \item{outcome}{a vector of 0 and 1. Outcome.}
  \item{t0}{a number between 0 and 1 that is the lower boundary of pAUC}
  \item{t1}{a number between 0 and 1 that is the upper boundary of pAUC}
  \item{reverse.sign.if.nece}{a boolean. If TRUE, score is multiplied by -1 if AUC is less than 0.5.}
  \item{x}{a list of two elements: sensitivity and specificity.}
  \item{diag.line}{boolean. If TRUE, a diagonal line is plotted}
  \item{add}{boolean. If TRUE, add to existing plot. If FALSE, create a new plot.}
  \item{quiet}{boolean}
  \item{cutpoints}{cutpoints}
  \item{threshold}{threshold}
  \item{verbose}{boolean}
  \item{type}{line type for \code{lines}}
  \item{...}{arguments passed to \code{plot} or \code{lines}}
}

\value{
  \code{computeRoc} returns a list of sensitivity and specificity. \cr
  \code{plotRoc}  and \code{addRoc} plots ROC curves.\cr
}
\details{
These functions originally come from Thomas Lumley and Tianxi Cai et al.\cr


}

\author{
  Shuxin Yin \email{}\cr
  Youyi Fong \email{youyifong@gmail.com}\cr
  Krisztian Sebestyen \email{}\cr
}

\examples{

n=1e2
score=c(rnorm(n/2,1), rnorm(n/2,0))
outcome=rep(1:0, each=n/2)
# cannot print due to r cmd check
#plotRoc(computeRoc(score, outcome))

# commented out b/c slower on pc and cause note when r cmd check
## test, fastauc2 is a version without all the checking
#score=rnorm(1e5)
#outcome=rbinom(1e5,1,.5)
#system.time(for (i in 1:1e2) fastauc(score,outcome)) # 4.9 sec
#system.time(for (i in 1:1e2) fastauc2(score,outcome)) # 3.8 sec


}
