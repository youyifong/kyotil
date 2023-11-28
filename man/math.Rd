\name{math.functions}
\alias{math.functions}
\alias{as.binary}
\alias{binom.coef}
\alias{expit}
\alias{logDiffExp}
\alias{logit}
\alias{logMeanExp}
\alias{logSumExp}
\alias{logSumExpFor2}
\alias{permn}
\alias{Stirling2}
\alias{interpolate}
\title{
Math Functions
}
\description{
H calculates entropy. 
}
\usage{

as.binary(n, base = 2, r = FALSE)

binom.coef(n, m)

expit(x)

logDiffExp(logx1, logx2)

logit(x)

logMeanExp(logx, B = NULL)

logSumExp(logx)

logSumExpFor2(logx, logy)

permn(x, fun = NULL, ...)

Stirling2(n, m)

interpolate(pt1, pt2, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{n}{
tbdn
}
  \item{base}{
tbdbase
}
  \item{r}{
tbdr
}

  \item{m}{
tbdm
}

  \item{pt1}{a vector of length 2}
  \item{pt2}{a vector of length 2}
  \item{x}{
tbdx
}


  \item{logx1}{
tbdlogx1
}
  \item{logx2}{
tbdlogx2
}


  \item{logx}{
tbdlogx
}
  \item{B}{
tbdB
}


  \item{logy}{
tbdlogy
}

  \item{fun}{
tbdfun
}
  \item{\dots}{
tbd\dots
}

}

\examples{
H(rep(1/5,5))
H(rep(3,5))
}
