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
%%     ~~Describe \code{n} here~~
}
  \item{base}{
%%     ~~Describe \code{base} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}

  \item{m}{
%%     ~~Describe \code{m} here~~
}

  \item{pt1}{a vector of length 2}
  \item{pt2}{a vector of length 2}
  \item{x}{
%%     ~~Describe \code{x} here~~
}


  \item{logx1}{
%%     ~~Describe \code{logx1} here~~
}
  \item{logx2}{
%%     ~~Describe \code{logx2} here~~
}


  \item{logx}{
%%     ~~Describe \code{logx} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}


  \item{logy}{
%%     ~~Describe \code{logy} here~~
}

  \item{fun}{
%%     ~~Describe \code{fun} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}

}

\examples{
H(rep(1/5,5))
H(rep(3,5))
}
