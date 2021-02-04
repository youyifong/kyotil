\name{stat.functions}
\alias{stat.functions}
\alias{mutual.info}
\alias{H}
\alias{cor.mixed}
\alias{cor.mixed.default}
\alias{cor.mixed.vector}
\alias{cor.mixed.formula}
\alias{info.cor}
\alias{yule.y}
\alias{kappa.cor}
\alias{l.measure}
\alias{skew}
\title{
Stat Functions
}
\description{
H calculates entropy. 
}
\usage{

H(p, logbase = c("e", "2"))

mutual.info(two.way.table, logbase = c("e", "2"))

cor.mixed(x, ...)

\method{cor.mixed}{default} (x, na.fun, method=c("pearson","spearman"), ...)
\method{cor.mixed}{vector} (x, y, na.fun, method=c("pearson","spearman"), ...) 
\method{cor.mixed}{formula} (formula, data, na.fun, method=c("pearson","spearman"), ...) 

skew (x, na.rm = FALSE) 

info.cor(two.way.table)

yule.y(two.by.two.matrix)

kappa.cor(two.by.two.matrix, weight = c(1, 1), maximum = FALSE)

l.measure(two.by.two.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{p}{
either a count vector or a probability vector, but can not be a vector of membership indicator
}
  \item{logbase}{
%%     ~~Describe \code{logbase} here~~
}
  \item{na.rm}{
%%     ~~Describe \code{logbase} here~~
}

  \item{two.way.table}{
%%     ~~Describe \code{two.way.table} here~~
}

  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}

  \item{na.fun}{
%%     ~~Describe \code{na.fun} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}

  \item{y}{
%%     ~~Describe \code{y} here~~
}

  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}


  \item{two.by.two.matrix}{
%%     ~~Describe \code{two.by.two.matrix} here~~
}

  \item{weight}{
%%     ~~Describe \code{weight} here~~
}
  \item{maximum}{
%%     ~~Describe \code{maximum} here~~
}

}
\examples{
H(rep(1/5,5))
H(rep(3,5))
}
