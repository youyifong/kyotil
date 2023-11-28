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
\alias{kappacor}
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

kappacor(two.by.two.matrix, weight = c(1, 1), maximum = FALSE)

l.measure(two.by.two.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{p}{
either a count vector or a probability vector, but can not be a vector of membership indicator
}
  \item{logbase}{
tbdlogbase
}
  \item{na.rm}{
tbdlogbase
}

  \item{two.way.table}{
tbdtwo.way.table
}

  \item{x}{
tbdx
}
  \item{\dots}{
tbd\dots
}

  \item{na.fun}{
tbdna.fun
}
  \item{method}{
tbdmethod
}

  \item{y}{
tbdy
}

  \item{formula}{
tbdformula
}
  \item{data}{
tbddata
}


  \item{two.by.two.matrix}{
tbdtwo.by.two.matrix
}

  \item{weight}{
tbdweight
}
  \item{maximum}{
tbdmaximum
}

}
\examples{
H(rep(1/5,5))
H(rep(3,5))
}
