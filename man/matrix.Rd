\name{matrix.array.functions}
\alias{matrix.array.functions}
\alias{matrix.functions}
\alias{array.functions}
\alias{AR1}
\alias{concatList}
\alias{EXCH}
\alias{fill.jagged.array}
\alias{getMidPoints}
\alias{getUpperRight}
\alias{last}
\alias{mix}
\alias{rep.data.frame}
\alias{rep.matrix}
\alias{rep.matrix.block}
\alias{shift.left}
\alias{shift.right}
\alias{thin.rows}
\alias{ThinRows}
\alias{tr}
\title{
Matrix and Array Functions
}
\description{
concatList returns a string that concatenates the elements of the input list or array

}
\usage{

AR1(p, w)

concatList(lis, sep = "")

EXCH(p, rho)

fill.jagged.array(a)

getMidPoints(x)

getUpperRight(matri, func = NULL)

last(x, n = 1, ...)

mix(a, b)

\method{rep}{data.frame}(x, times = 1, ...)

\method{rep}{matrix}(x, times = 1, each = 1, by.row = TRUE, ...)

\method{rep}{matrix.block}(x, times = 2, ...)

shift.left(x, k = 1)

shift.right(x, k = 1)

thin.rows(dat, thin.factor = 10)

ThinRows(dat, thin.factor = 10)

tr(m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}

  \item{lis}{
list or array
}
  \item{sep}{
%%     ~~Describe \code{sep} here~~
}

  \item{rho}{
%%     ~~Describe \code{rho} here~~
}

  \item{a}{
%%     ~~Describe \code{a} here~~
}

  \item{x}{
%%     ~~Describe \code{x} here~~
}

  \item{matri}{
%%     ~~Describe \code{matri} here~~
}
  \item{func}{
%%     ~~Describe \code{func} here~~
}

  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}

  \item{b}{
%%     ~~Describe \code{b} here~~
}

  \item{times}{
%%     ~~Describe \code{times} here~~
}

  \item{each}{
%%     ~~Describe \code{each} here~~
}
  \item{by.row}{
%%     ~~Describe \code{by.row} here~~
}


  \item{k}{
%%     ~~Describe \code{k} here~~
}


  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{thin.factor}{
%%     ~~Describe \code{thin.factor} here~~
}


  \item{m}{
%%     ~~Describe \code{m} here~~
}
}
\examples{

concatList(1:3,"_")
}
