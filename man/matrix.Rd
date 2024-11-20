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
\alias{mylast}
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

mylast(x, n = 1, ...)

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
tbdp
}
  \item{w}{
tbdw
}

  \item{lis}{
list or array
}
  \item{sep}{
tbdsep
}

  \item{rho}{
tbdrho
}

  \item{a}{
tbda
}

  \item{x}{
tbdx
}

  \item{matri}{
tbdmatri
}
  \item{func}{
tbdfunc
}

  \item{n}{
tbdn
}
  \item{\dots}{
tbd\dots
}

  \item{b}{
tbdb
}

  \item{times}{
tbdtimes
}

  \item{each}{
tbdeach
}
  \item{by.row}{
tbdby.row
}


  \item{k}{
tbdk
}


  \item{dat}{
tbddat
}
  \item{thin.factor}{
tbdthin.factor
}


  \item{m}{
tbdm
}
}
\examples{

concatList(1:3,"_")
}
