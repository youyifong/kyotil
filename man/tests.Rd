\name{testing.functions}
\alias{hosmerlem}
\alias{myfisher.test}
\alias{quick.t.test}
\alias{mycor.test}
\alias{signtest}
\alias{tukey.mtest}
\alias{vector.t.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing Functions
}
\description{
Testing functions.
}
\usage{

hosmerlem(y, yhat, g = 10)

quick.t.test(x, y, var.equal = FALSE)

signtest(x)

tukey.mtest(mu, ms, n)

vector.t.test(mean.x, mean.y, var.x, var.y, n)

myfisher.test(x,y,...)

mycor.test(x, method = c("pearson", "kendall", "spearman"), idx =
 NULL) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{yhat}{
%%     ~~Describe \code{yhat} here~~
}
  \item{g}{
%%     ~~Describe \code{g} here~~
}

  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{var.equal}{
%%     ~~Describe \code{var.equal} here~~
}

  \item{method}{
%%     ~~Describe \code{method} here~~
}


  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{ms}{
%%     ~~Describe \code{ms} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}

  \item{mean.x}{
%%     ~~Describe \code{mean.x} here~~
}
  \item{mean.y}{
%%     ~~Describe \code{mean.y} here~~
}
  \item{var.x}{
%%     ~~Describe \code{var.x} here~~
}
  \item{var.y}{
%%     ~~Describe \code{var.y} here~~
}
  \item{idx}{
%%     ~~Describe \code{var.y} here~~
}
}
\examples{

signtest(runif(10))
}
