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
tbd
}
  \item{y}{
tbdy
}
  \item{yhat}{
tbdyhat
}
  \item{g}{
tbdg
}

  \item{x}{
tbdx
}
  \item{var.equal}{
tbdvar.equal
}

  \item{method}{
tbdmethod
}


  \item{mu}{
tbdmu
}
  \item{ms}{
tbdms
}
  \item{n}{
tbdn
}

  \item{mean.x}{
tbdmean.x
}
  \item{mean.y}{
tbdmean.y
}
  \item{var.x}{
tbdvar.x
}
  \item{var.y}{
tbdvar.y
}
  \item{idx}{
tbdvar.y
}
}
\examples{

signtest(runif(10))
}
