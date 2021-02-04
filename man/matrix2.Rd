\name{matrix2}
\alias{DXD}
\alias{tXDX}
\alias{symprod}
\alias{txSy}
\alias{.as.double}
\title{
Matrix Functions that May Be Faster than %*%
}
\description{
DXD computes D \%*\% X  \%*\% D, where D is a diagonal matrix. tXDX computes t(X) \%*\% D  \%*\% X. symprod computes S \%*\% X for symmetric S. txSy computes t(x) \%*\% S \%*\% y for symmetric S.
}
\usage{

DXD(d1, X, d2)

tXDX(X,D)

symprod(S, X)

txSy(x, S, y)

.as.double(x, stripAttributes = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{d1}{a diagonal matrix or an array}
  \item{d2}{a diagonal matrix or an array}
  \item{x}{array}
  \item{y}{array}
  \item{S}{symmetric matrix}
  \item{X}{matix}
  \item{D}{matix}
  \item{stripAttributes}{boolean}
}
\details{
.as.double does not copying whereas as.double(x) for older versions of R when using .C(DUP = FALSE) make duplicate copy of x. 
In addition, even if x is a 'double', since x has attributes (dim(x)) as.double(x) duplicates

The functions do not check whether S is symmetric. If it is not symmetric, then the result will be wrong.
DXD offers a big gain, while symprod and txSy gains are more incremental.
}
\author{
Krisztian Sebestyen
}
\examples{

d1=1:3
d2=4:6
X=matrix(1:9,3,3)
all(DXD(d1, X, d2) == diag(d1) \%*\% X \%*\% diag(d2))

S=matrix(c(1,2,3,2,4,5,3,5,8),3,3)
X=matrix(1:9,3,3)
all( symprod(S, X) == S \%*\% X )

x=1:3
y=4:6
S=matrix(c(1,2,3,2,4,5,3,5,8),3,3)
txSy(x, S, y) == drop(t(x)\%*\%S\%*\%y)



}
