\name{getK}
\alias{getK}
\alias{getK}
\title{getK}
\description{
\code{getK} calculates the kernel matrix between X and itself and returns a n by n matrix. Alternatively, it calculates the kernel matrix between X and X2 and returns a n by n2 matrix. 
}
\usage{
getK (X,kernel,para=NULL,X2=NULL,C = NULL)
}
\arguments{
  \item{X}{covariate matrix with dimension n by d. Note this is not the paired difference of covariate matrix.}
  \item{kernel}{string specifying type of kernel:
    polynomial or p (1 + <x,y>)^para, \cr
    rbf or r exp(-para*||x-y||^2), \cr
    linear or l <x,y>, \cr
    ibs or i 0.5*mean(2.0 - |x-y|) or sum(w*(2.0 - |x-y|))/sum(w), with x[i],y[i] in \{0,1,2\} and weights 'w' given in 'para'. \cr
    hamming or h for sum(x == y) with x[i],y[i] binary, \cr
    no default.
  }
  \item{para}{parameter of the kernel fucntion. for \code{ibs} or \code{hamming}, para can be a vector of weights.}
  \item{X2}{optional second covariate matrix with dimension n2 by d}
  \item{C}{logical. If TRUE, kernels are computed by custom routines in C, which may be more memory efficient, and faster too for ibs and hamming kernels.}
}

\details{
    IBS stands for 'Identical By State'. If 'x','y' are in in \{0,1,2\} then \cr
    IBS(x,y) = 0 if |x-y|=2, 1 if |x-y|=1, 2 if |x-y|=0, or IBS(x,y) = 2.0 - |x-y|. \cr
    K(u,v) = sum(IBS(u[i],v[i])) / 2K where K = length(u). \cr
    The 'hamming' kernel is the equivalent of the 'ibs' kernel for binary data. 
    Note that 'hamming' kernel is based on hamming similarity(!), not on dissimilarity distance.
    
    Within in the code, C is default to TRUE for ibs and hamming kernels and FALSE otherwise.
}
\value{
    A kernel matrix.
}

\author{
  Youyi Fong \email{youyifong@gmail.com}\cr
  Krisztian Sebestyen \email{ksebestyen@gmail.com}\cr
  Shuxin Yin \email{}\cr
}

\examples{

X = cbind(x1=rnorm(n=5), x2=rnorm(n=5))
dim(X)
X2 = cbind(x1=rnorm(n=3), x2=rnorm(n=3))
dim(X2)

K = getK(X,"linear")
dim(K)

K = getK(X,"linear",X2=X2)
dim(K)
K1 = getK(X2,"l",X2=X)
dim(K1)
all(K==t(K1))


# RBF kernel
K = getK(X,"rbf",para=1,X2=X2)
K1 = getK(X2,"r",para=1,X2=X)
all(K==t(K1))


# IBS kernel for ternary data 
X <- as.matrix(expand.grid(0:2,0:2))
K = getK(X,kernel = 'ibs')

# add weight
w = runif(ncol(X))
K = getK(X,kernel = 'ibs',para = w) 


# IBS kernel for binary data via option 'h' for 'hamming similarity measure'
X <- as.matrix(expand.grid(0:1,0:1))
K=getK(X,kernel = 'h')

}
