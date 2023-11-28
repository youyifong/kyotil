\name{crossvalidation}
\alias{sample.for.cv}
\alias{get.kfold.splits}
\alias{get.splits}
\alias{kfold.split}
\alias{lpo.split}
\alias{ran.kfold.split}
\title{
Cross Validation Functions
}
\description{
Cross validation utility functions. 
}
\usage{

sample.for.cv (dat, v, seed)
get.kfold.splits (dat, k, seed)
kfold.split (k, n1, n0)
ran.kfold.split(k, n1, n0, replicates)
lpo.split(n1, n0)
get.splits (dat, cv.scheme=c("LPO","5fold","50xrandom4:1"), seed) 


}
\arguments{

  \item{dat}{a data frame. One of the columns must be named y and y should be 0/1 with 1 for case and 0 for control}
  \item{v}{v-fold cross validation}
  \item{seed}{seed for random number generators}
  \item{k}{
    var.equal
    }
  \item{n1}{
    var.equal
    }
  \item{n0}{
    var.equal
    }
  \item{replicates}{
    var.equal
    }
  \item{cv.scheme}{
    var.equal
    }
}

\value{
  sample.for.cv returns a list of two vector of integers: train and test, which refer to the rows of dat
}

\details{
 sample.for.cv: case and controls are sampled separately.
}
