\name{misc}
\alias{pava}
\alias{summ}
\alias{sample.for.cv}
\alias{empty2na}
\alias{predict.pcc}
\alias{rank.inv.norm}
\alias{INT}
\title{
Misc Functions
}
\description{
Misc functions. summ computes iterative sum, sort of like diff.
}
\usage{

pava (x, wt = rep(1, length(x)))
summ(x)
sample.for.cv (dat, v, seed)
empty2na(x) 
\method{predict}{pcc}(object, newdat, ...)
rank.inv.norm(x)
INT(x)

}
\arguments{

  \item{dat}{a data frame. One of the columns must be named y and y should be 0/1 with 1 for case and 0 for control}
  \item{v}{v-fold cross validation}
  \item{seed}{seed for random number generators}
  \item{x}{
    %%     ~~Describe \code{x} here~~
    }
  \item{wt}{
    %%     ~~Describe \code{var.equal} here~~
    }
  \item{object}{
    %%     ~~Describe \code{var.equal} here~~
    }
  \item{newdat}{
    %%     ~~Describe \code{var.equal} here~~
    }
  \item{...}{
    %%     ~~Describe \code{var.equal} here~~
    }
}

\value{
  sample.for.cv returns a list of two vector of integers: train and test, which refer to the rows of dat
}

\details{
 sample.for.cv: case and controls are sampled separately.
 rank.inv.norm: rank-based inverse normal/gaussian transformation
}
