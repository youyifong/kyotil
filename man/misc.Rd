\name{misc}
\alias{pava}
\alias{summ}
\alias{empty2na}
\alias{predict.pcc}
\alias{rank.inv.norm}
\alias{INT}
\alias{dec_to_binary}
\title{
Misc Functions
}
\description{
Misc functions. summ computes iterative sum, sort of like diff.
}
\usage{

pava (x, wt = rep(1, length(x)))
summ(x)
empty2na(x) 
\method{predict}{pcc}(object, newdat, ...)
rank.inv.norm(x)
INT(x)
dec_to_binary (x,d) 

}
\arguments{

  \item{x}{
    %%     ~~Describe \code{x} here~~
    }
  \item{d}{
    number of digits in the returned binary representation, including leading 0's
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
		summ returns 
}

\details{
 rank.inv.norm: rank-based inverse normal/gaussian transformation

 
 dec_to_binary covert a decimal number to a binary representation with d digits
}
