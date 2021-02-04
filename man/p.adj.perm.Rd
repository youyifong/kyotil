\name{p.adj.perm}
\alias{p.adj.perm}
\title{
Permutation-based Multitesting P Values Adjustment
}
\description{
An implementation of Westfall and Young
}
\usage{
p.adj.perm(p.unadj, p.perms, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.unadj}{
%%     ~~Describe \code{p.unadj} here~~
}
  \item{p.perms}{
%%     ~~Describe \code{p.perms} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
}
\details{
This implementation is not as fast as the implementation from the package multtest. But ususally the step to create p.perms is the rate-limiting step.  \cr
\cr

The smallest of the Westfall and Young FWER-controlling multitesting adjusted p values coincides with the p value for testing a global null without any assumptions. 
But for the multitesting adjustment to hold, it requires the subset pivotality condition.


}
\references{
Westfall, P. H., & Young, S. S. (1993). Resampling-based multiple testing: Examples and methods for p-value adjustment (Vol. 279). John Wiley & Sons.
\cr
Westfall, P. H., & Troendle, J. F. (2008). Multiple testing with minimal assumptions. Biometrical Journal: Journal of Mathematical Methods in Biosciences, 50(5), 745-755.
}
\author{
Sue Li, sli@fredhutch.org
}
