\name{misc}
\alias{pava}
\alias{summ}
\alias{empty2na}
\alias{predict.pcc}
\alias{rank.inv.norm}
\alias{INT}
\title{
Misc Functions
}
\description{
Collapse sample strata for COVE boost correlates study
}
\usage{

cove.boost.collapse.strata (dat.ph1, n.demo)

}
\arguments{

  \item{dat.ph1}{
    %%     ~~Describe \code{x} here~~
    }
  \item{n.demo}{
    number of demographics strata, e.g. 6 for COVE correlates
    }
}

\value{
		dat.ph1, whose Wstratum has been updated
}

\details{
  This function is used by both correlates_processing repo and correlates_reporting3 repo
}
