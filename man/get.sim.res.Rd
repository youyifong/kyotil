\name{get.sim.res}
\alias{get.sim.res}
\alias{MCsummary}
\alias{getFormattedMCSummary}
\title{
Read simulation results
}
\description{
Go through a folder and read all files and combine the results into a multidimensional array.
}
\usage{
get.sim.res (dir, res.name="res", verbose=TRUE)
MCsummary (dir, res.name = "res", exclude.some = TRUE,
                 exclude.col = 1, verbose = TRUE)
getFormattedMCSummary (path, sim, nn, fit.method, exclude.some = TRUE,
                 exclude.col = 1, verbose = TRUE, coef.0 = NULL, digit1
                 = 2, sum.est = c("mean", "median"), sum.sd =
                 c("median", "mean"), style = 1, keep.intercept =
                 FALSE)
}
\arguments{
  \item{dir}{directory of MC result files}
  \item{path}{partial path to the directory of MC result files}
  \item{res.name}{name of the R object saved in the files, default is res, but may be others}
  \item{verbose}{Boolean}
  \item{sim}{a string to denote simulation setting}
  \item{nn}{a vector of sample sizes}
  \item{fit.method}{a string to denote fitting method. sim, nn and fit.method together forms the name of the directory containing MC result files}
  \item{exclude.col}{column number}
  \item{exclude.some}{whether to exclude MC results that are extreme}
  \item{coef.0}{simulation truth}
  \item{digit1}{digits}
  \item{sum.est}{use mean or median as location estimate summary}
  \item{sum.sd}{use mean or median as sd estimate summary}
  \item{style}{integer}
  \item{keep.intercept}{whether to include intercept in the table}
}
\details{
Depends on package abind to combine arrays from files.
}
\value{
A multidimensional array.
}
