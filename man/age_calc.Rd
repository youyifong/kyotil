\name{age_calc}
\alias{age_calc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Age Calculation
}
\description{
Calculate age, by Jason P Becker, modified very slightly in how arguments are passed to the function.
}
\usage{
age_calc(dob, enddate = Sys.Date(), units = c("days","months","years"), precise = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dob}{
POSIXlt or Date. Birthday 
}
  \item{enddate}{
POSIXlt or Date. Date to compute age
}
  \item{units}{
string. Choose a unit.
}
  \item{precise}{
Boolean. 
}
}
\references{
http://blog.jsonbecker.com/2013/12/calculating-age-with-precision-in-r.html
}
\author{
Jason P Becker
}
\examples{

age_calc (dob=strptime("29OCT2002", format="\%d\%b\%Y"), 
    enddate=strptime("30OCT2003", format="\%d\%b\%Y"), units='years', precise=TRUE)
age_calc (dob=strptime("29OCT2002", format="\%d\%b\%Y"), 
    enddate=strptime("30DEC2003", format="\%d\%b\%Y"), units='years', precise=FALSE)

}
