\name{kyotil}
\alias{kyotil}
\title{kyotil}
\description{
Utility functions by Youyi Fong and Krisz Sebestyen, and some functions copied from other packages for convenience (acknowledged on their manual pages).

Most useful functions: mypostscript/mypdf, mytex, %.% (concatenate strings, e.g. "a" %.% "b"), getFormattedSummary. 

See the Index link below for a list of available functions.

The package depends on Hmisc. The main reason for that, besides the usefulness of the package, is Hmisc depends on ggplot2, which also define %.%. By introducing this dependency, ggplot2::%.% is overwritten by kyotil::%.%

}
