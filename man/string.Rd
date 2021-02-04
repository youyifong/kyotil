\name{string.functions}
\alias{string.functions}
\alias{\%.\%}
\alias{contain}
\alias{escapeUnderline}
\alias{fileStem}
\alias{firstIndex}
\alias{getExt}
\alias{getFileStem}
\alias{getStem}
\alias{lastIndex}
\alias{remove.prefix}
\alias{trim}
\title{String Functions}
\description{ \%+\% concatenates its arguments and returns a string.}
\usage{


a \%.\% b

contain(s1, s2)
trim (x, trim.trailing=TRUE, trim.leading=TRUE)  

escapeUnderline(name)

fileStem(file.name)

firstIndex(s1, s2)

getExt(file.name)

getFileStem(file.name)

getStem(file.name)

lastIndex(s1, s2)

remove.prefix(s, sep = "_")

}



\arguments{

  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}

  \item{s1}{
%%     ~~Describe \code{s1} here~~
}
  \item{s2}{
%%     ~~Describe \code{s2} here~~
}


  \item{name}{
%%     ~~Describe \code{name} here~~
}

  \item{file.name}{
%%     ~~Describe \code{file.name} here~~
}






  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{sep}{
%%     ~~Describe \code{sep} here~~
}
  \item{x}{
%%     ~~Describe \code{sep} here~~
}
  \item{trim.leading}{
%%     ~~Describe \code{sep} here~~
}
  \item{trim.trailing}{
%%     ~~Describe \code{sep} here~~
}

}
\examples{

x=1
x \%.\% "b" \%.\% "c"

}
