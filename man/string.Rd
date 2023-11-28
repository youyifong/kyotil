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
a
}
  \item{b}{
b
}

  \item{s1}{
s1
}
  \item{s2}{
s2
}


  \item{name}{
name
}

  \item{file.name}{
file.name
}






  \item{s}{
s
}
  \item{sep}{
sep
}
  \item{x}{
sep
}
  \item{trim.leading}{
sep
}
  \item{trim.trailing}{
sep
}

}
\examples{

x=1
x \%.\% "b" \%.\% "c"

}
