\name{imaging_spatial_analysis}
\alias{imaging_spatial_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imaging analysis for spatial region
}
\description{
Counting the number of masks in a rectangular region
}
\usage{
get_count_from_xy_coor(file, topleft, bottomright)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
_sizes_coordinates.txt
}
  \item{topleft}{
topleft (x,y) coordiate for a rectangular box
}
  \item{bottoright}{
bottomright: bottomright (x,y) coordiate for a rectangular box  
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The number of masks inside of the rectangular box 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{imaging_spatial_analysis}}.
}
\examples{
get_count_from_xy_coor(file='M926910_Position1_CD3-BUV395_sizes_coordinates.txt', topleft=c(0,0), bottomright=c(1392,600))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
