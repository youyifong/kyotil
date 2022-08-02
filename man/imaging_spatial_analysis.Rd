\name{get_count_from_xy_coor}
\alias{get_count_from_xy_coor}
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
This function counts cells inside of rectangular box made by the topleft and bottomright xy-coordinates.
}
\value{
The number of masks inside of the rectangular box 
}
\author{
Sunwoo Han
}
\examples{
#get_count_from_xy_coor(file='M926910_Position1_CD3-BUV395_sizes_coordinates.txt', topleft=c(0,0), bottomright=c(1392,600))
}
