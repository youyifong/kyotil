\name{DMHeatMap}
\alias{DMHeatMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Better Heatmap Function
}
\description{
Makes a heatmap representation of correaltion coefficients easier.
}
\usage{
DMHeatMap(x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
 distfun = dist, hclustfun = hclust, dendrogram =
 c("both", "row", "column", "none"), symm = FALSE,
 scale = c("none", "row", "column"), na.rm = TRUE, revC
 = identical(Colv, "Rowv"), add.expr, breaks, symbreaks
 = min(x < 0, na.rm = TRUE) || scale != "none", col =
 "heat.colors", colsep, rowsep, sepcolor = "white",
 sepwidth = c(0.05, 0.05), cellnote, notecex = 1,
 notecol = "cyan", na.color = par("bg"), trace =
 c("column", "row", "both", "none"), tracecol = "cyan",
 hline = median(breaks), vline = median(breaks),
 linecol = tracecol, margins = c(5, 5), ColSideColors,
 RowSideColors, cexRow = 0.2 + 1/log10(nr), cexCol =
 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL,
 labColor = NULL, axis = TRUE, heatmapOnly = FALSE, key
 = TRUE, keysize = 1.5, density.info = c("histogram",
 "density", "none"), denscol = tracecol, symkey = min(x
 < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, main
 = NULL, xlab = NULL, ylab = NULL, lmat = NULL, lhei =
 NULL, lwid = NULL, lower.left.only = TRUE, legend =
 TRUE, legend.x = "topright", verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{axis}{
%%     ~~Describe \code{x} here~~
}
  \item{heatmapOnly}{
%%     ~~Describe \code{x} here~~
}
  \item{verbose}{
%%     ~~Describe \code{x} here~~
}
  \item{legend.x}{
%%     ~~Describe \code{x} here~~
}
  \item{legend}{
%%     ~~Describe \code{x} here~~
}
  \item{Rowv}{
%%     ~~Describe \code{Rowv} here~~
}
  \item{Colv}{
%%     ~~Describe \code{Colv} here~~
}
  \item{distfun}{
%%     ~~Describe \code{distfun} here~~
}
  \item{hclustfun}{
%%     ~~Describe \code{hclustfun} here~~
}
  \item{dendrogram}{
%%     ~~Describe \code{dendrogram} here~~
}
  \item{symm}{
%%     ~~Describe \code{symm} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{na.rm}{
%%     ~~Describe \code{na.rm} here~~
}
  \item{revC}{
%%     ~~Describe \code{revC} here~~
}
  \item{add.expr}{
%%     ~~Describe \code{add.expr} here~~
}
  \item{breaks}{
%%     ~~Describe \code{breaks} here~~
}
  \item{symbreaks}{
%%     ~~Describe \code{symbreaks} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{colsep}{
%%     ~~Describe \code{colsep} here~~
}
  \item{rowsep}{
%%     ~~Describe \code{rowsep} here~~
}
  \item{sepcolor}{
%%     ~~Describe \code{sepcolor} here~~
}
  \item{sepwidth}{
%%     ~~Describe \code{sepwidth} here~~
}
  \item{cellnote}{
%%     ~~Describe \code{cellnote} here~~
}
  \item{notecex}{
%%     ~~Describe \code{notecex} here~~
}
  \item{notecol}{
%%     ~~Describe \code{notecol} here~~
}
  \item{na.color}{
%%     ~~Describe \code{na.color} here~~
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
}
  \item{tracecol}{
%%     ~~Describe \code{tracecol} here~~
}
  \item{hline}{
%%     ~~Describe \code{hline} here~~
}
  \item{vline}{
%%     ~~Describe \code{vline} here~~
}
  \item{linecol}{
%%     ~~Describe \code{linecol} here~~
}
  \item{margins}{
%%     ~~Describe \code{margins} here~~
}
  \item{ColSideColors}{
%%     ~~Describe \code{ColSideColors} here~~
}
  \item{RowSideColors}{
%%     ~~Describe \code{RowSideColors} here~~
}
  \item{cexRow}{
%%     ~~Describe \code{cexRow} here~~
}
  \item{cexCol}{
%%     ~~Describe \code{cexCol} here~~
}
  \item{labRow}{
%%     ~~Describe \code{labRow} here~~
}
  \item{labCol}{
%%     ~~Describe \code{labCol} here~~
}
  \item{labColor}{
%%     ~~Describe \code{labColor} here~~
}
  \item{key}{
%%     ~~Describe \code{key} here~~
}
  \item{keysize}{
%%     ~~Describe \code{keysize} here~~
}
  \item{density.info}{
%%     ~~Describe \code{density.info} here~~
}
  \item{denscol}{
%%     ~~Describe \code{denscol} here~~
}
  \item{symkey}{
%%     ~~Describe \code{symkey} here~~
}
  \item{densadj}{
%%     ~~Describe \code{densadj} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{lmat}{
%%     ~~Describe \code{lmat} here~~
}
  \item{lhei}{
%%     ~~Describe \code{lhei} here~~
}
  \item{lwid}{
%%     ~~Describe \code{lwid} here~~
}
  \item{lower.left.only}{
%%     ~~Describe \code{lower.left.only} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\examples{

cor=matrix(runif(15),5,3)
breaks=c(-1,-.7,-.5,-.3,-.1,.1,.3,.5,.7,1)
hU=DMHeatMap(cor, trace="none", symm=FALSE,dendrogram="none", col=RColorBrewer::brewer.pal(
    length(breaks)-1,"RdYlGn"), distfun = function(c) as.dist(1 - c), cexRow =1.5, cexCol =1.5, 
    lmat=rbind( c(2, 1), c(4,3) ), lhei=c(4, 1 ), breaks=breaks, margins=c(2,2), key = FALSE, 
    Rowv=NA, lower.left.only=FALSE)


}
