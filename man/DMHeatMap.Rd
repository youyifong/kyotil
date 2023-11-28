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
tbd
}
  \item{axis}{
tbd
}
  \item{heatmapOnly}{
tbd
}
  \item{verbose}{
tbd
}
  \item{legend.x}{
tbd
}
  \item{legend}{
tbd
}
  \item{Rowv}{
tbdRowv
}
  \item{Colv}{
tbdColv
}
  \item{distfun}{
tbddistfun
}
  \item{hclustfun}{
tbdhclustfun
}
  \item{dendrogram}{
tbddendrogram
}
  \item{symm}{
tbdsymm
}
  \item{scale}{
tbdscale
}
  \item{na.rm}{
tbdna.rm
}
  \item{revC}{
tbdrevC
}
  \item{add.expr}{
tbdadd.expr
}
  \item{breaks}{
tbdbreaks
}
  \item{symbreaks}{
tbdsymbreaks
}
  \item{col}{
tbdcol
}
  \item{colsep}{
tbdcolsep
}
  \item{rowsep}{
tbdrowsep
}
  \item{sepcolor}{
tbdsepcolor
}
  \item{sepwidth}{
tbdsepwidth
}
  \item{cellnote}{
tbdcellnote
}
  \item{notecex}{
tbdnotecex
}
  \item{notecol}{
tbdnotecol
}
  \item{na.color}{
tbdna.color
}
  \item{trace}{
tbdtrace
}
  \item{tracecol}{
tbdtracecol
}
  \item{hline}{
tbdhline
}
  \item{vline}{
tbdvline
}
  \item{linecol}{
tbdlinecol
}
  \item{margins}{
tbdmargins
}
  \item{ColSideColors}{
tbdColSideColors
}
  \item{RowSideColors}{
tbdRowSideColors
}
  \item{cexRow}{
tbdcexRow
}
  \item{cexCol}{
tbdcexCol
}
  \item{labRow}{
tbdlabRow
}
  \item{labCol}{
tbdlabCol
}
  \item{labColor}{
tbdlabColor
}
  \item{key}{
tbdkey
}
  \item{keysize}{
tbdkeysize
}
  \item{density.info}{
tbddensity.info
}
  \item{denscol}{
tbddenscol
}
  \item{symkey}{
tbdsymkey
}
  \item{densadj}{
tbddensadj
}
  \item{main}{
tbdmain
}
  \item{xlab}{
tbdxlab
}
  \item{ylab}{
tbdylab
}
  \item{lmat}{
tbdlmat
}
  \item{lhei}{
tbdlhei
}
  \item{lwid}{
tbdlwid
}
  \item{lower.left.only}{
tbdlower.left.only
}
  \item{\dots}{
tbd\dots
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
