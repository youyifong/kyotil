\name{plotting}
\alias{abline.pt.slope}
\alias{wtd.hist}
\alias{abline.pts}
\alias{abline.shade}
\alias{abline.shade.2}
\alias{butterfly.plot}
\alias{empty.plot}
\alias{getMfrow}
\alias{whiskers}
\alias{myforestplot}
\alias{my.interaction.plot}
\alias{myboxplot}
\alias{myplot}
\alias{myboxplot.formula}
\alias{myboxplot.data.frame}
\alias{myboxplot.list}
\alias{mylegend}
\alias{mymatplot}
\alias{mypairs}
\alias{myhist}
\alias{mypdf}
\alias{mypng}
\alias{mypostscript}
\alias{panel.cor}
\alias{panel.hist}
\alias{panel.nothing}
\alias{corplot}
\alias{corplot.default}
\alias{corplot.formula}
\alias{add.mtext.label}
\alias{add.mtext.label.2}
\alias{shrink_margin}
\alias{shrink_whitespace}
\alias{reduce_margin}
\alias{reduce_whitespace}
\alias{mydev.off}
\alias{myfigure}
\alias{mylines}
\alias{mytiff}
\alias{myplot.loess}
\alias{smoothed.scaled.hist}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting Functions}
\description{
mypostscript and mypdf sets the width and height based on mfrow input. 
}
\usage{

smoothed.scaled.hist (dat.ls, bin_width, scale.factors=NULL, cols=NULL, 
  legend=NULL, cex.legend=1, ...) 


myplot (object, ...) 

\method{myplot}{loess}(object, xlab="x", ylab="fitted", ...)

whiskers (x, s, ...) 

abline.pt.slope(pt1, slope, x2=NULL, ...)

abline.pts(pt1, pt2 = NULL)

butterfly.plot(dat, dat2 = NULL, add = FALSE, xaxislabels = rep("", 4), x.ori = 0, 
    xlab = "", ylab = "", cex.axis = 1, ...)

empty.plot()

add.mtext.label (text, cex = 1.4, adj = -0.2) 
mydev.off(file = "temp", ext = c("pdf"), res = 200, mydev =
                 NULL, silent = TRUE)

shrink_margin()
shrink_whitespace()
reduce_margin()
reduce_whitespace()
add.mtext.label.2(text, side, line = 1, ...) 

getMfrow(len)

myhist (x, add.norm=TRUE, col.norm="blue", ...)

myforestplot(dat, xlim = NULL, xlab = "", main = "", col.1 = "red",
 col.2 = "blue", plot.labels = TRUE, order = FALSE,
 decreasing = FALSE, vline = TRUE, cols = NULL, log =
 "", null.val = NULL)

my.interaction.plot(dat, x.ori = 0, xaxislabels = rep("", 2), cex.axis = 1, add = FALSE,
    xlab = "", ylab = "", pcol = NULL, lcol = NULL, ...)

myboxplot(object, ...)

\method{myboxplot}{formula}(formula, data, cex = 0.5, xlab = "", ylab = NULL, main
                 = "", box = TRUE, at = NULL, na.action = NULL, p.val =
                 NULL, pch = 1, col = "white", col.points = 1, border =
                 1, test = "", friedman.test.formula = NULL,
                 reshape.formula = NULL, reshape.id = NULL, jitter =
                 TRUE, add.interaction = FALSE, drop.unused.levels =
                 TRUE, bg.pt = NULL, add = FALSE, seed = 1,
                 write.p.at.top = FALSE, ...)
    
\method{myboxplot}{data.frame}(object, cex = 0.5, ylab = "", xlab = "", main = "",
 box = TRUE, at = NULL, pch = 1, col = 1, test = "",
 paired = FALSE, ...)

\method{myboxplot}{list}(object, paired = FALSE, ...)

abline.shade.2(x, col=c(0,1,0))
abline.shade(pt, type = 5, col = c(0, 1, 0), alpha = 0.3)


mylegend(legend, x, y=NULL, lty = NULL, bty = "n", ...)

mymatplot(x, y, type = "b", lty = c(1, 2, 1, 2, 1, 2), pch =
                 NULL, col = rep(c("darkgray", "black"), each = 3),
                 xlab = NULL, ylab = "", draw.x.axis = TRUE, bg = NA,
                 lwd = 1, at = NULL, make.legend = TRUE, legend = NULL,
                 impute.missing.for.line = TRUE, legend.x = 9,
                 legend.title = NULL, legend.cex = 1, legend.lty = lty,
                 legend.inset = 0, xaxt = "s", y.intersp = 1.5,
                 x.intersp = 0.3, text.width = NULL, add = FALSE, ...
)
                          
mypairs(dat, ladder = FALSE, show.data.cloud = TRUE,
 ladder.add.line = T, ladder.add.text = T, ...)

wtd.hist (x, breaks = "Sturges", freq = NULL, probability = !freq, 
    include.lowest = TRUE, right = TRUE, density = NULL, angle = 45, 
    col = NULL, border = NULL, main = paste("Histogram of", xname), 
    xlim = range(breaks), ylim = NULL, xlab = xname, ylab, axes = TRUE, 
    plot = TRUE, labels = FALSE, nclass = NULL, weight = NULL, 
    ...) 

mylines(x, y, type = "l", ...)

myfigure(mfrow = c(1, 1), mfcol = NULL, width = NULL, 
    height = NULL, oma = NULL, mar = NULL, main.outer = FALSE, bg=NULL) 


mypdf(...)

mypng(...)
mytiff(...)

mypostscript(file = "temp", mfrow = c(1, 1), mfcol = NULL, width =
                 NULL, height = NULL, ext = c("eps", "pdf", "png",
                 "tiff"), oma = NULL, mar = NULL, main.outer = FALSE,
                 save2file = TRUE, res = 200, silent = TRUE, ...)

panel.cor(x, y, digits = 2, prefix = "", cex.cor, cor., leading0
 = FALSE, cex.cor.dep = TRUE, ...)

panel.hist(x, ...)

panel.nothing(x, ...)

corplot(object, ...)

\method{corplot}{default}(object, y, ...)

\method{corplot}{formula}(formula, data, main = "", method = c("pearson",
                 "spearman"), col = 1, cex = 0.5, add.diagonal.line =
                 TRUE, add.lm.fit = FALSE, add.loess.fit = FALSE,
                 col.lm = 2, add.deming.fit = FALSE, col.deming = 4,
                 add = FALSE, digit.cor = 2, log = "", same.xylim =
                 FALSE, xlim = NULL, ylim = NULL, ...)



}

\arguments{

\item{digit.cor}{number of digits to print correlation}

\item{col.points}{color of points}

  \item{dat.ls}{
named list of vectors. A histogram is made for each vector.
}
  \item{bin_width}{
width of bin for histograms
}
  \item{scale.factors}{
named vector of scale factors to scale the histogram counts by
}
  \item{cex.legend}{
cex for legend
}


  \item{silent}{
tbdadd
}
  \item{legend.lty}{
tbdadd
}
  \item{cex.cor.dep}{
tbdadd
}
  \item{add.loess.fit}{
tbdadd
}
  \item{leading0}{
tbdadd
}
  \item{null.val}{
tbdadd
}
  \item{write.p.at.top}{
tbdadd
}
  \item{text.width}{
tbdadd
}
  \item{text}{
tbdadd
}
  \item{cex}{
tbdcex
}
  \item{adj}{
tbdpt2
}
  \item{file}{
tbdfile
}
  \item{ext}{
tbdext
}
  \item{res}{
resolution.
}


  \item{add.norm}{Boolean, whether to add normal approximation density line}
  \item{col.norm}{string, color of added normal density line}
  \item{pt1}{
tbdpt1
}
  \item{s}{
tbdslope
}
  \item{ladder}{
tbdslope
}
  \item{slope}{
tbdslope
}
  \item{friedman.test.formula}{
tbdslope
}
  \item{reshape.id}{
tbdslope
}
  \item{impute.missing.for.line}{
tbdslope
}
  \item{cor.}{
tbdslope
}
  \item{mydev}{
tbdslope
}
  \item{jitter}{Boolean}
  \item{add.interaction}{Boolean}
  \item{\dots}{
tbd\dots
}

  \item{xaxt}{
tbdpt2
}

  \item{breaks}{
tbdpt2
}
  \item{freq}{
tbdpt2
}
  \item{bg.pt}{
tbdpt2
}
  \item{probability}{
tbdpt2
}
  \item{include.lowest}{
tbdpt2
}
  \item{right}{
tbdpt2
}
  \item{density}{
tbdpt2
}
  \item{angle}{
tbdpt2
}
  \item{border}{
tbdpt2
}
  \item{axes}{
tbdpt2
}
  \item{plot}{
tbdpt2
}
  \item{labels}{
tbdpt2
}
  \item{nclass}{
tbdpt2
}
  \item{weight}{
tbdpt2
}



  \item{pt2}{
tbdpt2
}
  \item{pt}{
tbdpt2
}
  \item{alpha}{
tbdpt2
}

  \item{dat}{
tbddat
}
  \item{lwd}{
line width.
}
  \item{x.intersp}{
controls the look of legend.
}
  \item{y.intersp}{
controls the look of legend.
}
  \item{legend.inset}{
legend inset
}
  \item{dat2}{
tbddat2
}
  \item{add}{
tbdadd
}
  \item{log}{
  log
}
  \item{add.lm.fit}{
  lm fit
}
  \item{add.deming.fit}{
  add 
}
  \item{col.lm}{
  col
}
  \item{col.deming}{
  col
}
  \item{reshape.formula}{
a formula object.
}

  \item{xaxislabels}{
tbdxaxislabels
}
  \item{x.ori}{
tbdx.ori
}
  \item{xlab}{
tbdxlab
}
  \item{ylab}{
tbdylab
}
  \item{cex.axis}{
tbdcex.axis
}



  \item{len}{
tbdlen
}

  \item{same.xylim}{Boolean. Whether xlim and ylim should be the same
}
  \item{xlim}{
tbdxlim
}
  \item{ylim}{
tbdxlim
}
  \item{main}{
tbdmain
}
  \item{col.1}{
tbdcol.1
}
  \item{col.2}{
tbdcol.2
}

  \item{pcol}{
tbdpcol
}
  \item{lcol}{
tbdlcol
}

  \item{object}{
tbdobject
}

  \item{formula}{
tbdformula
}
  \item{data}{
tbddata
}
  \item{box}{
tbdbox
}
  \item{at}{
tbdat
}
  \item{pch}{
tbdpch
}
  \item{col}{
tbdcol
}
  \item{test}{
string. For example, "t","w","f","k", "tw"
}


  \item{legend}{
tbdlegend
}
  \item{x}{
tbdx
}
  \item{lty}{
tbdlty
}
  \item{bty}{
tbdbty
}

  \item{type}{
tbdtype
}

  \item{make.legend}{
tbdmake.legend
}
  \item{legend.x}{
tbdlegend.x
}
  \item{legend.title}{
tbdlegend.title
}
  \item{legend.cex}{
tbdlegend.cex
}
  \item{draw.x.axis}{
tbddraw.x.axis
}
  \item{bg}{
tbdbg
}

  \item{method}{
tbdmethod
}


  \item{mfrow}{
tbdmfrow
}
  \item{mfcol}{
tbdmfcol
}
  \item{width}{
tbdwidth
}
  \item{height}{
tbdheight
}
  \item{oma}{
tbdoma
}
  \item{mar}{
tbdmar
}
  \item{main.outer}{
tbdmain.outer
}
  \item{save2file}{
tbdsave2file
}

  \item{y}{
tbdy
}
  \item{digits}{
tbddigits
}
  \item{prefix}{
tbdprefix
}
  \item{cex.cor}{
  cex cor
}
  \item{plot.labels}{Boolean}
  \item{order}{Boolean}
  \item{decreasing}{Boolean}
  \item{add.diagonal.line}{
tbdadd.diagonal.line
}
  \item{x2}{
tbdadd.diagonal.line
}
  \item{vline}{
tbdadd.diagonal.line
}
  \item{cols}{
tbdadd.diagonal.line
}
  \item{na.action}{
tbdadd.diagonal.line
}
  \item{drop.unused.levels}{
tbdadd.diagonal.line
}


  \item{p.val}{
tbdx
}
  \item{seed}{
tbdx
}
  \item{paired}{
tbdx
}
  \item{show.data.cloud}{
tbdx
}
  \item{ladder.add.line}{
tbdx
}
  \item{ladder.add.text}{
tbdx
}

\item{side}{
  tbdx
  }
\item{line}{
  tbdx
  }
\item{cex}{
}

\details{
myboxplot shows data points along with boxes. The data poins are jittered and the pattern of jittering is made reproducible in repeated calls. The test can only take one type of test currently.

myforestplot is modified from code from Allan deCamp/SCHARP.
dat should have three columns. first column should be point estimate, second and third lci and uci, fourth p value. 
col.1 is the color used for CIs that do not include null, col.2 is used for CIs that do include null.
If order is TRUE, the rows are ordered by the first column of dat. descreasing can be used to change the behavior of order.

corplot.formula uses MethComp::Deming by Bendix Carstensen to fit Deming regression.

wtd.hist is copied from weights package, author: Josh Pasek.

mymatplot will use na.approx (zoo) to fill in NA before plotting in order to draw continuous lines. The filled-in values will not be shown as points.

smoothed.scaled.hist draws histograms and overlay densities on top.
}

\examples{

set.seed(1)
x=1:50+rnorm(50,0,4)
y=1:50+rnorm(50,0,4)
dat=data.frame(x, y)
corplot(y~x,dat,add.lm.fit=TRUE,add.deming.fit=TRUE,col.lm="red",col.deming="blue")

dat=data.frame(y=c(1:10,2:11), x=rep(c("a","b"),each=10), ptid=c(1:10,1:10))
par(mfrow=c(1,2))
myboxplot(y~x, dat, test="w", jitter=FALSE)
myboxplot(y~x, dat, test="f", add.interaction=TRUE, reshape.formula=y~x, reshape.id="ptid")


myboxplot(list(jitter(1:10), jitter(3:12)), test="w")
myboxplot(list(jitter(1:10), jitter(3:12)), test="w", paired=TRUE)

smoothed.scaled.hist(list(A=rnorm(100,0,1)), bin_width=0.1, xlab="x")
smoothed.scaled.hist(list(A=rnorm(100,0,1), B=rnorm(500,10,2)), 
                  bin_width=0.1, xlab="x")



\dontrun{
myfigure(mfrow=c(1,2))
    plot(1:10)
    plot(1:10)
mydev.off(ext="png,pdf", file="tmp")
}



#myboxplot x axis may look weird if log="xy"

}
