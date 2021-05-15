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
\alias{mydev.off}
\alias{myfigure}
\alias{mylines}
\alias{mytiff}
\alias{myplot.loess}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting Functions}
\description{
mypostscript and mypdf sets the width and height based on mfrow input. 
}
\usage{

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
 NULL)

getMfrow(len)

myhist (x, add.norm=TRUE, col.norm="blue", ...)

myforestplot(dat, xlim = NULL, xlab = "", main = "", col.1 = "red",
 col.2 = "blue", plot.labels = TRUE, order = FALSE,
 decreasing = FALSE, vline = TRUE, cols = NULL, log =
 "", null.val = NULL)

my.interaction.plot(dat, x.ori = 0, xaxislabels = rep("", 2), cex.axis = 1, add = FALSE,
    xlab = "", ylab = "", pcol = NULL, lcol = NULL, ...)

myboxplot(object, ...)

\method{myboxplot}{formula}(formula, data, cex = 0.5, xlab = "", ylab = "", main =
 "", box = TRUE, at = NULL, na.action = NULL, p.val =
 NULL, pch = 1, col = 1, test = "",
 friedman.test.formula = NULL, reshape.formula = NULL,
 reshape.id = NULL, jitter = TRUE, add.interaction =
 FALSE, drop.unused.levels = TRUE, bg.pt = NULL, add =
 FALSE, seed = 1, write.p.at.top = FALSE, ...)
    
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
 legend.title = NULL, legend.cex = 1, legend.inset = 0,
 xaxt = "s", y.intersp = 1.5, x.intersp = 0.3,
 text.width = NULL, add = FALSE, ...)
         
                          
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

mypostscript(file = "temp", mfrow = c(1, 1), mfcol = NULL, width = NULL, 
    height = NULL, ext = c("eps", "pdf", "png", "tiff"), oma = NULL, 
    mar = NULL, main.outer = FALSE, save2file = TRUE, res = 200, 
    ...)

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
 add = FALSE, log = "", same.xylim = FALSE, xlim =
 NULL, ylim = NULL, ...)



}

\arguments{

  \item{cex.cor.dep}{
%%     ~~Describe \code{add} here~~
}
  \item{add.loess.fit}{
%%     ~~Describe \code{add} here~~
}
  \item{leading0}{
%%     ~~Describe \code{add} here~~
}
  \item{null.val}{
%%     ~~Describe \code{add} here~~
}
  \item{write.p.at.top}{
%%     ~~Describe \code{add} here~~
}
  \item{text.width}{
%%     ~~Describe \code{add} here~~
}
  \item{text}{
%%     ~~Describe \code{add} here~~
}
  \item{cex}{
%%     ~~Describe \code{cex} here~~
}
  \item{adj}{
%%     ~~Describe \code{pt2} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{ext}{
%%     ~~Describe \code{ext} here~~
}
  \item{res}{
resolution.
}


  \item{add.norm}{Boolean, whether to add normal approximation density line}
  \item{col.norm}{string, color of added normal density line}
  \item{pt1}{
%%     ~~Describe \code{pt1} here~~
}
  \item{s}{
%%     ~~Describe \code{slope} here~~
}
  \item{ladder}{
%%     ~~Describe \code{slope} here~~
}
  \item{slope}{
%%     ~~Describe \code{slope} here~~
}
  \item{friedman.test.formula}{
%%     ~~Describe \code{slope} here~~
}
  \item{reshape.id}{
%%     ~~Describe \code{slope} here~~
}
  \item{impute.missing.for.line}{
%%     ~~Describe \code{slope} here~~
}
  \item{cor.}{
%%     ~~Describe \code{slope} here~~
}
  \item{mydev}{
%%     ~~Describe \code{slope} here~~
}
  \item{jitter}{Boolean}
  \item{add.interaction}{Boolean}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}

  \item{xaxt}{
%%     ~~Describe \code{pt2} here~~
}

  \item{breaks}{
%%     ~~Describe \code{pt2} here~~
}
  \item{freq}{
%%     ~~Describe \code{pt2} here~~
}
  \item{bg.pt}{
%%     ~~Describe \code{pt2} here~~
}
  \item{probability}{
%%     ~~Describe \code{pt2} here~~
}
  \item{include.lowest}{
%%     ~~Describe \code{pt2} here~~
}
  \item{right}{
%%     ~~Describe \code{pt2} here~~
}
  \item{density}{
%%     ~~Describe \code{pt2} here~~
}
  \item{angle}{
%%     ~~Describe \code{pt2} here~~
}
  \item{border}{
%%     ~~Describe \code{pt2} here~~
}
  \item{axes}{
%%     ~~Describe \code{pt2} here~~
}
  \item{plot}{
%%     ~~Describe \code{pt2} here~~
}
  \item{labels}{
%%     ~~Describe \code{pt2} here~~
}
  \item{nclass}{
%%     ~~Describe \code{pt2} here~~
}
  \item{weight}{
%%     ~~Describe \code{pt2} here~~
}



  \item{pt2}{
%%     ~~Describe \code{pt2} here~~
}
  \item{pt}{
%%     ~~Describe \code{pt2} here~~
}
  \item{alpha}{
%%     ~~Describe \code{pt2} here~~
}

  \item{dat}{
%%     ~~Describe \code{dat} here~~
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
%%     ~~Describe \code{dat2} here~~
}
  \item{add}{
%%     ~~Describe \code{add} here~~
}
  \item{log}{
}
  \item{add.lm.fit}{
}
  \item{add.deming.fit}{
}
  \item{col.lm}{
}
  \item{col.deming}{
}
  \item{reshape.formula}{
a formula object.
}

  \item{xaxislabels}{
%%     ~~Describe \code{xaxislabels} here~~
}
  \item{x.ori}{
%%     ~~Describe \code{x.ori} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{cex.axis}{
%%     ~~Describe \code{cex.axis} here~~
}



  \item{len}{
%%     ~~Describe \code{len} here~~
}

  \item{same.xylim}{Boolean. Whether xlim and ylim should be the same
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{ylim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{col.1}{
%%     ~~Describe \code{col.1} here~~
}
  \item{col.2}{
%%     ~~Describe \code{col.2} here~~
}

  \item{pcol}{
%%     ~~Describe \code{pcol} here~~
}
  \item{lcol}{
%%     ~~Describe \code{lcol} here~~
}

  \item{object}{
%%     ~~Describe \code{object} here~~
}

  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{box}{
%%     ~~Describe \code{box} here~~
}
  \item{at}{
%%     ~~Describe \code{at} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{test}{
string. For example, "t","w","f","k", "tw"
}


  \item{legend}{
%%     ~~Describe \code{legend} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{lty}{
%%     ~~Describe \code{lty} here~~
}
  \item{bty}{
%%     ~~Describe \code{bty} here~~
}

  \item{type}{
%%     ~~Describe \code{type} here~~
}

  \item{make.legend}{
%%     ~~Describe \code{make.legend} here~~
}
  \item{legend.x}{
%%     ~~Describe \code{legend.x} here~~
}
  \item{legend.title}{
%%     ~~Describe \code{legend.title} here~~
}
  \item{legend.cex}{
%%     ~~Describe \code{legend.cex} here~~
}
  \item{draw.x.axis}{
%%     ~~Describe \code{draw.x.axis} here~~
}
  \item{bg}{
%%     ~~Describe \code{bg} here~~
}

  \item{method}{
%%     ~~Describe \code{method} here~~
}


  \item{mfrow}{
%%     ~~Describe \code{mfrow} here~~
}
  \item{mfcol}{
%%     ~~Describe \code{mfcol} here~~
}
  \item{width}{
%%     ~~Describe \code{width} here~~
}
  \item{height}{
%%     ~~Describe \code{height} here~~
}
  \item{oma}{
%%     ~~Describe \code{oma} here~~
}
  \item{mar}{
%%     ~~Describe \code{mar} here~~
}
  \item{main.outer}{
%%     ~~Describe \code{main.outer} here~~
}
  \item{save2file}{
%%     ~~Describe \code{save2file} here~~
}

  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{digits}{
%%     ~~Describe \code{digits} here~~
}
  \item{prefix}{
%%     ~~Describe \code{prefix} here~~
}
  \item{cex.cor}{
}
  \item{plot.labels}{Boolean}
  \item{order}{Boolean}
  \item{decreasing}{Boolean}
  \item{add.diagonal.line}{
%%     ~~Describe \code{add.diagonal.line} here~~
}
  \item{x2}{
%%     ~~Describe \code{add.diagonal.line} here~~
}
  \item{vline}{
%%     ~~Describe \code{add.diagonal.line} here~~
}
  \item{cols}{
%%     ~~Describe \code{add.diagonal.line} here~~
}
  \item{na.action}{
%%     ~~Describe \code{add.diagonal.line} here~~
}
  \item{drop.unused.levels}{
%%     ~~Describe \code{add.diagonal.line} here~~
}


  \item{p.val}{
%%     ~~Describe \code{x} here~~
}
  \item{seed}{
%%     ~~Describe \code{x} here~~
}
  \item{paired}{
%%     ~~Describe \code{x} here~~
}
  \item{show.data.cloud}{
%%     ~~Describe \code{x} here~~
}
  \item{ladder.add.line}{
%%     ~~Describe \code{x} here~~
}
  \item{ladder.add.text}{
%%     ~~Describe \code{x} here~~
}


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

\dontrun{
myfigure(mfrow=c(1,2))
    plot(1:10)
    plot(1:10)
mydev.off(ext="png,pdf", file="tmp")
}



#myboxplot x axis may look weird if log="xy"

}
