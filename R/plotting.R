myplot <- function(object, ...) UseMethod("myplot") 

# plot x versus fitted
myplot.loess = function(object, xlab="x", ylab="fitted", ...) {
    plot(object$x[order(object$x)], object$fitted[order(object$x)], xlab=xlab, ylab=ylab, ...)    
} 

# one issue with myfigure/mydev.off is that positioning of legend depends on the graphical window size in R
# bg is needed b/c by default the bg is transparent
myfigure=function (mfrow=c(1,1), mfcol=NULL, width=NULL, height=NULL, oma=NULL, mar=NULL, main.outer=FALSE, bg=NULL) {        
    if (!is.null(mfcol)) {
        nrow=mfcol[1]; ncol=mfcol[2]        
    } else {
        nrow=mfrow[1]; ncol=mfrow[2]
    }        
    if(is.null(width) | is.null(height))  tmp=get.width.height(nrow,ncol) else tmp=c(width,height)
 #   unlockBinding(".mydev", getNamespace("kyotil")) #won't pass rcmdcheck 
    eval(eval(substitute(expression(.mydev <<- list(width=tmp[1],height=tmp[2])))))             
    
    if (!is.null(mfcol)) par(mfcol=mfcol) else par(mfrow=mfrow)    
    #needed for dev.copy
    dev.control(displaylist = "enable")
    if (!is.null(oma)) par(oma=oma)
    if (!is.null(mar)) par(mar=mar)    
    if (main.outer) {
        tmp=par()$oma
        tmp[3]=tmp[3]+1
        par(oma=tmp)
    }    
    if(!is.null(bg)) par(bg=bg)
}
mydev.off=function(file="temp", ext=c("pdf"), res=200, mydev=NULL) {        
    if (!is.null(mydev)) .mydev=mydev
    exts=unlist(strsplit(ext, ","))
    tmp=strsplit(file,"/")[[1]]
    for (ext in exts) {
        if (ext=="pdf") {
            subfolder=concatList(c(tmp[-length(tmp)], "pdf"), sep="/")
            filename=if(file.exists(subfolder))  subfolder%.%"/"%.%last(tmp) else file
            dev.copy(pdf,        file=filename%.%"."%.%ext, width=.mydev$width, height=.mydev$height, paper="special")
            cat("Saving figure to "%.%paste(filename,sep="")%.%"."%.%ext%.%"\n")        
        } else if (ext=="eps") {
            subfolder=concatList(c(tmp[-length(tmp)], "eps"), sep="/")
            filename=if(file.exists(subfolder))  subfolder%.%"/"%.%last(tmp) else file
            dev.copy(postscript, file=filename%.%"."%.%ext, width=.mydev$width, height=.mydev$height, paper="special", horizontal=FALSE)
            cat("Saving figure to "%.%paste(filename,sep="")%.%"."%.%ext%.%"\n")        
        } else if (ext=="png") {
            subfolder=concatList(c(tmp[-length(tmp)], "png"), sep="/")
            filename=if(file.exists(subfolder))  subfolder%.%"/"%.%last(tmp) else file
            dev.copy(png,    filename=filename%.%"."%.%ext, width=.mydev$width, height=.mydev$height, units="in", res=res)
            cat("Saving figure to "%.%paste(filename,sep="")%.%"."%.%ext%.%"\n")        
        } else if (ext=="tiff") {
            subfolder=concatList(c(tmp[-length(tmp)], "tiff"), sep="/")
            filename=if(file.exists(subfolder))  subfolder%.%"/"%.%last(tmp) else file
            dev.copy(tiff,   filename=filename%.%"."%.%ext, width=.mydev$width, height=.mydev$height, units="in", res=res, compression="jpeg")
            cat("Saving figure to "%.%paste(filename,sep="")%.%"."%.%ext%.%"\n")        
        }
        dev.off()
    }
    # this resets all pars, important when myfigure contains oma or mar
    resetPar <- function() {
        dev.new()
        op <- par(no.readonly = TRUE)
        dev.off()
        op
    }    
    par(resetPar())
}

get.width.height=function(nrow,ncol){
    if (nrow==1 & ncol==1) {width=6.7; height=6.7
    } else if (nrow==1 & ncol==2) {width=9.7; height=5
    } else if (nrow==1 & ncol==3) {width=9.7; height=3.4
    } else if (nrow==1 & ncol==4) {width=14; height=3.4

    } else if (nrow==2 & ncol==3) {width=9.7; height=6.7
    } else if (nrow==2 & ncol==4) {width=13; height=6.7
    } else if (nrow==2 & ncol==2) {width=8; height=8
    } else if (nrow==2 & ncol==1) {width=6.7; height=9.7
    
    } else if (nrow==3 & ncol==6) {width=17.5; height=9
    } else if (nrow==3 & ncol==7) {width=17.5; height=7
    } else if (nrow==3 & ncol==5) {width=15; height=9.6
    } else if (nrow==3 & ncol==4) {width=12; height=9.6
    } else if (nrow==3 & ncol==3) {width=9.7; height=10.3
    } else if (nrow==3 & ncol==1) {width=6; height=9.7
    } else if (nrow==3 & ncol==2) {width=6.7; height=10.3

    } else if (nrow==4 & ncol==1) {width=4; height=13
    } else if (nrow==4 & ncol==2) {width=6; height=13
    } else if (nrow==4 & ncol==3) {width=9; height=12
    } else if (nrow==4 & ncol==4) {width=9.7; height=10.3
    } else if (nrow==4 & ncol==5) {width=15; height=12.5
    } else if (nrow==4 & ncol==6) {width=15; height=10
    } else if (nrow==4 & ncol==7) {width=17.5; height=9
    } else if (nrow==4 & ncol==9) {width=20; height=9
    } else if (nrow==4 & ncol==8) {width=17.5; height=9
    
    } else if (nrow==5 & ncol==1) {width=5; height=13
    } else if (nrow==5 & ncol==2) {width=7; height=15
    } else if (nrow==5 & ncol==3) {width=9; height=15
    } else if (nrow==5 & ncol==4) {width=12; height=15
    } else if (nrow==5 & ncol==5) {width=15; height=15
    } else if (nrow==5 & ncol==6) {width=9; height=8.3

    } else if (nrow==6 & ncol==5) {width=18; height=17
    } else if (nrow==6 & ncol==3) {width=9; height=19
    } else if (nrow==6 & ncol==4) {width=12; height=19

    } else if (nrow==7 & ncol==3) {width=9; height=22
    } else if (nrow==7 & ncol==5) {width=18; height=19

    } else if (nrow==8 & ncol==5) {width=10; height=16
    } else {
        print ("nrow x ncol not supported: "%.%nrow%.%" x "%.%ncol %.% ". Default to width 10, height 10")
        width=10; height=10
    }
    return(c(width,height))
}

##test

#mypdf(mfrow=c(1,3),file="test1x3");plot(1:10,main="LUMX",xlab="t",ylab="y");plot(1:10);plot(1:10);dev.off()
#mypdf(mfrow=c(2,3),file="test2x3");plot(1:10,main="LUMX",xlab="t",ylab="y");plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);dev.off()
#mypdf(mfrow=c(4,4),file="test4x4");plot(1:10,main="LUMX",xlab="t",ylab="y");plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10,main="Luminex");dev.off()
#mypdf(mfrow=c(1,1),file="test1x1", );plot(1:10,main="LUMX",xlab="t",ylab="y");dev.off()

#    # use convert.ext to convert one format to another
#    # get file name relative to getwd()
#    # .Devices sometimes have "null device", sometimes have "windows"
#    tmp=which(sapply(.Devices, function(x) x=="pdf"))
#    if(length(tmp)==1) {
#        filename=attr(.Devices[[tmp]],"filepath")         
#    }
#    # close pdf device or default device
#    dev.off()
#    # convert file if needed
#    if(png & length(tmp)==1) {
#        system('"C:/Program Files/ImageMagick-7.0.3-Q16/convert.exe" -resize 2000 -density 200 "'%.%getwd()%.%'/'%.%filename%.%'" "'%.%getwd()%.%'/'%.%fileStem(filename)%.%'.png"')
#    }        

# deprecated
# cannot print both to pdf and tiff or print both to screen and pdf, 
mypdf=function (...) mypostscript(ext="pdf",...)
mypng=function(...) mypostscript(ext="png",...)
mytiff=function(...) mypostscript(ext="tiff",...)
mypostscript=function (file="temp", mfrow=c(1,1), mfcol=NULL, width=NULL, height=NULL, ext=c("eps","pdf","png","tiff"), oma=NULL, mar=NULL,main.outer=FALSE, save2file=TRUE, res=200, ...) {    
    
    ext=match.arg(ext)
    
    if (!is.null(mfcol)) {
        nrow=mfcol[1]; ncol=mfcol[2]        
    } else {
        nrow=mfrow[1]; ncol=mfrow[2]
    }
    
    #if (nrow>4) warning ("nrow > 4 will not fit a page without making the figures hard to see")
    
    # sca controls how much to scale down for use in a paper
    if(is.null(width) | is.null(height))  tmp=get.width.height(nrow,ncol) else tmp=c(width,height)
    width=tmp[1]; height=tmp[2]
    
    if(save2file){      
        if (ext=="pdf") {
            pdf (paper="special", file=file%.%"."%.%ext, width=width, height=height, ...)
        } else if (ext=="eps") {
            postscript (paper="special", horizontal=FALSE, file=file%.%"."%.%ext, width=width, height=height, ...)
        } else if (ext=="png") {
            png (filename=file%.%"."%.%ext, width=width, height=height, units="in", res=res, ...)
        } else if (ext=="tiff") {
            tiff (filename=file%.%"."%.%ext, width=width, height=height, units="in", res=res, compression="jpeg", ...)
        }
        cat("Saving figure to "%.%paste(file,sep="")%.%"\n")        
    } else {
        print("not saving to file")
    }
    
    if (!is.null(mfcol)) par(mfcol=mfcol)
    else par(mfrow=mfrow)    
    
    if (!is.null(oma)) par(oma=oma)
    if (!is.null(mar)) par(mar=mar)
    
    if (main.outer) {
        tmp=par()$oma
        tmp[3]=tmp[3]+1
        par(oma=tmp)
    }
    
}



# if lty is specified, a line will be drawn
mylegend=function(legend, x, y=NULL, lty=NULL,bty="n", ...) {
    if(is.null(y)) x=switch(x, "topleft", "top", "topright", "left", "center" , "right", "bottomleft", "bottom", "bottomright")
    legend(bty=bty,x=x, y=y, legend=legend, lty=lty, ...)
}


# copied from pairs help page
# if cex.cor is negative, the sign is reversed and the font of cor is fixed. otherwise, by default, the font of cor is proportional to cor
# allow cor to be spearman or pearson
# will generating lots of warnings, ignore them
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, cor., leading0=FALSE, cex.cor.dep=TRUE, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, method=ifelse(missing(cor.), "spearman", cor.), use="pairwise.complete.obs")
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(!leading0) txt = sub("0","",txt)
    if(missing(cex.cor)) cex.cor <- 2.5
    if(cex.cor.dep) {
        text(0.5, 0.5, txt, cex = cex.cor*ifelse(abs(r)<0.1, sqrt(0.1), sqrt(abs(r)) ))
    } else {
        text(0.5, 0.5, txt, cex = cex.cor) # do this if we don't want cex to depend on correlations
    }
#    print(txt); text(.1, .1, "a"); text(.9, .9, "a")
#    abline(v=1e3)
}
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cornflowerblue", ...)
#    abline(v=1e3)
#    abline(h=.5)
}
panel.smooth.only=function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
    #points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)
}

panel.nothing=function(x, ...) {}
# panel.ladder is a copy of panel.smooth with modifications
panel.ladder=function (x, y, col = par("col"), bg = NA, 
    cex = 1, col.smooth = "red", span = 2/3, iter = 3, cex.text=1.2, cex.star=2,
    cor., add.line=T, add.text=T, ...) 
{
    points(x, y, pch = 20, col = col, bg = bg, cex = cex)
    dxy=as.data.frame(cbind(x,y)); dxy=dxy[complete.cases(dxy),]
    r.lm=lm(y~x, dxy) 
    if(add.line) lines(r.lm, col="1", pred.level=NA)#, args.pband=list(col=SetAlpha("blue",1)) )# 
        
    r <- cor(x, y, method=ifelse(missing(cor.), "spearman", cor.), use="pairwise.complete.obs")
    x.1=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
    y.1=(y-min(y,na.rm=T))/(max(y,na.rm=T)-min(y,na.rm=T))
    p=suppressWarnings(cor.test(x,y, method = ifelse(missing(cor.), "spearman", cor.))$p.value)
    sig=ifelse (round(p,2)<=0.05,ifelse (p<0.01,  ifelse (p<0.01,"***","**")  ,"*"),"")
    if(add.text) text(0*(max(x,na.rm=T)-min(x,na.rm=T))+min(x,na.rm=T), y = 1*(max(y,na.rm=T)-min(y,na.rm=T))+min(y,na.rm=T), labels = "r = "%.%round(r,2), adj = c(0,1), pos = NULL, offset = 0.5, vfont = NULL, cex = cex.text, col = 2, font = 2)
    if(add.text) text(1*(max(x,na.rm=T)-min(x,na.rm=T))+min(x,na.rm=T), y = 1*(max(y,na.rm=T)-min(y,na.rm=T))+min(y,na.rm=T), labels = sig,                 adj = c(1,1), pos = NULL, offset = 0.5, vfont = NULL, cex = cex.star, col = 2, font = 2)

#    ok <- is.finite(x) & is.finite(y)
#    if (any(ok)) 
#        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
#            col = col.smooth, ...)
}
# when log="xy" is passed in, diag and upper panels do not print properly
# cex.labels controls the cex of diagonal panels text size
mypairs=function(dat, ladder=FALSE, show.data.cloud=TRUE, ladder.add.line=T, ladder.add.text=T, ...){
    if(ladder) { # ladder plot
        .pairs(dat, lower.panel=panel.ladder, upper.panel=NULL, diag.panel=NULL, xaxt="n", yaxt="n", gap=0, add.line=ladder.add.line, add.text=ladder.add.text)
    } else {
        .pairs(dat, lower.panel=if (show.data.cloud) panel.smooth else panel.smooth.only, upper.panel=panel.cor, diag.panel=panel.hist, ...)
    }    
}
# a copy of pairs with only one change that is needed to not draw boxes long the diagonal line
.pairs=function (x, labels, panel = points, ...,
          lower.panel = panel, upper.panel = panel,
          diag.panel = NULL, text.panel = textPanel,
          label.pos = 0.5 + has.diag/3, line.main = 3,
          cex.labels = NULL, font.labels = 1, 
          row1attop = TRUE, gap = 1, log = "", show.axis=TRUE){
    if(doText <- missing(text.panel) || is.function(text.panel))
    textPanel <-
        function(x = 0.5, y = 0.5, txt, cex, font)
        text(x, y, txt, cex = cex, font = font)# replacing txt with parse(text=txt) allows us to print math, but it chokes at space
        #text(x, y, parse(text=txt), cex = cex, font = font)# parse(text=txt) allows us to print math
    
    localAxis <- function(side, x, y, xpd, bg, col=NULL, main, oma, ...) {
      ## Explicitly ignore any color argument passed in as
      ## it was most likely meant for the data points and
      ## not for the axis.
        xpd <- NA
        if(side %% 2L == 1L && xl[j]) xpd <- FALSE
        if(side %% 2L == 0L && yl[i]) xpd <- FALSE
        if(side %% 2L == 1L) Axis(x, side = side, xpd = xpd, ...)
        else Axis(y, side = side, xpd = xpd, ...)
    }
    
    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main)
        lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main)
        upper.panel(...)
    
    localDiagPanel <- function(..., main, oma, font.main, cex.main)
        diag.panel(...)
    
    dots <- list(...); nmdots <- names(dots)
    if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for(i in seq_along(names(x))) {
            if(is.factor(x[[i]]) || is.logical(x[[i]]))
               x[[i]] <- as.numeric(x[[i]])
            if(!is.numeric(unclass(x[[i]])))
                stop("non-numeric argument to 'pairs'")
        }
    } else if (!is.numeric(x)) stop("non-numeric argument to 'pairs'")
    panel <- match.fun(panel)
    if((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
        lower.panel <- match.fun(lower.panel)
    if((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
        upper.panel <- match.fun(upper.panel)
    if((has.diag  <- !is.null( diag.panel)) && !missing( diag.panel))
        diag.panel <- match.fun( diag.panel)
    
    if(row1attop) {
        tmp <- lower.panel; lower.panel <- upper.panel; upper.panel <- tmp
        tmp <- has.lower; has.lower <- has.upper; has.upper <- tmp
    }
    
    nc <- ncol(x)
    if (nc < 2) stop("only one column in the argument to 'pairs'")
    
    # names to print on the diagonal
    if(doText) {
        if (missing(labels)) {
            labels <- colnames(x)
            if (is.null(labels)) labels <- paste("var", 1L:nc)
        }
        else if(is.null(labels)) doText <- FALSE
    }
    
    oma <- if("oma" %in% nmdots) dots$oma
    main <- if("main" %in% nmdots) dots$main
    if (is.null(oma))
    oma <- c(4, 4, if(!is.null(main)) 6 else 4, 4)
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    #myprint(oma,gap)
    on.exit(par(opar))
    dev.hold(); on.exit(dev.flush(), add = TRUE)
    
    xl <- yl <- logical(nc)
    if (is.numeric(log)) xl[log] <- yl[log] <- TRUE
    else {xl[] <- grepl("x", log); yl[] <- grepl("y", log)}
    for (i in if(row1attop) 1L:nc else nc:1L)
        for (j in 1L:nc) {
            l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", ""))
            localPlot(x[, j], x[, i], xlab = "", ylab = "",
                      axes = FALSE, type = "n", ..., log = l)
            if( i==j | (i < j && has.lower) || (i > j && has.upper) ) { 
                if (i!=j) box()# YF: add i==j 
                if(show.axis){
                    if(i == 1  && (!(j %% 2L) || !has.upper || !has.lower ))
                        localAxis(1L + 2L*row1attop, x[, j], x[, i], ...)
                    if(i == nc && (  j %% 2L  || !has.upper || !has.lower ))
                        localAxis(3L - 2L*row1attop, x[, j], x[, i], ...)
                    if(j == 1  && (!(i %% 2L) || !has.upper || !has.lower ))
                        localAxis(2L, x[, j], x[, i], ...)
                    if(j == nc && (  i %% 2L  || !has.upper || !has.lower ))
                        localAxis(4L, x[, j], x[, i], ...)
                }
                mfg <- par("mfg")
                if(i == j) {
                    if (has.diag) localDiagPanel(as.vector(x[, i]), ...)
            if (doText) {
                        par(usr = c(0, 1, 0, 1))
                        if(is.null(cex.labels)) {
                            l.wid <- strwidth(labels, "user")
                            cex.labels <- max(0.8, min(2, .9 / max(l.wid)))
                        }
                        xlp <- if(xl[i]) 10^0.5 else 0.5
                        ylp <- if(yl[j]) 10^label.pos else label.pos
                        text.panel(xlp, ylp, labels[i],
                                   cex = cex.labels, font = font.labels)
                    }
                } else if(i < j)
                    localLowerPanel(as.vector(x[, j]), as.vector(x[, i]), ...)
                else
                    localUpperPanel(as.vector(x[, j]), as.vector(x[, i]), ...)
                if (any(par("mfg") != mfg))
                    stop("the 'panel' function made a new plot")
            } else par(new = FALSE)
    
        }
    if (!is.null(main)) {
        font.main <- if("font.main" %in% nmdots) dots$font.main else par("font.main")
        cex.main <- if("cex.main" %in% nmdots) dots$cex.main else par("cex.main")
        mtext(main, 3, line.main, outer=TRUE, at = 0.5, cex = cex.main, font = font.main)
    }
    invisible(NULL)
}


getMfrow=function (len) {
    ret=NULL
    if (len==1) { 
        ret=c(1,1)
    } else if (len==2) { 
        ret=c(1,2)
    } else if (len==3) { 
        ret=c(1,3)
    } else if (len<=4) { 
        ret=c(2,2)
    } else if (len<=6) { 
        ret=c(2,3)
    } else if (len<=9) { 
        ret=c(3,3)
    } else if (len<=12) { 
        ret=c(3,4)
    } else if (len<=16) { 
        ret=c(4,4)
    } else if (len<=20) { 
        ret=c(4,5)
    } else if (len<=25) { 
        ret=c(5,5)
    }
    ret
}


empty.plot=function () {
    plot(1,1,type="n",xlab="",ylab="",xaxt="n", yaxt="n", bty="n")
}

# dat
myforestplot=function(dat, xlim=NULL, xlab="", main="", col.1="red", col.2="blue",plot.labels=TRUE,order=FALSE,decreasing=FALSE, vline=TRUE,cols=NULL,log="",null.val=NULL) {
    if (order) dat=dat[order(dat[,1],decreasing=decreasing),] 
    p=nrow(dat)    
    # makes no plot, but helps set the x axis later
    plot(c(dat[,2], dat[,3]),rep(1:p,2), xlim=xlim, yaxt="n", xaxt="s", xlab=xlab, ylab="", main="", type="n", cex.main=1.4, axes=F, log=log)
    mtext(side=3, line=3, adj=0, text=main, cex=1.4, font=2, xpd=NA)
    if(is.null(null.val)) {
        if (range(dat[,2:3])[1]>0) null.val=1 else null.val=0 # if all values are greater than 0, 1 is probably the null, otherwise, 0 is probably the null
    }
    if(vline) abline(v=null.val, col="gray")
    if(is.null(cols)) cols=ifelse(dat[,4]<0.05, col.1, col.2)
    points(dat[,1], nrow(dat):1, pch=19, col=cols)
    segments(dat[,2], nrow(dat):1, dat[,3], nrow(dat):1, lwd=2, col=cols)
    axis(1, cex.axis=1.4)
    # add labels
    if (plot.labels) axis(4, at=p:1, rownames(dat), tick=F, las=2, col=1, cex.axis=1, xpd=NA, line=-.5)
}


myboxplot <- function(object, ...) UseMethod("myboxplot") 


# this function may fail sometimes, most likely at eval
# myboxplot.formula and myboxplot.list make a boxplot with data points and do inferences for two group comparions. 
# cex=.5; ylab=""; xlab=""; main=""; box=FALSE; highlight.list=NULL; at=NULL;pch=1;col=1;
# friedman.test.formula is of the form a ~ b | c
myboxplot.formula=function(formula, data, cex=.5, xlab="", ylab="", main="", box=TRUE, at=NULL, na.action=NULL, p.val=NULL,
    pch=1, col=1, test="", friedman.test.formula=NULL, reshape.formula=NULL, reshape.id=NULL, jitter=TRUE, add.interaction=FALSE,  drop.unused.levels = TRUE, bg.pt=NULL, add=FALSE, seed=1, write.p.at.top=FALSE, ...){
    
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (inherits(save.seed,"try-error")) {        
        set.seed(1)
        save.seed <- get(".Random.seed", .GlobalEnv)
    }                        
    set.seed(seed)
    
     # removes empty groups formed through model.frame
     mf=model.frame(formula, data)
     response <- attr(attr(mf, "terms"), "response") 
     tmp.dat=split(mf[[response]], mf[-response])
     if(drop.unused.levels) {
         len.n=sapply(tmp.dat, length)
         tmp.dat=tmp.dat[len.n!=0]
     }
     res=boxplot(tmp.dat, range=0, xlab=xlab, at=at, col=NULL, cex=cex, 
        boxlty=if(!box) 0 else NULL,whisklty=if(!box) 0 else NULL,staplelty=if(!box) 0 else NULL,
        #pars = list(boxwex = if(box) 0.8 else 0, staplewex = if(box) 0.5 else 0, outwex = if(box) 0.5 else 0), 
        main=main, ylab=ylab, add=add, ...)
    
    # na.action is key below b/c otherwise pch vector will be out of sync with data when there are missing data
    dat.tmp=model.frame(formula, data, na.action=NULL);# str(dat.tmp); str(data)
    xx=interaction(dat.tmp[,-1]); #str(xx); print(table(xx))
    if(drop.unused.levels) xx=droplevels(xx)
    if(is.null(at)){        
        xx=as.numeric(xx); #print(table(xx))
    } else{
        xx=at[xx]
    }      
    if (add.interaction) jitter=FALSE
    if (jitter) xx=jitter(xx)  
    points(xx, dat.tmp[[1]], cex=cex,pch=pch,col=col, bg=bg.pt)
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)     
    
    # inference
    # if p.val is passed in, then use that p.val
    x.unique=unique(dat.tmp[[2]])
    if (length(test)>0) {
        sub=""
        pvals=NULL
        if ("t" %in% test) {
            if (is.null(p.val)) p.val=t.test(formula, data)$p.value
            pvals=c(pvals, Student=p.val)
            if (write.p.at.top) sub=paste0("p=",signif(p.val,2)) else sub=sub%.%" Student's t "%.%ifelse(length(test)==1,"p-val ","")%.%signif(p.val,2) 
        }
        if ("w" %in% test) {
            if (is.null(p.val)) p.val=suppressWarnings(wilcox.test(formula, data)$p.value)
            pvals=c(pvals, Wilcoxon=p.val)
            if (write.p.at.top) sub=paste0("p=",signif(p.val,2)) else sub=sub%.%" Wilcoxon "%.%ifelse(length(test)==1,"p-val ","")%.%signif(p.val,2)
        }
        if ("k" %in% test) {
            if (is.null(p.val)) p.val=kruskal.test(formula, data)$p.value
            pvals=c(pvals, Kruskal=p.val)
            if (write.p.at.top) sub=paste0("p=",signif(p.val,2)) else sub=sub%.%" Kruskal "%.%ifelse(length(test)==1,"p-val ","")%.%signif(p.val,2)
        }
        if ("f" %in% test) {
            if (!is.null(friedman.test.formula)) {
            # if there is missing data, this won't work, try the else and supply reshape.formula and reshape.id
                if (is.null(p.val)) p.val=friedman.test(friedman.test.formula, data)$p.value
                pvals=c(pvals, Friedman=p.val)
                if (write.p.at.top) sub=paste0("p=",signif(p.val,2)) else sub=sub%.%" Friedman "%.%ifelse(length(test)==1,"p-val ","")%.%signif(p.val,2)
            } else if (!is.null(reshape.formula) & !is.null(reshape.id)) {
                dat.wide=myreshapewide (reshape.formula, data, idvar = reshape.id)
                #str(dat.wide)# show this so that we know we are using the right data to do the test
                ftest = try(friedman.test (as.matrix(dat.wide[,-(1)])),silent=T)
                if (is.null(p.val)) if (!inherits(ftest,"try-error")) p.val=ftest$p.value else p.val=NA
                pvals=c(pvals, Friedman=p.val)
                if (add.interaction) my.interaction.plot(as.matrix(dat.wide[,-1]), add=T)
                if (write.p.at.top) sub=paste0("p=",signif(p.val,2)) else sub=sub%.%" Friedman "%.%ifelse(length(test)==1,"p-val ","")%.% ifelse (is.na(p.val), "NA", signif(p.val,2))
            } else warning("cannot perform Friedman test without friedman.test.formula or reshape.formula,reshape.id")
        }
        if (write.p.at.top) {
            if (pvals[1]<0.05) title(main=sub, line=.5, font.main=3)
        } else title(sub=sub, line=2.5)
        res$pvals=pvals
    }
    
    invisible(res)
    
}

myboxplot.data.frame=function(object, cex=.5, ylab="", xlab="", main="", box=TRUE, at=NULL, pch=1, col=1, test="", paired=FALSE, ...){
    myboxplot.list(as.list(object), cex=cex, ylab=ylab, xlab=xlab, main=main, box=box, at=at, pch=pch, col=col, test=test, ...)
}
myboxplot.matrix=function(object, cex=.5, ylab="", xlab="", main="", box=TRUE, at=NULL, pch=1, col=1, test="", paired=FALSE, ...){
    myboxplot.list(as.list(as.data.frame(object)), cex=cex, ylab=ylab, xlab=xlab, main=main, box=box, at=at, pch=pch, col=col, test=test, ...)
}

myboxplot.list=function(object, paired=FALSE, ...){
    
    # make a dataframe out of list object
    dat=NULL
    if(is.null(names(object))) names(object)=1:length(object)
    for (i in 1:length(object)) {
        dat=rbind(dat,data.frame(y=object[[i]], x=names(object)[i]))
    }
    
    p.val=NULL
    if (paired) {
        if(length(object)==2) {
            p.val = wilcox.test(object[[1]], object[[2]], paired=TRUE)$p.value
        } else {
            stop("when it is paired, we expect a list of two")
        }
    }
        
    dat$x=factor(dat$x, levels=names(object))
    myboxplot(y~x, dat, p.val=p.val, ...)
    
}


# can use after myboxplot
# both dat must have two columns, each row is dat from one subject
# x.ori=0; xaxislabels=rep("",2); cex.axis=1; add=FALSE; xlab=""; ylab=""; pcol=NULL; lcol=NULL
my.interaction.plot=function(dat, x.ori=0, xaxislabels=rep("",2), cex.axis=1, add=FALSE, xlab="", ylab="", pcol=NULL, lcol=NULL, ...){
    if (!add) plot(0,0,type="n",xlim=c(1,2),ylim=range(dat), ylab=ylab, xlab=xlab, xaxt="n", ...)
    cex=.25; pch=19
    if (is.null(lcol)) lcol=ifelse(dat[,1]>dat[,2],"red","black") else if (length(lcol)==1) lcol=rep(lcol,nrow(dat))
    if (!is.null(pcol)) if (length(pcol)==1) pcol=matrix(pcol,nrow(dat),2)
    for (i in 1:nrow(dat)) {
        points (1+x.ori, dat[i,1], cex=cex, pch=pch, col=ifelse(is.null(pcol), 1, pcol[i,1]))
        points (2+x.ori, dat[i,2], cex=cex, pch=pch, col=ifelse(is.null(pcol), 1, pcol[i,2]))
        lines (1:2+x.ori, dat[i,], lwd=.25, col=lcol[i])
    }
    axis(side=1, at=1:2+x.ori, labels=xaxislabels, cex.axis=cex.axis)
}

# called butterfly.plot, because it is meant to plot two treatment arms at two time points, the two arms are plotted in a mirror fashion, see "by analyte.pdf" for an example
# if dat2 is null: dat is matrix with four columns. each row is one subject, the columns will be plotted side by side, with lines connecting values from one ptid
# if dat2 is not null, dat has two columns, which are plotted side by side with lines connecting them, same for dat2
# if add is true, no plot function will be called
butterfly.plot=function (dat, dat2=NULL, add=FALSE, xaxislabels=rep("",4), x.ori=0, xlab="", ylab="", cex.axis=1, ...){
    if (!add) plot(0,0,type="n",xlim=c(1,4),ylim=range(dat), xaxt="n", xlab=xlab, ylab=ylab, ...)
    for (i in 1:nrow(dat)) {
        lines (1:2+x.ori, dat[i,1:2], lwd=.25, col=ifelse(dat[i,1]<=dat[i,2],"red","black"))
        if (is.null(dat2)) {
            lines (2:3+x.ori, dat[i,2:3], lwd=.25, col="lightgreen")
            lines (3:4+x.ori, dat[i,3:4], lwd=.25, col=ifelse(dat[i,3]<=dat[i,4],"black","red"))
        }
    }
    if (!is.null(dat2)) {
        for (i in 1:nrow(dat2)) {
            lines (3:4+x.ori, dat2[i,1:2], lwd=.25, col=ifelse(dat2[i,1]<=dat2[i,2],"black","red"))
        }
    }
    axis(side=1, at=1:4+x.ori, labels=xaxislabels, cex.axis=cex.axis)
}


corplot <- function(object, ...) UseMethod("corplot") 

corplot.default=function(object,y,...){
    dat=data.frame(object,y)
    names(dat)=c("x1", "x2")
    corplot(x2~x1, dat, ...)
}

# col can be used to highlight some points
corplot.formula=function(formula,data,main="",method=c("pearson","spearman"),col=1,cex=.5,add.diagonal.line=TRUE,add.lm.fit=FALSE,add.loess.fit=FALSE,col.lm=2,add.deming.fit=FALSE,col.deming=4,add=FALSE,
    log="",same.xylim=FALSE,xlim=NULL,ylim=NULL, ...){
    vars=dimnames(attr(terms(formula),"factors"))[[1]]
    cor.=NULL
    if (length(method)>0) {
        cor.=sapply (method, function (method) {
            cor(data[,vars[1]],data[,vars[2]],method=method,use="p")
        })
        tmp=main==""
        main=main%.%ifelse(tmp, "", " (")
        main=main%.%"cor: "%.%concatList(round(cor.,2),"/")
        main=main%.%ifelse(tmp, "", ")")
    }

    if (!add) {
        if (same.xylim) {
            xlim=range(model.frame(formula, data))
            ylim=range(model.frame(formula, data))
        }
        plot(formula,data,main=main,col=col,cex=cex,log=log,xlim=xlim,ylim=ylim,...)
    } else {
        points(formula,data,main=main,col=col,cex=cex,log=log,...)
    }
    
    if(add.diagonal.line) abline(0,1)
    if(add.lm.fit) {
        fit=lm(formula, data)
        abline(fit,untf=log=="xy", col=col.lm)
    }
    if(add.loess.fit) {
        fit=loess(formula, data)
        mylines(fit$x,fit$fitted,col=col.lm)
    }
    if(add.deming.fit) {
        # this implementation is faster than the one by Therneau, Terry M.
        fit=Deming(model.frame(formula, data)[[2]], model.frame(formula, data)[[1]]) # this function is in Deming.R copied from MethComp package by Bendix Carstensen
        abline(coef(fit)["Intercept"], coef(fit)["Slope"], untf=log=="xy", col=col.deming)   
        # Therneau, Terry M.'s implementation in a loose R file that is in 3software folder, slower than Deming, but may be more generalized?
        #fit <- deming(model.frame(formula, data)[[2]], model.frame(formula, data)[[1]], xstd=c(1,0), ystd=c(1,0))
        #abline(fit, untf=log=="xy", col=col.deming)        
    }
    
    invisible(cor.)
}

abline.pts=function(pt1, pt2=NULL){
    if (is.null(pt2)) {
        if (nrow(pt1)>=2) {
            pt2=pt1[2,]
            pt1=pt1[1,]
        } else {
            stop("wrong input")
        }
    }
    slope=(pt2-pt1)[2]/(pt2-pt1)[1]
    intercept = pt1[2]-slope*pt1[1]
    abline(intercept, slope)
}
#abline.pts(c(1,1), c(2,2))

abline.pt.slope=function(pt1, slope, x2=NULL, ...){
    intercept = pt1[2]-slope*pt1[1]
    if (is.null(x2)) {
        abline(intercept, slope, ...)
    } else {
        pt2=c(x2, intercept+slope*x2)
        lines(c(pt1[1],pt2[1]), c(pt1[2],pt2[2]), ...)
    }
    
}

# put a shade in a rectangle between a point and one of the four quadrants
# pt is a vector of two values
# col is red blue gree
abline.shade=function(pt, type=5, col=c(0,1,0), alpha=0.3){
    usr <- par('usr')   
    # rec: xleft, ybottom, xright, ytop
    # usr: xleft, xright, ybottom, ytop
    # the first four types are quadrant
    if (type==1) {
        rect(pt[1], pt[2], usr[2], usr[4], col=rgb(red=col[1], blue=col[2], green=col[3], alpha=alpha), border=NA) 
    } else if (type==2) {
        rect(pt[1], usr[3], usr[2], pt[2], col=rgb(red=col[1], blue=col[2], green=col[3], alpha=alpha), border=NA) 
    } else if (type==3) {
        rect(usr[1], usr[3], pt[1], pt[2], col=rgb(red=col[1], blue=col[2], green=col[3], alpha=alpha), border=NA) 
    } else if (type==4) {
        rect(usr[1], pt[2], pt[1], usr[4], col=rgb(red=col[1], blue=col[2], green=col[3], alpha=alpha), border=NA) 
    } else if (type==5) {
        # between pt[1] and pt[2] horizontally
        rect(pt[1], usr[3], pt[2], usr[4], col=rgb(red=col[1], blue=col[2], green=col[3], alpha=alpha), border=NA) 
    } else stop("type not recognized")
    
}

# put a shade between two lines
# x is a vector of two values
# col is red blue gree
abline.shade.2=function(x, col=c(0,1,0)){
    usr <- par('usr') 
    rect(x[1], usr[3], x[2], usr[4], col=rgb(red=col[1], blue=col[2], green=col[3], alpha=.5), border=NA) 
}


#abline.pt.slope(c(1,1), 1)
# When impute.missing.for.line is TRUE, lines are drawn even when there are missing values in between two observations
mymatplot=function(x, y, type="b", lty=c(1,2,1,2,1,2), pch=NULL, col=rep(c("darkgray","black"),each=3), xlab=NULL, ylab="", 
    draw.x.axis=TRUE, bg=NA, lwd=1, at=NULL, make.legend=TRUE, legend=NULL, impute.missing.for.line=TRUE,
    legend.x=9, legend.title=NULL, legend.cex=1, legend.lty=lty, legend.inset=0, xaxt="s", y.intersp=1.5, x.intersp=0.3, text.width=NULL, 
    add=FALSE, ...) {
    
    missing.y=FALSE
    if (missing(y)) {
        missing.y=TRUE
        y=x
        x=1:nrow(y)
    } 
    
    # fill in missing values if necessary, i.e. when there are lines to draw
    if (impute.missing.for.line & any(is.na(y)) & type %in% c("l","b")) {
        y.imputed=zoo::na.approx(y, x=x, na.rm=FALSE); rownames(y.imputed)=rownames(y)
        imputed=TRUE
        cat("imputing data ...\n")
    } else {
        imputed=FALSE
        y.imputed=y
    }
    
    if (is.null(xlab)) xlab=names(dimnames(y))[1]
    if (is.null(legend.title)) legend.title=names(dimnames(y))[2]
    # draw line first, then points
    if(type %in% c("l","b")) matplot(x, y.imputed, lty=lty, pch=pch, col=col, xlab=xlab, xaxt=xaxt, ylab=ylab, bg=bg, lwd=lwd, type="l", add=add, ...)
    if(type %in% c("p","b")) matplot(x, y, lty=lty, pch=pch, col=col, xlab=xlab, xaxt=xaxt, ylab=ylab, bg=bg, lwd=lwd, type="p", add=add | type!="p", ...)
    if (xaxt=="n") 
        if(missing.y & draw.x.axis) axis(side=1, at=if(is.null(at)) x else at, labels=if(is.null(at)) rownames(y) else at) else if (draw.x.axis) {
            axis(side=1, at=if(is.null(at)) x else at, labels=if(is.null(at)) x else at)
        }
    if (make.legend) {
        if (is.null(legend)) legend=colnames(y)
        if (length(unique(pch))>1) {
            mylegend(legend, x=legend.x, lty=legend.lty, title=legend.title, col=col, pt.bg=bg, cex=legend.cex, lwd=lwd, inset=legend.inset, y.intersp=y.intersp, x.intersp=x.intersp, text.width=text.width, pch=pch)
        } else {
            mylegend(legend, x=legend.x, lty=legend.lty, title=legend.title, col=col, pt.bg=bg, cex=legend.cex, lwd=lwd, inset=legend.inset, y.intersp=y.intersp, x.intersp=x.intersp, text.width=text.width)
        }
    }
}


myhist=function(x, add.norm=TRUE, col.norm="blue", ...){
    if (!add.norm) hist(x, ...) else {
        hist=hist(x,breaks=30, plot=F)
        dnorm=dnorm(seq(range(x)[1],range(x)[2], length=100),mean(x),sd(x))
        hist(x, freq=F, ylim=range(hist$density, dnorm), ...)
        lines(seq(range(x)[1],range(x)[2], length=100), dnorm, col=col.norm)
    }    
}


# eclipse
plot.ellipse=function(x0,y0,a=1,b=1,theta=0,alpha=0,add=TRUE,...) {
    theta <- seq(0, 2 * pi, length=500)
#    x <- x0 + a * cos(theta)
#    y <- y0 + b * sin(theta)    
    x <- x0 + a * cos(theta) * cos(alpha) - b * sin(theta) * sin(alpha)
    y <- y0 + a * cos(theta) * sin(alpha) + b * sin(theta) * cos(alpha)
    if(add) {
        lines(x,y,...)
    } else plot(x, y, type = "l",...)
}

add.mtext.label=function(text, cex=1.4, adj=-0.2) mtext(side=3, line=2, adj=adj, text=text, cex=cex, font=2, xpd=NA)


# copied from DescTools, Andri Signorell
# just for check not to bark!
utils::globalVariables(c("hred","hblue"))

lines.lm <- function (x, col = Pal()[1], lwd = 2, lty = "solid",
                      type = "l", n = 100, conf.level = 0.95, args.cband = NULL,
                      pred.level = NA, args.pband = NULL, ...) {
    
  mod <- x$model

  # we take simply the second column of the model data.frame to identify the x variable
  # this will crash, if there are several resps and yield nonsense if there is
  # more than one pred,
  # so check for a simple regression model y ~ x (just one resp, just one pred)

  # Note:
  # The following will not work, because predict does not correctly recognise the newdata data.frame:
  # lines(lm(d.pizza$temperature ~ d.pizza$delivery_min), col=hred, lwd=3)
  # see what happens to the data.frame colnames in: predict(x, newdata=data.frame("d.pizza$delivery_min"=1:20))
  # this predict won't work.
  # always provide data:    y ~ x, data

  # thiss is not a really new problem:
  # http://faustusnotes.wordpress.com/2012/02/16/problems-with-out-of-sample-prediction-using-r/

  # we would only plot lines if there's only one predictor

  pred <- all.vars(formula(x)[[3]])
  if(length(pred) > 1) {
    stop("Can't plot a linear model with more than 1 predictor.")
  }

  # the values of the predictor
  xpred <- x$model[, pred] # modified by YF
  #xpred <- eval(x$call$data)[, pred]

  newx <- data.frame(seq(from = min(xpred, na.rm = TRUE),
                         to = max(xpred, na.rm = TRUE), length = n))

  colnames(newx) <- pred
  fit <- predict(x, newdata = newx)

  if (!(is.na(pred.level) || identical(args.pband, NA)) ) {
    args.pband1 <- list(col = SetAlpha(col, 0.12), border = NA)
    if (!is.null(args.pband))
      args.pband1[names(args.pband)] <- args.pband

    ci <- predict(x, interval="prediction", newdata=newx, level=pred.level) # Vorhersageband
    do.call("DrawBand", c(args.pband1, list(x = c(unlist(newx), rev(unlist(newx)))),
                          list(y = c(ci[,2], rev(ci[,3])))))
  }

  if (!(is.na(conf.level) || identical(args.cband, NA)) ) {
    args.cband1 <- list(col = SetAlpha(col, .2), border = NA)# modified by YF
    if (!is.null(args.cband))
      args.cband1[names(args.cband)] <- args.cband

    ci <- predict(x, interval="confidence", newdata=newx, level=conf.level) # Vertrauensband
    do.call("DrawBand", c(args.cband1, list(x = c(unlist(newx), rev(unlist(newx)))),
                          list(y = c(ci[,2], rev(ci[,3])))))
  }

  lines(y = fit, x = unlist(newx), col = col, lwd = lwd, lty = lty,
        type = type)
}
DrawBand <- function(x, y, col = SetAlpha("grey", 0.5), border = NA) {

  # accept matrices but then only n x y
  if(!identical(dim(y), dim(x))){
    x <- as.matrix(x)
    y <- as.matrix(y)

    if(dim(x)[2] == 1 && dim(y)[2] == 2)
      x <- x[, c(1,1)]
    else if(dim(x)[2] == 2 && dim(y)[2] == 1)
      y <- y[, c(1,1)]
    else
      stop("incompatible dimensions for matrices x and y")

    x <- c(x[,1], rev(x[,2]))
    y <- c(y[,1], rev(y[,2]))

  }

  # adds a band to a plot, normally used for plotting confidence bands
  polygon(x=x, y=y, col = col, border = border)
}
SetAlpha <- function(col, alpha=0.5) {

  if (length(alpha) < length(col)) alpha <- rep(alpha, length.out = length(col))
  if (length(col) < length(alpha)) col <- rep(col, length.out = length(alpha))

  acol <- substr(ColToHex(col), 1, 7)
  acol[!is.na(alpha)] <- paste(acol[!is.na(alpha)], DecToHex(round(alpha[!is.na(alpha)]*255,0)), sep="")
  acol[is.na(col)] <- NA
  return(acol)
}
DecToHex <- function(x) as.hexmode(as.numeric(x))
ColToHex <- function(col, alpha=1) {
  col.rgb <- col2rgb(col)
  col <- apply( col.rgb, 2, function(x) sprintf("#%02X%02X%02X", x[1], x[2], x[3]) )
  if(alpha != 1 ) col <- paste( col, DecToHex( round( alpha * 255, 0)), sep="")
  return(col)
  # old: sprintf("#%02X%02X%02X", col.rgb[1], col.rgb[2], col.rgb[3])
}
Pal <- function(pal, n=100, alpha=1) {

  if(missing(pal)) {
    res <- getOption("palette", default = structure(Pal("Helsana")[c(6,1:5,7:10)] ,
                     name = "Helsana", class = c("palette", "character")) )

  } else {

    palnames <- c("RedToBlack","RedBlackGreen","SteeblueWhite","RedWhiteGreen",
                  "RedWhiteBlue0","RedWhiteBlue1","RedWhiteBlue2","RedWhiteBlue3","Helsana","Tibco","RedGreen1",
                  "Spring","Soap","Maiden","Dark","Accent","Pastel","Fragile","Big","Long","Night","Dawn","Noon","Light")

    if(is.numeric(pal)){
      pal <- palnames[pal]
    }
    big <- c("#800000", "#C00000", "#FF0000", "#FFC0C0",
            "#008000","#00C000","#00FF00","#C0FFC0",
            "#000080","#0000C0", "#0000FF","#C0C0FF",
            "#808000","#C0C000","#FFFF00","#FFFFC0",
            "#008080","#00C0C0","#00FFFF","#C0FFFF",
            "#800080","#C000C0","#FF00FF","#FFC0FF",
            "#C39004","#FF8000","#FFA858","#FFDCA8")

    switch(pal
           , RedToBlack    = res <- colorRampPalette(c("red","yellow","green","blue","black"), space = "rgb")(n)
           , RedBlackGreen = res <- colorRampPalette(c("red", "black", "green"), space = "rgb")(n)
           , SteeblueWhite = res <- colorRampPalette(c("steelblue","white"), space = "rgb")(n)
           , RedWhiteGreen = res <- colorRampPalette(c("red", "white", "green"), space = "rgb")(n)
           , RedWhiteBlue0 = res <- colorRampPalette(c("red", "white", "blue"))(n)
           , RedWhiteBlue1 = res <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                                              "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(n)
           , RedWhiteBlue2 = res <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))(n)
           , RedWhiteBlue3 = res <- colorRampPalette(c(hred, "white", hblue))(n)
           , Helsana       = res <- c("rot"="#9A0941", "orange"="#F08100", "gelb"="#FED037"
                                       , "ecru"="#CAB790", "hellrot"="#D35186", "hellblau"="#8296C4", "hellgruen"="#B3BA12"
                                       , "hellgrau"="#CCCCCC", "dunkelgrau"="#666666", "weiss"="#FFFFFF")
           , Tibco         =  res <- apply( mcol <- matrix(c(
                                       0,91,0, 0,157,69, 253,1,97, 60,120,177,
                           156,205,36, 244,198,7, 254,130,1,
                           96,138,138, 178,113,60
                            ), ncol=3, byrow=TRUE), 1, function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
           , RedGreen1 =  res <- c(rgb(227,0,11, maxColorValue=255), rgb(227,0,11, maxColorValue=255),
                                  rgb(230,56,8, maxColorValue=255), rgb(234,89,1, maxColorValue=255),
                       rgb(236,103,0, maxColorValue=255), rgb(241,132,0, maxColorValue=255),
                       rgb(245,158,0, maxColorValue=255), rgb(251,184,0, maxColorValue=255),
                       rgb(253,195,0, maxColorValue=255), rgb(255,217,0, maxColorValue=255),
                       rgb(203,198,57, maxColorValue=255), rgb(150,172,98, maxColorValue=255),
                       rgb(118,147,108, maxColorValue=255))

           , Spring =  res <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3","#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
           , Soap =  res <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3","#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
           , Maiden =  res <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072","#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9","#BC80BD","#CCEBC5")
           , Dark =  res <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A","#66A61E", "#E6AB02", "#A6761D", "#666666")
           , Accent =  res <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99","#386CB0", "#F0027F", "#BF5B17", "#666666")
           , Pastel =  res <- c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4","#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2")
           , Fragile =  res <- c("#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4","#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC")
           , Big =  res <- big
           , Long =  res <- big[c(12,16,25,24,
                         2,11,6,15,18,26,23,
                         3,10,7,14,19,27,22,
                         4,8,20,28)]
           , Night =  res <- big[seq(1, 28, by=4)]
           , Dawn =  res <- big[seq(2, 28, by=4)]
           , Noon =  res <- big[seq(3, 28, by=4)]
           , Light = res <- big[seq(4, 28, by=4)]

           , GrandBudapest = res < c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
           , Moonrise1 = res <- c("#F3DF6C", "#CEAB07", "#D5D5D3", "#24281A")
           , Royal1 = res <- c("#899DA4", "#C93312", "#FAEFD1", "#DC863B")
           , Moonrise2 = res <- c("#798E87","#C27D38", "#CCC591", "#29211F")
           , Cavalcanti = res <- c("#D8B70A", "#02401B","#A2A475", "#81A88D", "#972D15")
           , Royal2 = res <- c("#9A8822", "#F5CDB4", "#F8AFA8", "#FDDDA0", "#74A089")
           , GrandBudapest2 = res <- c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
           , Moonrise3 = res <- c("#85D4E3", "#F4B5BD", "#9C964A", "#CDC08C", "#FAD77B")
           , Chevalier = res <- c("#446455", "#FDD262", "#D3DDDC", "#C7B19C")
           , Zissou = res <- c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")
           , FantasticFox = res <- c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20")
           , Darjeeling = res <- c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6")
           , Rushmore = res <- c("#E1BD6D", "#EABE94", "#0B775E", "#35274A", "#F2300F")
           , BottleRocket = res <- c("#A42820", "#5F5647", "#9B110E", "#3F5151", "#4E2A1E", "#550307", "#0C1707")
           , Darjeeling2 = res <- c("#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE",  "#000000")
    )

    attr(res, "name") <- pal
    class(res) <- append(class(res), "palette")

  }

  if(alpha != 1)
    res <- SetAlpha(res, alpha = alpha)

  return(res)

}

mylines=function(x, y, type = "l", ...) {    
    plot.xy(xy.coords(x[order(x)], y[order(x)]), type=type, ...)
}


# from princurve
whiskers <- function(x, s, ...) {
  graphics::segments(x[, 1], x[, 2], s[, 1], s[, 2], ...)
}
