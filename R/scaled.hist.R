scaled.hist = function(dat.ls, scale.factors, bin_width=100, cols=NULL, legend=NULL, cex.legend=1, xlim=NULL, ylim=NULL, ...) {

  # dat.ls is a list of lists
  # scale.factors is a list of data frames
  # bin_width is the width of the bins in um

  labels=names(dat.ls)
  names(scale.factors)=labels
  names(labels)=labels
  
  breaks = range(sapply(scale.factors, function(x) x$depth))
  breaks = seq(breaks[1], breaks[2]+bin_width, by=bin_width) # add bin_width to the end to include the last bin

  # normalize counts by area
  hist_info <- lapply(labels, function (b) {
    tmp = hist(dat.ls[[b]], breaks = breaks, plot = FALSE) # get histogram info without plotting
    
    m = merge(data.frame(depth=tmp$breaks[-length(tmp$breaks)], tmp$counts), scale.factors[[b]], all=T)
    
    if (!all(tmp$breaks[-length(tmp$breaks)]==m$depth)) {
      stop("scale.factors is longer than counts")
    }
    
    tmp$counts <- m$tmp.counts / m$area
    tmp$counts[is.na(tmp$counts)] = 0

    tmp
  })
  
  # plot histograms
  if (is.null(cols)) cols=1:length(labels)
  for (i in 1:length(labels)) {
    if (is.null(xlim)) xlim=range(breaks)
    if (is.null(ylim)) ylim=range(unlist(lapply(hist_info, function(a) a$counts)))
    if (i==1) {
      hist.col = col2rgb(cols[i])
      hist.col = rgb(hist.col[1], hist.col[2], hist.col[3], alpha=255*.15, maxColorValue=255)
      plot(hist_info[[i]], freq = TRUE, border="white",  col=hist.col, add=i>1, ylim=ylim, xlim=xlim , ...)
    }
    lines(hist_info[[i]]$mids, hist_info[[i]]$counts, type = "l", col = cols[i], lwd = 2)
  }
  mylegend(x=3, legend=if(is.null(legend)) labels else legend, lty=1, col=cols, y.intersp=1, text.width=700, cex=cex.legend)

  # return a data frame
  cbind(depth_start=breaks[-length(breaks)], sapply(labels, function(x) hist_info[[x]]$counts))
  
}

