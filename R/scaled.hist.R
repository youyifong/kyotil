scaled.hist = function(dat.ls, scale.factors, bin_width=100, cols=NULL, legend=NULL, cex.legend=1, xlim=NULL, ylim=NULL, span=NULL, add.minus.50=F, hist.to.show=2, ...) {

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
  # print(hist_info) # if counts are 0, the code below will complain: Error in plot.window(xlim, ylim, "", ...) : need finite 'ylim' values
  
  # plot histograms
  if (is.null(cols)) cols=1:length(labels)
  
  for (i in 1:length(labels)) {
    
    if (is.null(xlim)) xlim=c(min(breaks) - ifelse(add.minus.50, 50, 0), max(breaks))
    if (is.null(ylim)) ylim=range(unlist(lapply(hist_info, function(a) a$counts)))
    
    hist.col = col2rgb(cols[i])
    # make visible only for the hist.to.show element of the list
    if (i!=hist.to.show) hist.col=col2rgb(NA)
    hist.col = rgb(hist.col[1], hist.col[2], hist.col[3], alpha=255*.15, maxColorValue=255)
    plot(hist_info[[i]], freq = TRUE, border=if (i==hist.to.show) "white" else NA,  col=hist.col, add=i>1, ylim=ylim, xlim=xlim, ...)

    if (is.null(span)) {
      # non-smoothed version
      lines(hist_info[[i]]$mids, hist_info[[i]]$counts, type = "l", col = cols[i], lwd = 2)      
    } else {
      # smoothed version
      dat.tmp = data.frame(mids = hist_info[[i]]$mids, counts = hist_info[[i]]$counts)
      
      # optionally, add add (-50,0) as the first point
      if(add.minus.50) {
        if (!all(hist_info[[i]]$counts==Inf))
          dat.tmp=rbind(data.frame(mids=min(hist_info[[i]]$mids)-100, counts=0), dat.tmp)
      }
      # print(dat.tmp)

      smooth_fit = loess(counts ~ mids, data = dat.tmp, span = span)
      # use 90 instead of 100 below so remove some edge artifacts
      smooth_x = seq(min(hist_info[[i]]$mids) - ifelse(add.minus.50, 90, 0), max(hist_info[[i]]$mids), length.out = 500)
      lines(smooth_x, pmax(predict(smooth_fit, newdata = data.frame(mids = smooth_x)), 0), type = "l", col = cols[i], lwd = 2)
    }

  }
  if (!is.null(legend)) {
    mylegend(x=3, legend=legend, lty=1, col=cols, y.intersp=1, text.width=700, cex=cex.legend, lwd=1.8)
  }

  # return a data frame
  cbind(depth_start=breaks[-length(breaks)], sapply(labels, function(x) hist_info[[x]]$counts))
  
}

