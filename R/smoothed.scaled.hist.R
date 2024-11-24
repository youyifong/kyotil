smoothed.scaled.hist = function(dat.ls, bin_width, scale.factors=NULL, cols=NULL, legend=NULL, cex.legend=1, ...) {
  
  # bin_width <- 100
  min_break <- floor(min(unlist(dat.ls)) / bin_width) * bin_width
  max_break <- ceiling(max(unlist(dat.ls)) / bin_width) * bin_width
  breaks <- seq(min_break, max_break, by = bin_width)
  nbin_total=length(breaks)-1
  
  labels=names(dat.ls)
  names(labels)=labels
  
  if (is.null(scale.factors)) {
    nn = lapply(dat.ls, function(x) length(x))
    # default scale.factors are to make density comparable
    scale.factors = nn[[1]]/nn
  }
  names(scale.factors)=labels
  
  if (is.null(cols)) cols=1:length(dat.ls)
  
  # get histogram info without plotting
  hist_info <- lapply(labels, function (a) hist(dat.ls[[a]], breaks = breaks, plot = FALSE) )
  
  # normalize counts against size of biopsy at each depth
  for (a in labels) {
    if (length(scale.factors[[a]])==1) {
      # a single scale factor for all bins
      hist_info[[a]]$counts <- hist_info[[a]]$counts * scale.factors[[a]]
    
    } else {
      # make sure the depth and the breaks match up
      nbin = length(scale.factors[[a]]$depth)
      
      stopifnot(all(hist_info[[a]]$breaks[1:nbin]==scale.factors[[a]]$depth[1:nbin]))
      if (nbin_total>nbin) stopifnot(all(hist_info[[a]]$counts[(nbin+1):nbin_total]==0))
      # there are 1600 patches in every 100x100 um^2 area. 
      hist_info[[a]]$counts <- hist_info[[a]]$counts / (c(scale.factors[[a]]$count, rep(1,nbin_total-nbin)) / 1600)
      
    }
  }
  
  # # estimate density curve
  # densities <- lapply(labels, function (a) density(dat.ls[[a]]) )
  # for (a in labels) {
  #   # Scale the height of the density curve to match the histogram's frequency
  #   densities[[a]]$y <- densities[[a]]$y * diff(hist_info[[a]]$mids[1:2]) * length(dat.ls[[a]])
  #   # normalize the density curve height to 100 um width
  #   densities[[a]]$y <- densities[[a]]$y * scale.factors[[a]]
  # }
  
  
  hist.col = col2rgb(cols[1])
  hist.col = rgb(hist.col[1], hist.col[2], hist.col[3], alpha=255*.15, maxColorValue=255)

  for (i in 1:length(labels)) plot(hist_info[[i]], freq = TRUE, border="white",  col=hist.col, ...)
  # for (i in 1:length(labels)) lines(densities[[i]], col = cols[i])
  
  # mylegend(x=3, legend=if(is.null(legend)) labels else legend, lty=1, col=cols, y.intersp=1, text.width=700, cex=cex.legend)
  
  hist_info
  
}