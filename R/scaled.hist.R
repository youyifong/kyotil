scaled.hist = function(dat.ls, scale.factors, bin_width=100, patch_size=5, cols=NULL, legend=NULL, cex.legend=1, ...) {

# dat.ls is a list of data frames
# bin_width is the width of the bins in um
# scale.factors is a list of data frames

  labels=names(dat.ls)
  names(labels)=labels
  names(scale.factors)=labels

  min_break <- floor(min(unlist(dat.ls)) / bin_width) * bin_width
  max_break <- ceiling(max(unlist(dat.ls)) / bin_width) * bin_width
  breaks <- seq(min_break, max_break, by = bin_width)
  nbin_total=length(breaks)-1
  
  # get histogram info without plotting
  hist_info <- lapply(labels, function (a) hist(dat.ls[[a]], breaks = breaks, plot = FALSE) )
  
  # normalize counts against size of biopsy at each depth
  for (a in labels) {
    # nbin_total may be greater than nbin
    nbin = length(scale.factors[[a]]$depth)
    
    # make sure the depth and the breaks match up
    stopifnot(all(hist_info[[a]]$breaks[1:nbin]==scale.factors[[a]]$depth[1:nbin]))
    if (nbin_total>nbin) stopifnot(all(hist_info[[a]]$counts[(nbin+1):nbin_total]==0))

    # get average cell count per patch
    hist_info[[a]]$counts <- hist_info[[a]]$counts / c(scale.factors[[a]]$count, rep(1,nbin_total-nbin))
    # there are 1600 patches in every 100x100 um^2 area if each patch is 5px by 5px.
    hist_info[[a]]$counts <- hist_info[[a]]$counts * (bin_width/patch_size)^2
  }

  # plot histograms
  if (is.null(cols)) cols=1:length(labels)
  for (i in 1:length(labels)) {
    hist.col = col2rgb(cols[i])
    hist.col = rgb(hist.col[1], hist.col[2], hist.col[3], alpha=255*.15, maxColorValue=255)
    plot(hist_info[[i]], freq = TRUE, border="white",  col=hist.col, add=i>1, ylim=range(unlist(lapply(hist_info, function(a) a$counts))), ...)
  }
  mylegend(x=3, legend=if(is.null(legend)) labels else legend, lty=1, col=cols, y.intersp=1, text.width=700, cex=cex.legend)

  # return a data frame
  cbind(depth_start=hist_info[[1]]$breaks[1:nbin_total], sapply(labels, function(x) hist_info[[x]]$counts))
  
}

