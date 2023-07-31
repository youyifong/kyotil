# collapse Wstratum if there are empty cells
# n.demo is the number of demographics strata
cove.boost.collapse.strata = function(dat.b, n.demo) {
  
  # dat.b is expected to contain these columns
  #   Ptid
  #   sampling_bucket, which holds demo strata
  #   ph2
  #   Wstratum
  #   CalendarBD1Interval
  #   sampling_bucket_formergingstrata
  
  # this will ensure that there are always two columns even if there are no ph2 samples or there are only ph2 samples
  # make a copy of ph2 because we don't want to change ph2 into a factor. otherwise we need to change it back in several places
  dat.b$ph2f=factor(dat.b$ph2, levels=c(FALSE,TRUE)) 
  
  tab=with(dat.b, table(Wstratum, ph2f)); tab
  
  if (all(tab[,2]!=0)) return (dat.b)
  
  # collapsing is 3-step process
  
  # 1. do it across demo strata within each sampling bucket 
  sampling_buckets = unique(dat.b$sampling_bucket)
  for (i in sampling_buckets) {
    select = dat.b$sampling_bucket==i 
    
    # make sure there are such ph1 samples
    if (sum(select, na.rm=T)>0) {
      
      # make sure there is at least 1 ph2 sample and at least 2 Wstratum
      if (sum(dat.b$ph2f[select]==TRUE, na.rm=T)>=1 & length(unique(dat.b[select,"Wstratum"]))>=2) {
        tab = with(dat.b[select,], table(Wstratum, ph2f)); tab
        # merging
        if (ncol(tab)>1) { # if all samples are ph2, tab is only 1-column
          if (any(tab[,2]==0)) {
            dat.b[select,"Wstratum"] = min(dat.b[select,"Wstratum"], na.rm=T)
          }
        }
      } else {
        # if there are no ph2 samples, will collapse in the next level
      }
      
    } # okay if there are no samples in the bucket
  }
  tab=with(dat.b, table(Wstratum, ph2f)); tab
  
  if (all(tab[,2]!=0)) return (dat.b)
  
  # 2. do it across the 4 calendar periods
  # merge a period with the next period if not the last, merge with the last period with the previous if needed
  sampling_buckets.2 = unique(dat.b$sampling_bucket_formergingstrata)
  # need this in the merging process as it needs to be updated
  dat.b$tmpCalendarBD1Interval=dat.b$CalendarBD1Interval
  
  for (i in sampling_buckets.2) {
    select = dat.b$sampling_bucket_formergingstrata==i 
    
    # make sure there are such ph1 samples
    if (sum(select, na.rm=T)>0) {
      
      # make sure there is at least 1 ph2 sample
      if (sum(dat.b$ph2f[select]==TRUE, na.rm=T)>=1) {
        
        tab = with(dat.b[select,], table(tmpCalendarBD1Interval, ph2f)); tab
        jj = which(tab[,2]==0)
        while (length(jj)>0) {
          j=jj[1]
          if (j==nrow(tab)) {
            # the last is empty, set to the previous
            chng=-1
          } else {
            # set to the next
            chng=1
          }
          
          select2 = select & dat.b$tmpCalendarBD1Interval==(rownames(tab)[j])
          
          # merge Wstratum
          dat.b[select2,"Wstratum"] = dat.b[select2,"Wstratum"] + n.demo * chng
          # need to update sampling bucket so that step 3 can work
          dat.b[select2,"sampling_bucket"] = dat.b[select2,"sampling_bucket"] + chng
          
          # re-tabulate after updating tmpCalendarBD1Interval
          dat.b$tmpCalendarBD1Interval[select2] = dat.b$tmpCalendarBD1Interval[select2] + chng
          tab = with(dat.b[select,], table(tmpCalendarBD1Interval, ph2f)); tab
          jj = which(tab[,2]==0)
          
        } # end while
        
      } else {
        stop("something wrong - no ph2 samples in sampling_buckets.2 "%.%i)
      }
    } else {
      # no ph1 samples in this sampling_buckets.2, probably case
    }
    
  }
  tab=with(dat.b, table(Wstratum, ph2f)); tab
  
  if (all(tab[,2]!=0)) return (dat.b)
  
  # 3. merge across demo strata within each sampling bucket one more time 
  #       because step 2 collapsing time periods may empty demo strata
  sampling_buckets = unique(dat.b$sampling_bucket)
  for (i in sampling_buckets) {
    
    select = dat.b$sampling_bucket==i 
    
    # make sure there are such ph1 samples
    if (sum(select, na.rm=T)>0) {
      
      # make sure there is at least 1 ph2 sample and at least 2 Wstratum
      if (sum(dat.b$ph2[select]==TRUE, na.rm=T)>=1 & length(unique(dat.b[select,"Wstratum"]))>=2) {
        tab = with(dat.b[select,], table(Wstratum, ph2f)); tab
        # merging
        if (ncol(tab)>1) { # if all samples are ph2, tab is only 1-column
          if (any(tab[,2]==0)) {
            dat.b[select,"Wstratum"] = min(dat.b[select,"Wstratum"], na.rm=T)
          }
        }
      } else {
        # if there are no ph2 samples, will collapse in the next level
      }
    } # okay if there are no samples in the bucket
  }
  tab=with(dat.b, table(Wstratum, ph2f)); tab
  # make sure no empty cells after all this
  if(!all(tab[,2]>0)) {
    stop("error in cove.boost.collapse.strata, not all ph2 cells are positive !!!!!!!!!!!!!!!!!!")
  }
  

  return (dat.b)
  
}




# extract assay name from marker names, which include Day, e.g.
marker.name.to.assay=function(a) {
  
  if (startsWith(a,"Day")) {
    # Day22pseudoneutid50 => pseudoneutid50 
    sub("Day[[0123456789]+", "", a)
    
  } else if (startsWith(a,"BD")) {
    # BD29pseudoneutid50 => pseudoneutid50
    sub("BD[[0123456789]+", "", a)
    
  } else if (contain(a,"overBD1")) {
    # DeltaBD29overBD1pseudoneutid50 => pseudoneutid50
    sub("DeltaBD[[0123456789]+overBD1", "", a)    
    
  } else if (contain(a,"overB")) {
    # Delta22overBpseudoneutid50 => pseudoneutid50
    sub("Delta[[0123456789]+overB", "", a)
    
  } else if (contain(a,"over")) {
    sub("Delta[[0123456789]+over[[0123456789]+", "", a)
    
  } else stop("marker.name.to.assay: not sure what to do")
}


draw.x.axis.cor=function(xlim, llox, llox.label, for.ggplot=FALSE){
  
  xx=seq(ceiling(xlim[1]), floor(xlim[2]))        
  if (is.na(llox)) {
    labels = sapply (xx, function(x) if (x>=3) bquote(10^.(x)) else 10^x )
    
  } else if (llox.label=="delta") {
    labels = sapply (xx, function(x) if (x>=3 | x<=-3) bquote(10^.(x)) else 10^x )
    
  } else {
    
    xx=xx[xx>log10(llox*1.8)]
    labels = sapply (xx, function(x) if(x>=3) bquote(10^.(x)) else 10^x)
    xx=c(log10(llox), xx)
    labels=c(llox.label, labels)
  }
  
  # add e.g. 30 between 10 and 100
  if (length(xx)<4) {
    for (i in 1:length(xx)) {
      x=xx[i]
      xx=c(xx, x+log10(3))
      labels=c(labels, if (x>=3) bquote(3%*%10^.(x)) else 3*10^(x) )
    }

    labels=c(if (min(xx)-1>=3) bquote(3%*%10^.(min(xx)-1)) else 3*10^(min(xx)-1), labels)    
    xx=c(min(xx)-1+log10(3), xx)
  }
  print(xx)

  if (for.ggplot) {
    return(list(ticks = xx, labels = labels))
  } else {
    axis(1, at=xx, labels=sapply(labels, function (label) as.expression(label)))
  }
  
}


# get plotting range
get.xlim=function(dat, marker) {
  assay=marker.name.to.assay(a)
  
  # the default
  ret=range(dat[[marker]], log10(lloxs[assay]/2), na.rm=T)
  
  # may be customized, e.g. to have the same xlim for different variants in the same type of assay
  # if (TRIAL=="moderna_boost") {
  #   if(assay %in% c("bindSpike", "bindRBD")) {
  #     ret=range(dat[["Day"%.%time%.%"bindSpike"]], 
  #               dat[["Day"%.%time%.%"bindRBD"]], 
  #               log10(lloxs[c("bindSpike","bindRBD")]/2), na.rm=T)
  #     
  #   } 
  # }

  delta=(ret[2]-ret[1])/20     
  c(ret[1]-delta, ret[2]+delta)
}


