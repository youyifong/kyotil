# collapse Wstratum if there are empty cells
cove.boost.collapse.strata = function(dat.b, n.demo) {
  
  # dat.b is expected to contain these columns
  #   Ptid
  #   sampling_bucket, which holds demo strata
  #   ph2
  #   Wstratum
  #   CalendarBD1Interval
  #   sampling_bucket_formergingstrata
  
  tab=with(dat.b, table(Wstratum, ph2)); tab
  
  # collapsing is 3-step process
  
  # 1. do it across demo strata within each sampling bucket 
  sampling_buckets = unique(dat.b$sampling_bucket)
  for (i in sampling_buckets) {
    select = dat.b$sampling_bucket==i 
    
    # make sure there are such ph1 samples
    if (sum(select, na.rm=T)>0) {
      
      # make sure there is at least 1 ph2 sample and at least 2 Wstratum
      if (sum(dat.b[select,"ph2"], na.rm=T)>=1 & length(unique(dat.b[select,"Wstratum"]))>=2) {
        tab = with(dat.b[select,], table(Wstratum, ph2)); tab
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
  tab=with(dat.b, table(Wstratum, ph2)); tab
  
  
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
      if (sum(dat.b[select,"ph2"], na.rm=T)>=1) {
        
        tab = with(dat.b[select,], table(tmpCalendarBD1Interval, ph2)); tab
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
          tab = with(dat.b[select,], table(tmpCalendarBD1Interval, ph2)); tab
          jj = which(tab[,2]==0)
          
        } # end while
        
      } else {
        stop("something wrong - no ph2 samples in sampling_buckets.2 "%.%i)
      }
    } else {
      # no ph1 samples in this sampling_buckets.2, probably case
    }
    
  }
  tab=with(dat.b, table(Wstratum, ph2)); tab
  
  
  # 3. merge across demo strata within each sampling bucket one more time 
  #       because step 2 collapsing time periods may empty demo strata
  sampling_buckets = unique(dat.b$sampling_bucket)
  for (i in sampling_buckets) {
    
    select = dat.b$sampling_bucket==i 
    
    # make sure there are such ph1 samples
    if (sum(select, na.rm=T)>0) {
      
      # make sure there is at least 1 ph2 sample and at least 2 Wstratum
      if (sum(dat.b[select,"ph2"], na.rm=T)>=1 & length(unique(dat.b[select,"Wstratum"]))>=2) {
        tab = with(dat.b[select,], table(Wstratum, ph2)); tab
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
  tab=with(dat.b, table(Wstratum, ph2)); tab
  # make sure no empty cells after all this
  if(!all(tab[,2]>0)) {
    stop("error in cove.boost.collapse.strata, not all ph2 cells are positive !!!!!!!!!!!!!!!!!!")
  }
  
  return (dat.b)

}

