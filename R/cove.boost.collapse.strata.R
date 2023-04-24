# collapse Wstratum if there are empty cells
cove.boost.collapse.strata = function(dat.ph1, n.demo) {
  
  # dat.ph1 is expected to contain these columns
  #   Ptid
  #   sampling_bucket, which holds demo strata
  #   ph2
  #   Wstratum
  #   CalendarBD1Interval
  #   sampling_bucket_formergingstrata

  # collapsing is 3-step process
  
  # 1. do it across demo strata within each sampling bucket 
  
  sampling_buckets = 0:max(dat.ph1$sampling_bucket, na.rm=T)
  
  for (i in sampling_buckets) {
    
    select = dat.ph1$sampling_bucket==i 
    
    # make sure there are such ph1 samples
    if (sum(select, na.rm=T)>0) {
      
      # make sure there is at least 1 ph2 sample
      if (sum(dat.ph1[select,"ph2"], na.rm=T)>=1) {
        tab = with(dat.ph1[select,], table(Wstratum, ph2))
        # merging
        if (any(tab[,2]==0)) {
          dat.ph1[select,"Wstratum"] = min(dat.ph1[select,"Wstratum"], na.rm=T)
        }
      } else {
        # if there are no ph2 samples, will collapse in the next level
      }
    } # okay if there are no samples in the bucket
  }
  
  
  # 2. do it across the 4 calendar periods
  # merge a period with the next period if not the last, merge with the last period with the previous if needed
  
  sampling_buckets.2 = 0:max(dat.ph1$sampling_bucket_formergingstrata, na.rm=T)
  # need this in the merging process as it needs to be updated
  dat.ph1$tmpCalendarBD1Interval=dat.ph1$CalendarBD1Interval
  
  for (i in sampling_buckets.2) {
    select = dat.ph1$sampling_bucket_formergingstrata==i 
    
    # make sure there are such ph1 samples
    if (sum(select, na.rm=T)>0) {
      
      # make sure there is at least 1 ph2 sample
      if (sum(dat.ph1[select,"ph2"], na.rm=T)>=1) {
        tab = with(dat.ph1[select,], table(tmpCalendarBD1Interval, ph2))
        # merging
        jj = which(tab[,2]==0)
        for (j in jj) {
          select2 = select & dat.ph1$tmpCalendarBD1Interval==j
          # if the last is empty, set all to the first
          dat.ph1[select2,"Wstratum"] = dat.ph1[select2,"Wstratum"] + n.demo * ifelse(j<4, 1, -3)
          # also need to update sampling bucket so that step 3 can work
          dat.ph1[select2,"sampling_bucket"] = dat.ph1[select2,"sampling_bucket"] + ifelse(j<4, 1, -3)
          # need to update this so that select2 will pick these up
          dat.ph1$tmpCalendarBD1Interval[select2] = dat.ph1$tmpCalendarBD1Interval[select2] + 1
        }
      } else {
        stop("something wrong - no ph2 samples in sampling_buckets.2 "%.%i)
      }
    } else {
      # no ph1 samples in this sampling_buckets.2, probably case
    }
  }
  
  
  # 3. merge across demo strata within each sampling bucket one more time 
  #       because step 2 collapsing time periods may empty demo strata
  
  sampling_buckets = 0:max(dat.ph1$sampling_bucket, na.rm=T)
  for (i in sampling_buckets) {
    
    select = dat.ph1$sampling_bucket==i 
    
    # make sure there are such ph1 samples
    if (sum(select, na.rm=T)>0) {
      
      # make sure there is at least 1 ph2 sample
      if (sum(dat.ph1[select,"ph2"], na.rm=T)>=1) {
        tab = with(dat.ph1[select,], table(Wstratum, ph2))
        (tab)
        # merging
        if (any(tab[,2]==0)) {
          dat.ph1[select,"Wstratum"] = min(dat.ph1[select,"Wstratum"], na.rm=T)
        }
      } else {
        # if there are no ph2 samples, will collapse in the next level
      }
    } # okay if there are no samples in the bucket
  }
  
  
  # make sure no empty cells after all this
  tab=with(dat.ph1, table(Wstratum, ph2)); tab
  stopifnot(all(tab[,2]>0))
  
  return (dat.ph1)

}

