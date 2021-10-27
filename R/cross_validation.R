# k-fold cross-validation (CV) #
get.kfold.splits=function(dat, k, seed) {
  # save rng state before set.seed in order to restore before exiting this function
  save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
  if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      
  set.seed(seed)    
  
  n0=nrow(dat$control)
  n1=nrow(dat$case)
  
  # k-fold CV #
  training.subsets=list()
  test.subsets=list()
  tmp1=sample(1:n1)
  tmp0=sample(1:n0)
  splits=list()
  for (ki in 1:k) {
    splits[[ki]]=list(training=list(case=tmp1[(1:n1)%%k!=ki-1], control=tmp0[(1:n0)%%k!=ki-1]),
                      test=list(case=tmp1[(1:n1)%%k==ki-1], control=tmp0[(1:n0)%%k==ki-1]))
  }
  splits
    
  # restore rng state 
  assign(".Random.seed", save.seed, .GlobalEnv)     
  splits
}


# k-fold CV #
kfold.split=function(k, n1, n0){
  training.subsets=list()
  test.subsets=list()
  tmp1=sample(1:n1)
  tmp0=sample(1:n0)
  splits=list()
  for (ki in 1:k) {
    splits[[ki]]=list(training=list(case=tmp1[(1:n1)%%k!=ki-1], control=tmp0[(1:n0)%%k!=ki-1]),
                      test=list(case=tmp1[(1:n1)%%k==ki-1], control=tmp0[(1:n0)%%k==ki-1]))
  }        
  splits
}

# Random 4:1 #
ran.kfold.split=function(k, n1, n0, replicates){
  training.subsets=list()
  test.subsets=list()
  splits=list()
  for (r in 1:replicates) {
    tmp1=sample(1:n1)
    tmp0=sample(1:n0)
    splits[[r]]= list(training=list(case=tmp1[(1:n1)%%k!=1], control=tmp0[(1:n0)%%k!=1]),
                      test=list(case=tmp1[(1:n1)%%k==1], control=tmp0[(1:n0)%%k==1]))
  }        
  splits
}

# Leave pair out CV #
lpo.split=function(n1, n0){
  training.subsets=list()
  test.subsets=list()
  splits=list()
  ind=0
  for (i in 1:n1) {
    for (j in 1:n0) {
      ind=ind+1
      splits[[ind]]=list(training=list(case=setdiff(1:n1,i), control=setdiff(1:n0,j)),
                         test=    list(case=i,               control=j))
    }        
  }        
  splits
}

# CV schemes #
get.splits=function(dat, cv.scheme=c("LPO","5fold","50xrandom4:1"), seed) {
  # save rng state before set.seed in order to restore before exiting this function
  save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
  if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      
  set.seed(seed)    
  
  n0=nrow(dat$control)
  n1=nrow(dat$case)
  
  cv.scheme<-match.arg(cv.scheme)  
  if (cv.scheme=="LPO") {
    splits=lpo.split(n1, n0)
  } else if (cv.scheme=="5fold") {
    splits=kfold.split(5, n1, n0)
  } else if (cv.scheme=="50xrandom4:1") {
    splits=ran.kfold.split(5, n1, n0, 50)
  }
  # restore rng state 
  assign(".Random.seed", save.seed, .GlobalEnv)     
  splits
}
