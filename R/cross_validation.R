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
