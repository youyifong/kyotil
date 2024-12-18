# weightec Pearson correlation in the style of Hmisc
wtd.cor=function(x, y, weights=NULL, normwt=FALSE, na.rm=TRUE, method=c('unbiased', 'ML')) {
  var.x=wtd.var(x, weights=weights, normwt=normwt, na.rm=na.rm, method=method)
  var.y=wtd.var(y, weights=weights, normwt=normwt, na.rm=na.rm, method=method)
  mu.x=wtd.mean(x, weights=weights, normwt=normwt, na.rm=na.rm)
  mu.y=wtd.mean(y, weights=weights, normwt=normwt, na.rm=na.rm)
  wtd.mean((x-mu.x)*(y-mu.y), weights=weights, normwt=normwt, na.rm=na.rm)/sqrt(var.x*var.y)
}


# replace empty strings with NAs
mytable=function (..., exclude = if (useNA == "no") c(NA, NaN), useNA = "ifany", dnn = list.names(...), deparse.level = 1) {
  
  list.names <- function(...) {
    l <- as.list(substitute(list(...)))[-1L]
    if (length(l) == 1L && is.list(..1) && !is.null(nm <- names(..1))) 
      return(nm)
    nm <- names(l)
    fixup <- if (is.null(nm)) 
      seq_along(l)
    else nm == ""
    dep <- vapply(l[fixup], function(x) switch(deparse.level + 
                                                 1, "", if (is.symbol(x)) as.character(x) else "", 
                                               deparse(x, nlines = 1)[1L]), "")
    if (is.null(nm)) 
      dep
    else {
      nm[fixup] <- dep
      nm
    }
  }
  miss.use <- missing(useNA)
  miss.exc <- missing(exclude)
  useNA <- if (miss.use && !miss.exc && !match(NA, exclude, 
                                               nomatch = 0L)) 
    "ifany"
  else match.arg(useNA)
  doNA <- useNA != "no"
  if (!miss.use && !miss.exc && doNA && match(NA, exclude, 
                                              nomatch = 0L)) 
    warning("'exclude' containing NA and 'useNA' != \"no\"' are a bit contradicting")
  args <- list(...)
  if (length(args) == 1L && is.list(args[[1L]])) {
    args <- args[[1L]]
    if (length(dnn) != length(args)) 
      dnn <- paste(dnn[1L], seq_along(args), sep = ".")
  }
  if (!length(args)) 
    stop("nothing to tabulate")
  bin <- 0L
  lens <- NULL
  dims <- integer()
  pd <- 1L
  dn <- NULL
  for (a in args) {
    if (is.null(lens)) 
      lens <- length(a)
    else if (length(a) != lens) 
      stop("all arguments must have the same length")
    fact.a <- is.factor(a)
    if (doNA) 
      aNA <- anyNA(a)
    if (!fact.a) {
      a0 <- a
      op <- options(warn = 2)
      on.exit(options(op))
      a <- factor(a, exclude = exclude)
      options(op)
    }
    add.na <- doNA
    if (add.na) {
      ifany <- (useNA == "ifany")
      anNAc <- anyNA(a)
      add.na <- if (!ifany || anNAc) {
        ll <- levels(a)
        if (add.ll <- !anyNA(ll)) {
          ll <- c(ll, NA)
          TRUE
        }
        else if (!ifany && !anNAc) 
          FALSE
        else TRUE
      }
      else FALSE
    }
    if (add.na) 
      a <- factor(a, levels = ll, exclude = NULL)
    else ll <- levels(a)
    a <- as.integer(a)
    if (fact.a && !miss.exc) {
      ll <- ll[keep <- which(match(ll, exclude, nomatch = 0L) == 
                               0L)]
      a <- match(a, keep)
    }
    else if (!fact.a && add.na) {
      if (ifany && !aNA && add.ll) {
        ll <- ll[!is.na(ll)]
        is.na(a) <- match(a0, c(exclude, NA), nomatch = 0L) > 
          0L
      }
      else {
        is.na(a) <- match(a0, exclude, nomatch = 0L) > 
          0L
      }
    }
    nl <- length(ll)
    dims <- c(dims, nl)
    if (prod(dims) > .Machine$integer.max) 
      stop("attempt to make a table with >= 2^31 elements")
    dn <- c(dn, list(ll))
    bin <- bin + pd * (a - 1L)
    pd <- pd * nl
  }
  names(dn) <- dnn
  bin <- bin[!is.na(bin)]
  if (length(bin)) 
    bin <- bin + 1L
  y <- array(tabulate(bin, pd), dims, dimnames = dn)
  class(y) <- "table"
  y
  # table(..., exclude = exclude, useNA = "ifany", dnn = dnn, deparse.level = deparse.level)
}
  

# replace empty strings with NAs
empty2na=function(x) {x[x==""]=NA; x = as.factor(as.character(x))}

##  From the R-help e-mail by Ted Harding: http://tolstoy.newcastle.edu.au/R/e2/help/07/03/12853.html
##  See also http://tolstoy.newcastle.edu.au/R/help/05/05/4254.html
pava <- function(x, wt = rep(1, length(x)))
{
    n <- length(x)
    if (n <= 1) return(x)
    lvlsets <- 1:n
    repeat 
    {
        viol <- (as.vector(diff(x)) < 0)
        if (!(any(viol))) break
        i <- min((1:(n-1))[viol])
    
        lvl1 <- lvlsets[i]
        lvl2 <- lvlsets[i+1]
        ilvl <- ( (lvlsets == lvl1) | (lvlsets == lvl2) )
    
        x[ilvl] <- sum(x[ilvl] * wt[ilvl]) / sum(wt[ilvl])     
        lvlsets[ilvl] <- lvl1
    }
    x
} 


# from Cai Tianxi
## count how many YY's are smaller or equal to yy
N.L.E <- function(yy, YY)  ## sum I(YY <= yy[i])
{
   rank(c(yy+1e-8,YY))[1:length(yy)] - rank(yy)  ### add a small pertubation to avoid ties when calculating rank
}
N.L <- function(yy, YY)  ## sum I(YY < yy[i])
{
   rank(c(yy-1e-8,YY))[1:length(yy)] - rank(yy)  ### add a small pertubation to avoid ties when calculating rank
}
N.G.E <- function(yy, YY)  ## sum I(YY >= yy[i])
{
   length(YY)-(rank(c(yy-1e-8,YY))[1:length(yy)] - rank(yy))
}

# get first prinipal component
pr.1=function(x){
    x.s=scale(x)
    pr.s = prcomp(x.s)
    out=c(x.s %*% pr.s$rotation[,1])
    attr(out, "rotation")=pr.s$rotation[,1]
    out
}

# compute covariability as defined in Clarke 1995 page 2271, second formula
# x is a vector of characters, so is y
covariability=function(x,y){
    tmp=table(x)
    p.x=data.frame(tmp/(sum(tmp)))
    tmp=table(y)
    p.y=data.frame(tmp/(sum(tmp)))
    
    # tabulate all pairwise combination of aa
    pair = paste(x,y,sep="")
    tmp=table(pair)
    dat=data.frame(tmp/(sum(tmp)))
    
    ps=numeric(nrow(dat))
    for (i in 1:nrow(dat)) {
        row.=dat[i,]
        a2=as.character(row.[1,1])
        p.i=p.x$Freq[p.x$x==substr(a2,1,1)]
        p.j=p.y$Freq[p.y$y==substr(a2,2,2)]
        p.ij=row.[1,2]
        ps[i]=p.ij**2*log(p.ij/p.i/p.j)        
    }
    sum(ps)
}
## test
#dat=readFastaFile ("D:/1CHMM/domains/SET/SETpfamseed/seq/SETpfamseed_aligned.fasta")
#seqs.mat.a=mysapply(dat, s2c)
#covar=matrix(0,ncol(seqs.mat.a),ncol(seqs.mat.a))
#for (i in 1:(ncol(seqs.mat.a)-1)){
#    for (j in (i+1):ncol(seqs.mat.a)){
#        covar[i,j]<-covar[j,i]<-covariability(seqs.mat.a[,i], seqs.mat.a[,j])
#    }
#}
#sort(round(covar,3), decre=TRUE)[1:100]

# output format:
#+1 1:0.708333 2:1 3:1 4:-0.320755 5:-0.105023 6:-1 7:1 8:-0.419847 9:-1 10:-0.225806 12:1 13:-1 
# y: vector of 1 and non-1
# z: matrix or data.frame of covariates
# file.name: e.g. test. There should be no file extension
# ws: vector of weights
write.svm.input=function(y, z, file.name="svm_input", ws=NULL){
    rows=sapply(1:length(y), function(i){
        paste(c(ifelse(y[i]==1,"+1","-1"),(1:ncol(z))%.%":"%.%z[i,]), collapse=" ")
    })
    write(paste(rows,collapse="\n"), file = file.name%.%".dat", append = FALSE, sep = " ")        
    # write weight file
    if (!is.null(ws)) write(paste(ws,collapse="\n"), file = file.name%.%".wgt", append = FALSE, sep = " ")    
}

# m is data frame or matrix, each row is an observation
# w is a named vector, the names may not match m, may not have same number even
# return a vector
weighted.ave=function(m, w){
    if(is.data.frame(m)){
        nam=names(m)
        m=as.matrix(m)
    } else {
        nam=colnames(m)
    }
    ret=w[match(nam, names(w))]
    ret[is.na(ret)]=0
    names(ret)=nam
    c(m %*% ret )
}

# st.fit is an object returned by sn::st.mle, skew normal
get.st.mean=function(st.fit) {
    p=st.fit$dp
    
    xi=p["location"];
    alpha=p["shape"];
    df=p["df"]
    #omega=sd(X);
    omega=p["scale"]
    delta=alpha/sqrt(1+alpha^2);
    mu=delta*sqrt(df/pi)*gamma(0.5*(df-1))/gamma(0.5*df);
    
    # the 4 first moments
    moment1=xi+omega*mu;
    moment2=xi^2 + 2*omega*mu*xi + omega^2 * df/(df-2);
    moment3=xi^3 + 3*omega*mu*xi^2 + 3*omega^2*df/(df-2)*xi +
    omega^3*mu*(3-delta^2)*df/(df-3);
    moment4=xi^4 + 4*omega*mu*xi^3 + 6*omega^2*df/(df-2)*xi^2 +
    4*omega^3*mu*(3-delta^2)*df/(df-3)*xi+omega^4*3*df^2/((df-2)*(df-4));
    
    # the 4 useful measures
    mean=moment1;
    var=moment2;
    skew=moment3/var^(3/2);
    kurt=moment4/var^2 - 3;
    
    c(mean)
}


#iterated sum, like diff
summ=function(x) {
    x[-1]+x[-length(x)]
}

list_args <- Vectorize( function(a,b) c( as.list(a), as.list(b) ), SIMPLIFY = FALSE)
make_args_mtx <- function( alist ) Reduce(function(x, y) outer(x, y, list_args), alist)
multi.outer <- function(f, ... ) {
  args <- make_args_mtx(list(...))
  apply(args, 1:length(dim(args)), function(a) do.call(f, a[[1]] ) )
}


mylast = function (x, n=1, ...) {
    if (length(x)==1 & is.character(x)) tail (readLines(x), n=n, ...) # read file, length(x)==1 is needed b/c otherwise mylast(c("a","b")) won't behave properly
    else if (is.vector(x)) x[length(x)]
    else if (is.matrix(x)) x[nrow(x),]
    else if (is.array(x)) x[length(x)]
    else if (is.list(x)) x[[length(x)]]
    else stop ("mylast(): x not supported")
}


sample.for.cv=function(dat, v, seed){
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (inherits(save.seed,"try-error")) {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      
    set.seed(seed)
    n1=sum(dat$y==1)
    n2=sum(!dat$y==1)
    test.set = c(sample(which (dat$y==1), round(n1/v)), sample(which (dat$y==0), round(n2/v)))
    train.set = setdiff(1:nrow(dat), test.set)
    assign(".Random.seed", save.seed, .GlobalEnv)     
    list("train"=train.set, "test"=test.set)
}


get.sim.res = function(dir, res.name="res", verbose=TRUE) {
    
    if (!requireNamespace("abind")) {print("abind does not load successfully"); return (NULL) }
    
    ## get file names
    tmp=list.files(path=dir, pattern="batch[0-9]+.*.Rdata")
    if (length(tmp)==0) {
        # no files yet
        stop (paste0("no files in ",dir,"\n"))
    } else {
        fileNames = dir %.%"/"%.%tmp
    }
    
    # remove empty files, which sometimes appear, e.g. when disk is full
    empty.files=NULL
    for (x in fileNames) {
        if(file.info(x)$size==0) empty.files=c(empty.files, x)
    }
    fileNames=setdiff(fileNames, empty.files)
    
    ## read files
    res=lapply(fileNames, function(x) {
        if(verbose>=2) print(x); 
        load(file=x); 
        get(res.name)
    })
    # check dimension
    dims=lapply(res, dim)
    if (!is.null(dims[[1]])) {
        dim.same=sapply(dims, function(x) x==dims[[1]])
        if(!all(dim.same)) {
            print(dims)
            stop("not all files have the same dimensional results")
        }
    }
    
    res.all = do.call(abind::abind, res)
    names(dimnames(res.all))=names(dimnames(res[[1]]))
    
    if(verbose) str(res.all)
    
    if(verbose>=2) {
        cat("number of files: "%.%length(fileNames),"\n")
        print(dimnames(res.all)[[1]])
    }
    
    res.all    
}

MCsummary=function(dir, res.name="res", exclude.some=TRUE, exclude.col=1, verbose=TRUE){
    
    if(verbose) cat("\nsummarizing simulation results from ",dir,"\n",sep="")
    res.all=get.sim.res(dir,res.name,verbose=verbose)
    
    # exclude extreme values
    if(exclude.some) {
        q.1=0.05
        subset=res.all[exclude.col,1,]>quantile(res.all[exclude.col,1,],1-q.1,na.rm=T) 
        myprint(mean(subset)) # percent to be excluded
        res.all=res.all[,,!subset]
    }    
        
    mean.=      apply(res.all, 1:(length(dim(res.all))-1), function(x) mean(x, trim=0, na.rm=T))
    median.=    apply(res.all, 1:(length(dim(res.all))-1), function(x) median(x, na.rm=T)) 
    sd.=        apply(res.all, 1:(length(dim(res.all))-1), function(x) sd(x, na.rm=T)) 
    skew.=      apply(res.all, 1:(length(dim(res.all))-1), function(x) skew(x, na.rm=T))    
    width.mc=   apply(res.all, 1:(length(dim(res.all))-1), function(x) diff(quantile(x, c(.025,.975), na.rm=T))) # MC 95% interval
    # if res.all is two-dimensional, mean. etc is only one-dimensional, needs fix
    if(is.null(dim(mean.))) {dim(mean.)<-dim(median.)<-dim(sd.)<-dim(skew.)<-dim(width.mc)<-dim(bias)<-dim(rel.bias)<-dim(medbias)<-dim(rel.medbias)<-c(length(mean.),1)} # 
    
#    # coverage of MC.sd, not implemented for now
#    mc.sd=apply(res.all, 1:(length(dim(res.all))-1), function(x) sd(x, na.rm=T))[1:0]
#    cvg.mc.sd=apply(res.all[1:p,]>coef.0-1.96*mc.sd & res.all[1:p,]<coef.0+1.96*mc.sd, 1, mean, na.rm=T)
    
    list(mean=mean., median=median., sd=sd., skew=skew., width.mc=width.mc) 
        
}

# even if some of the n in nn have not results, it will result a row with NA for that n
# sum.sd determines whether we take mean of sd or median of sd, e.g. when computing summary of CI width
#exclude.some=T; verbose=T; coef.0=NULL; digit1=2; sum.est="mean"; sum.sd="median"; style=1; keep.intercept=FALSE
getFormattedMCSummary=function(path, sim, nn, fit.method, exclude.some=TRUE, exclude.col=1, verbose=TRUE, coef.0=NULL, digit1=2, sum.est=c("mean","median"), sum.sd=c("median","mean"),style=1, keep.intercept=FALSE) {
  
    sum.est<-match.arg(sum.est)    
    sum.sd<-match.arg(sum.sd)    
    names(nn)=nn
    stats=lapply(nn, function(n) {
        dir=path%.%"/"%.% sim%.%"_"%.%n%.%"_"%.%fit.method    
        try(MCsummary(dir, exclude.some=exclude.some, exclude.col=exclude.col, verbose=verbose))
    })
    subset.=!sapply(stats, function(x) inherits(x,"try-error"))
    stats=stats[subset.]
    nn=nn[subset.]    
    
    stat.names=mylast(dimnames(stats[[1]]$mean))[[1]] 
    sd.methods=sub("sd.","",stat.names[startsWith(stat.names,"sd.")  ])
    if (verbose) myprint(sd.methods)        
    
    # in the following objects, stats from different sample sizes are stacked one upon another
    # c() allows us to assign names 
    tab.1=cbind(
         mean.est=   c(do.call(rbind, lapply(stats, function (x) x$mean       [,"est",drop=F])))
        ,median.est= c(do.call(rbind, lapply(stats, function (x) x$median     [,"est",drop=F]))) 
        ,sd.est=     c(do.call(rbind, lapply(stats, function (x) x$sd         [,"est",drop=F]))) 
        ,skew.est=   c(do.call(rbind, lapply(stats, function (x) x$skew       [,"est",drop=F]))) 
        ,mcwdth.est= c(do.call(rbind, lapply(stats, function (x) x$width.mc   [,"est",drop=F]))) 
        ,              do.call(rbind, lapply(stats, function (x) x[[sum.sd]]  [,startsWith(stat.names,"sd."),drop=F]))
        ,              do.call(rbind, lapply(stats, function (x) x$mean       [,startsWith(stat.names,"covered."),drop=F]))
    )    
    p=nrow(tab.1)/length(nn)
    
    # add width of CI based on sd
    tmp=               do.call(rbind, lapply(stats, function (x) x[[sum.sd]]  [,startsWith(stat.names,"sd."),drop=F])) * 2 * 1.96
    colnames(tmp)=sub("sd.","wdth.",colnames(tmp))
    tab.1=cbind(tab.1, tmp)
    
    # compute bias based on coef.0
    if(!is.null(coef.0)) {
        tab.1=cbind(tab.1
            ,meanbias=        tab.1[,"mean.est"]-coef.0
            ,rel.meanbias=   (tab.1[,"mean.est"]-coef.0)/coef.0
            ,medianbias=      tab.1[,"median.est"]-coef.0
            ,rel.medianbias= (tab.1[,"median.est"]-coef.0)/coef.0
        )
    }
    
    # ro-oder the rows such that stats from different parameters are stacked one upon another    
    tab.2=apply(tab.1, 2, function(x) t(matrix(x,nrow=p)))
    
    #if (verbose) print(tab.2)
    
    if(style==1) {
        out=cbind(
              "$n$"=          c(rep(nn,p))
            , "est"=          c(formatDouble(tab.2[,sum.est%.%".est"] ,digit1))
            , "est(\\%bias)"= c(formatDouble(tab.2[,sum.est%.%".est"] ,digit1)%.%"("%.%round( tab.2[,"rel."%.%sum.est%.%"bias"] *100)%.%")")
            , "range"=        c(formatDouble(tab.2[,"mcwdth.est"],digit1))
            # put width and coverage prob together, e.g. 0.17(93)
            ,                 matrix(formatDouble(tab.2[,startsWith(colnames(tab.2),"wdth.")],digit1) %.% ifelse (is.nan(tab.2[,startsWith(colnames(tab.2),"covered.")]),"", "("%.%round( tab.2[,startsWith(colnames(tab.2),"covered.")] *100)%.%")"), 
                                  nrow=p*length(nn), dimnames=list(NULL,sd.methods))
        ) 
        rownames(out)=outer(nn, names(coef.0), paste)
        if(!keep.intercept) out=out[-(1:length(nn)),] 
    
    } else stop("style not supported")
    
    out
        
}


#Need to re-implement this more efficiently!!!
#    library(maptools)
nearestPointOnSegment=function (s, p) {
    ap = c(p[1] - s[1, 1], p[2] - s[1, 2])
    ab = c(s[2, 1] - s[1, 1], s[2, 2] - s[1, 2])
    t = sum(ap * ab)/sum(ab * ab)
    t = ifelse(t < 0, 0, ifelse(t > 1, 1, t))
    x = s[1, 1] + ab[1] * t
    y = s[1, 2] + ab[2] * t
    result = c(x, y, sqrt((x - p[1])^2 + (y - p[2])^2))
    names(result) = c("X", "Y", "distance")
    result
}


# neither project_to_curve or earlier version of get.lam works correctly to project,  nor pcurve package works
predict.pcc=function(object, newdat, ...) {
    s=object$s[object[["ord"]],]
    res=sapply (1:nrow(newdat), function(i) {
        tmp=sapply (2:nrow(s), function(j) {
            fit=nearestPointOnSegment(s[(j-1):j,], newdat[i,])
        })
        fit=tmp[,which.min(tmp[3,])]
        #segments(newdat[i,1], newdat[i,2], fit[1], fit[2])
    })
    #sum(res["distance",])
    res
}


# rank-based inverse normal transformation
rank.inv.norm = function(x) {
    qnorm(rank(x)/(1+length(x)))
}

INT=rank.inv.norm
    
# covert a decimal number to a binary representation with d digits
dec_to_binary <- function(x,d) sapply(x, function(xx) ifelse(is.na(xx), NA, paste(rev(as.integer(intToBits(xx))[1:d]), collapse="")))
