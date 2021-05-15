cbinduneven=function(li) {
    # bind a list of data frame or named vector that are not of the same length
    allnames=lapply(li, rownames)
    alllen=lapply(allnames, length)
    #nams = allnames[[which.max(alllen)]]
    nams = allnames[[1]]# so that the rows are added consecutively
    nams= c(nams, setdiff(unique(unlist(allnames)), nams)) # append additional names
    if (any(nams=="")) stop("cbinduneven: empty rownames are not allowed:\n"%.%concatList(nams,"|"))
    #myprint(nams)
    
    res=NULL
    for (i in 1:length(li)){
        a=li[[i]]
        if (is.table(a)) a=as.matrix(a,ncol=1)
        p=ncol(a)
        toadd = matrix(NA, nrow=length(nams), ncol=p, dimnames=list(nams,colnames(a)))
        toadd[rownames(a),]=as.matrix(a)
        tmp=as.data.frame(toadd, stringsAsFactors=FALSE)
        if(ncol(tmp)==1) names(tmp)=names(li)[i] else names(tmp)=names(li)[i]%.%names(tmp) 
        if (i==1) {
            res=tmp # cbind NULL and tmp doesn't work
        } else {
            res=cbind(res,tmp)
        }
    }
    res
}



# returns binary representation of an integer
binary<-function(i) if (i) paste(binary(i %/% 2), i %% 2, sep="") else "" 

# returns binary representatin of an integer with leading 0, the length of string is n
binary2<-function(i, n) {
    a<-2^((n-1):0)
    b<-2*a
    sapply(i,function(x) paste(as.integer((x %% b)>=a),collapse=""))
} 

unix=function (){
    substr(Sys.getenv("R_HOME") , 1,1)=="/"
}


#mysystem can call any exe file
mysystem = function (cmd, ...) {
    system ( paste(Sys.getenv("COMSPEC")," /c ",cmd) , invisible =TRUE, intern=FALSE, ...)
}


# convert temp from f to c
f2c = function (f) {
    (f-32)*5/9
}

# convert a factor to integer using its value, e.g. 1 to 1, 2 to 2
ftoi = function (f) {
    as.integer (as.character (f) )
}


# DONT CHANGE THIS, USED By RUMINEX!
# if ret.mat is set to true, always return a matrix
# in the output, each row corresponds to one element of X, instead of each column
mysapply=function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, ret.mat=TRUE) 
{
    if (is.null(names(X)) & is.numeric(X)) names(X)=X%.%""
    FUN <- match.fun(FUN)
    answer <- lapply(X, FUN, ...)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
        names(answer) <- X
    if (simplify && length(answer) && length(common.len <- unique(unlist(lapply(answer, 
        length)))) == 1) {
         if (common.len >= 1) 
            if (common.len == 1 & !ret.mat) 
                unlist(answer, recursive = FALSE)
            else 
                t(array(unlist(answer, recursive = FALSE), dim = c(common.len, 
                    length(X)), dimnames = if (!(is.null(n1 <- names(answer[[1]])) & 
                    is.null(n2 <- names(answer)))) 
                    list(n1, n2)))
        else t(answer)
    }
    else t(answer)
}
## test mysapply
#sapply(1:3, function (i) if (i==2) rep(NA,2) else 1:3 )
#mysapply(1:3, function (i) if (i==2) rep(NA,2) else 1:3 )


# This function process all columns of x together instead of processing them one at a time
# FUN can return an array or a list. It does not have to return a scalar. This
#    saves from having to redo grouping for every col that has to be returned
#    also this eliminates the necessity to process a column of x at a time
myaggregate = function (x, by, FUN, new.col.name="aggregate.value", ...) 
{
    if (!is.data.frame(x)) 
        x <- as.data.frame(x)
    if (!is.list(by)) 
        stop(paste(sQuote("by"), "must be a list"))
    if (is.null(names(by))) 
        names(by) <- paste("Group", seq(along = by), sep = ".")
    else {
        nam <- names(by)
        ind <- which(nchar(nam) == 0)
        names(by)[ind] <- paste("Group", ind, sep = ".")
    }
    
    #original
    #y <- lapply(x, tapply, by, FUN, ..., simplify = FALSE)
    #modified
    z=mytapply(x, by, FUN, ...)
    
    #original
    #if (any(sapply(unlist(y, recursive = FALSE), length) > 1)) 
    #    stop(paste(sQuote("FUN"), "must always return a scalar"))
    #z <- y[[1]]
    
    d <- dim(z)
    w <- NULL
    for (i in seq(along = d)) {
        j <- rep.int(rep.int(seq(1:d[i]), prod(d[seq(length = i - 
            1)]) * rep.int(1, d[i])), prod(d[seq(from = i + 1, 
            length = length(d) - i)]))
        w <- cbind(w, dimnames(z)[[i]][j])
    }
    w <- w[which(!unlist(lapply(z, is.null))), , drop = FALSE]
        
    #original
    #y <- data.frame(w, lapply(y, unlist, use.names = FALSE))
    #modified
    stopifnot(length(unlist(z))%%nrow(w)==0)
    y <- data.frame(w, matrix(unlist(z), nrow=nrow(w), byrow=TRUE), stringsAsFactors=FALSE) # note stringsAsFactors option is hard coded here
    #original
    #names(y) <- c(names(by), names(x))
    #modified
    names(y) <- c(names(by), new.col.name)
    y
}

# This function can handle X being matrix instead of just a vector
mytapply = function (X, INDEX, FUN = NULL, ..., simplify = TRUE) 
{
    FUN <- if (!is.null(FUN)) 
        match.fun(FUN)
    if (!is.list(INDEX)) 
        INDEX <- list(INDEX)
    nI <- length(INDEX)
    namelist <- vector("list", nI)
    names(namelist) <- names(INDEX)
    extent <- integer(nI)
    
    #original
    #nx <- length(X)
    #modified
    nx = ifelse(!(is.data.frame(X) | is.matrix(X)), length(X), length(X[,1]) )
    
    one <- as.integer(1)   
    group <- rep.int(one, nx)
    ngroup <- one
    for (i in seq(INDEX)) {
        index <- as.factor(INDEX[[i]])
        if (length(index) != nx) 
            stop("arguments must have same length")
        namelist[[i]] <- levels(index)
        extent[i] <- nlevels(index)
        group <- group + ngroup * (as.integer(index) - one)
        ngroup <- ngroup * nlevels(index)
    }
    if (is.null(FUN)) 
        return(group)
    ans <- lapply(split(X, group), FUN, ...)
    index <- as.numeric(names(ans))
    if (simplify && all(unlist(lapply(ans, length)) == 1)) {
        ansmat <- array(dim = extent, dimnames = namelist)
        ans <- unlist(ans, recursive = FALSE)
    }
    else {
        ansmat <- array(vector("list", prod(extent)), dim = extent, 
            dimnames = namelist)
    }
    names(ans) <- NULL
    ansmat[index] <- ans
    ansmat
}

# idvar used to be optional: it was taken to be idvar=setdiff(names(dat), c(category.var,outcome.var))
myreshapewide=function(formula, dat, idvar, keep.extra.col=FALSE){
    tmp = as.character(formula)
    outcome.var=tmp[2]
    category.var=tmp[3]
    stopifnot(category.var %in% names(dat))
    
#    if (is.null(idvar)) {
#        idvar=setdiff(names(dat), c(category.var,outcome.var))
#        # if any of idvar has NA then it is a problem if not treated, here we opt to remove such columns
#        idvar=idvar[sapply(idvar, function (idvar.) all(!is.na(dat[,idvar.])) )]
#    } else {
    
    # remove those columns with changing values within an id
    # need [,-(1:length(idvar)),drop=FALSE] because the first columns are idvar
    if(keep.extra.col) {
        tmp=apply(aggregate (x=dat[,!names(dat) %in% c(idvar,category.var,outcome.var),drop=FALSE], by=dat[,names(dat) %in% idvar,drop=FALSE], function(y) {
            if(is.factor(y)) y=as.character(y)
            length(rle(y[!is.na(y)])$values)<2
        })[,-(1:length(idvar)),drop=FALSE], 
            2, all)
        varying.var=names(tmp)[!tmp]
        dat=dat[,!names(dat) %in% setdiff(varying.var, c(category.var,outcome.var)),drop=FALSE]
    } else {
        dat=dat[,names(dat) %in% c(category.var,outcome.var,idvar)]
    }
    
#    }
    reshape(dat, direction="wide", v.names=outcome.var, timevar=category.var, idvar=idvar )
}


myreshapelong=function(dat, cols.to.be.stacked, label.cols.to.be.stacked, new.col.name){
    reshape(dat, direction="long", varying = cols.to.be.stacked, v.names = new.col.name, timevar=label.cols.to.be.stacked, times=cols.to.be.stacked)
}


# keep column names
# return data frame if input is data frame
myscale=function(x){
    flag=is.data.frame(x)
    nam=dimnames(x)[[2]]
    aux=scale(x)
    dimnames(aux)[[2]]=nam
    if(flag) {
        tmp=as.data.frame(aux)
        mostattributes(tmp)=attributes(x)
        attr(tmp,"scaled:center")=attr(aux,"scaled:center")
        attr(tmp,"scaled:scale")=attr(aux,"scaled:scale")
#        
#        # add old attr
#        tmp1=attributes(aux)
#        for(a in tmp1){
#            #print(a)
#            print(names(a))
#            print("a")
#            #attr(tmp,"scaled:scale")=attr(x,"scaled:scale")
#        }
        aux=tmp
        names(aux)=names(x)
    }
    aux
}

meanmed=function(x, na.rm=FALSE){
    c(mean=mean(x, na.rm=na.rm), sd=sd(x, na.rm=na.rm), median=median(x, na.rm=na.rm), iqr=IQR(x, na.rm=na.rm), var=var(x, na.rm=na.rm))
}
read.tsv=function (file, header = TRUE, sep = "\t", ...) {
    read.csv(file, header = header, sep = sep, ...)
}

read.sv=function (file, header = TRUE, ...) {
    sep=","
    if (getExt(file)=="tsv") sep="\t"
    read.csv(file, header = header, sep = sep, ...)
}

keepWarnings <- function(expr) { 
    localWarnings <- list() 
    value <- withCallingHandlers(expr, 
    warning = function(w) { 
    localWarnings[[length(localWarnings)+1]] <<- w 
    invokeRestart("muffleWarning") 
    }) 
    list(value=value, warnings=localWarnings) 
} 

# make table that shows both counts/frequency and proportions
# style 1: count only; 2: count + percentage; 3: percentage only
table.prop=function (x,y=NULL,digit=1,style=2,whole.table.add.to.1=FALSE,useNA="ifany",add.perc=FALSE, add.total.column=FALSE) {
    if (is.null(y)) {
        tbl=table(x,useNA=useNA)
        whole.table.add.to.1=TRUE  # to trick the computation of prop
    } else {
        tbl=table(x,y,useNA=useNA)
    }
    if (whole.table.add.to.1) prop = tbl/sum(tbl) else prop = apply (tbl, 2, prop.table)
    if (style==2) {
        res = tbl %.% " (" %.% round(prop*100,digit) %.%ifelse(add.perc,"%","")%.% ")"
    } else if (style==3) {
        res = prop*100 # no need to do formating here
    } else if (style==4) {
        res = tbl %.% " (" %.% round(prop*100,digit) %.% "%)"
    } else if (style==5) {
        res = round(prop*100,digit)%.%"%" # no need to do formating here
    } else res=tbl
    res=matrix(res, nrow=nrow(tbl))    
    dimnames(res)=dimnames(tbl)
    names(dimnames(res))=NULL
    
    if (is.null(y)) res=res[,1]
    
    # add a total column 
    if (!is.null(y) & add.total.column) {
        tab2=table.prop (x,y=NULL,digit=digit,style=style,whole.table.add.to.1=whole.table.add.to.1,useNA=useNA,add.perc=add.perc) 
        res=cbind(" "=tab2, res)
    }
    
    res
}

# case is vector of 0/1 and group is vector of multi-group indicators
# the second row is taken in the computation, so be careful about NA's
table.cases=function (case,group,include.all=TRUE,desc="cases") {
    tbl=table(case, group)
#    tab=binom::binom.confint(x=c(tbl[2,],if(include.all) sum(tbl[2,])), n=c(colSums(tbl),if(include.all) sum(tbl)), alpha=0.05, method=c("wilson"))[,-1] # remove method column
    tab=Hmisc::binconf(x=c(tbl[2,],if(include.all) sum(tbl[2,])), n=c(colSums(tbl),if(include.all) sum(tbl)), alpha=0.05, method=c("wilson"), include.x=TRUE, include.n=TRUE)
    tab[,3:5]=tab[,3:5]*100 # get percentage
    colnames(tab)=c("# "%.%desc,"# total","% "%.%desc,"lb","ub")
    rownames(tab)=c(colnames(tbl),if(include.all) "all")
    tab
}

# case is vector of 0/1, group1 and group2 are vectors of multi-group indicators
# output:
#             low medium  high
#low 48% (12/25) 5% (1/19)   0% (0/6)
#medium  20% (3/15)  21% (4/19)  0% (0/16)
#high    20% (2/10)  25% (3/12)  0% (0/28)
# e.g. table.cases.3(HIVwk28preunbl, cut2(dat.merge[[a]],g=3), cut2(dat.merge[["cd8.env.poly"]],g=3))
table.cases.3=function(case,group1,group2){
    n.x=length(table(group1))
    n.y=length(table(group2))
    tab=table(group1, group2, case)
    #print(tab);str(tab);print(rownames(tab));print(colnames(tab))
    res=matrix(round(100*tab[,,2]/(tab[,,1]+tab[,,2])) %.%  "% (" %.% tab[,,2] %.% "/" %.%(tab[,,1]+tab[,,2])%.% ")",nrow=n.x)
    if (!is.null(rownames(tab))) rownames(res)=rownames(tab)# if(n.x==2) c("low","high") else c("low","medium","high")
    if (!is.null(colnames(tab))) colnames(res)=colnames(tab)# if(n.y==2) c("low","high") else c("low","medium","high")
    res
}

# from Thomas, found on R user group
methods4<-function(classes, super=FALSE, ANY=FALSE){ 
    if (super) classes<-unlist(sapply(classes, function(cl) getAllSuperClasses(getClass(cl)))) 
    if (ANY) classes<-c(classes,"ANY") 
    gens<-getGenerics()@.Data 
    sigs<-lapply(gens, function(g) linearizeMlist(getMethods(g))@classes) 
    names(sigs)<-gens@.Data 
    sigs<-lapply(sigs, function(gen){ gen[unlist(sapply(gen, function(sig) any(sig %in% classes)))]}) 
    sigs[sapply(sigs,length)>0] 
} 
 

# p1 and p2 are two points. return y that corresponds to x on the line between p1 and p2
interpolate=function(pt1, pt2, x){
    slope=(pt2-pt1)[2]/(pt2-pt1)[1]
    intercept = pt1[2]-slope*pt1[1]
    intercept + slope * x    
}

# x is a data frame or matrix. correlation between columns and p values are returned.
mycor=function(x, use = "everything", method = c("pearson", "kendall", "spearman"), 
    alternative = c("two.sided", "less", "greater"), exact = NULL, conf.level = 0.95, continuity = FALSE, 
    digits.coef=2, digits.pval=3,
    ...) {
    
    na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", "everything", "na.or.complete"))
    if (is.na(na.method)) stop("invalid 'use' argument")
    method <- match.arg(method)
    alternative <- match.arg(alternative)

    cor.coef=cor(x, y = NULL, use = use, method = method)
    
    cor.pval=diag(ncol(cor.coef))
    for (i in 1:ncol(x))
    for (j in 2:ncol(x)) {
        cor.pval[i,j]<-cor.pval[j,i]<-suppressWarnings(cor.test(x[,i], x[,j], alternative = alternative, method = method, exact = exact, conf.level = conf.level, continuity = continuity, ...)$p.value)
    }
    
    name=if (is.data.frame(x)) names(x) else colnames(x)
    out=matrix(NA,nrow=nrow(cor.coef),ncol=ncol(cor.coef), dimnames=list(name,name))
    for (i in 2:nrow(cor.coef))
    for (j in 1:(i-1)) {
        out[i,j]=formatDouble(cor.coef[i,j], digits.coef) %.% " (" %.% formatDouble(cor.pval[i,j], digits.pval) %.% ")" %.%
            ifelse (round(cor.pval[i,j],2)<=0.05,ifelse (cor.pval[i,j]<0.01,  ifelse (cor.pval[i,j]<0.01,"***","**")  ,"*"),"")
    }
    #out[-1,-ncol(out)]
    out
}
