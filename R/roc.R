computeRoc<-function(score, outcome, reverse.sign.if.nece=TRUE, cutpoints=NULL){
    # removing na
    outcome=outcome[!is.na(score)]
    score=score[!is.na(score)]
    
    # reverse sign if necessary
    if (reverse.sign.if.nece) {
        if(fastauc(score, outcome, FALSE)<.5) score=-score
    }
    
    if(is.null(cutpoints)) cutpoints<-c(-Inf,sort(unique(score)))
    sensitivity<-sapply(cutpoints, function(cutp) mean(score>cutp & outcome)/mean(outcome))
    specificity<-sapply(cutpoints, function(cutp) mean(score<=cutp & !outcome)/mean(!outcome))
    list(sensitivity=sensitivity, specificity=specificity, threshold=cutpoints)
}



plotRoc<-function(x,add=FALSE,type="l",diag.line=TRUE,...){
    if(add) {
        lines(1-x$specificity, x$sensitivity,xlab="1-specificity",ylab="sensitivity", type=type, ...) 
    } else {
        plot(1-x$specificity, x$sensitivity,xlab="1-specificity",ylab="sensitivity", type=type, ...)
        if (diag.line) abline(0,1, col="gray",lwd=.5)
    }
}

addRoc<-function(x,...){
    lines(1-x$specificity,x$sensitivity,...)
}


classification.error=function(score, outcome, threshold=NULL, verbose=FALSE){
    
   if (all(is.na(score)) | all(is.na(outcome))) {
        res=NA
        warning("all score NA or all outcome NA")
    } else {
        if(is.null(threshold)) threshold=quantile(score, 1-mean(outcome), na.rm=TRUE)
        if (verbose) myprint(threshold)        
        pred=as.numeric(score>threshold)
        res=mean(pred!=outcome)
    }
    
    res
}
