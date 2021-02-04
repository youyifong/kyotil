

binaryloess <- function(x, y, scale=c("logit","linear"), span=0.7, weights=NULL, ...) {
    
    if (missing(scale)) stop("scale mssing")
    scale=match.arg(scale)
    
    complete=!is.na(x) & !is.na(y)
    x=x[complete]
    y=y[complete]    
    if (is.null(weights)) weights=rep(1, length(x)) else weights=weights[complete]
    
    loessfit <- predict(loess(y~x,span=span, weights=weights))
    # truncate
    pi <- pmax(pmin(loessfit,0.9999),0.0001)
    
    if(scale=="logit") {
        yy <- logit(pi)    
    } else if (scale=="linear") {
        yy =  pi    
    }
    
    plot(x, yy, ...)
    #lines (lowess(y~x,f=span))    
}
