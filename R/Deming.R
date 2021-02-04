# modified from MethComp
Deming <- function( x, y, vr=sdr^2, sdr=sqrt(vr), boot=TRUE, keep.boot=FALSE, alpha=0.05 ) {

    if( missing( vr ) & missing( sdr ) ) var.ratio <- 1
    else var.ratio <- vr
    vn <- c( deparse( substitute( x ) ),
             deparse( substitute( y ) ) )
    pn <- c( "Intercept", "Slope", paste( "sigma", vn, sep="." ) )
    
    alfa <- alpha
    dfr <- data.frame( x=x, y=y )
    dfr <- dfr[complete.cases(dfr),]
    x <- dfr$x
    y <- dfr$y
    n <- nrow( dfr )
    SSDy <- var( y )*(n-1)
    SSDx <- var( x )*(n-1)
    SPDxy <- cov( x, y )*(n-1)
    beta <- ( SSDy - var.ratio*SSDx +
              sqrt( ( SSDy - var.ratio*SSDx )^2 +
                    4*var.ratio*SPDxy^2 ) ) / ( 2*SPDxy)
    alpha <- mean( y ) - mean( x ) * beta
    ksi <- ( var.ratio*x + beta*(y-alpha) )/(var.ratio+beta^2)
    sigma.x <- ( var.ratio*sum( (x-ksi)^2 ) +
                           sum( (y-alpha-beta*ksi)^2 ) ) /
    # The ML-estiamtes requires 2*n at this point bu we do not think we have that
    # many observations so we stick to (n-2). Any corroboation from litterature?
    
    # YF: See Seber and Wild (2003) page 495 for justification 
               ( (n-2)*var.ratio )
    sigma.y <- var.ratio*sigma.x
    sigma.x <- sqrt( sigma.x )
    sigma.y <- sqrt( sigma.y )
    res <- c( alpha, beta, sigma.x, sigma.y )
    names( res ) <- pn
    out=list(coef=res)
    if(boot) {
        if( is.numeric( boot ) ) N <- boot else N <- 1000
        boot.samples <- matrix( NA, N, 4 )
        colnames( boot.samples ) <- pn
        for( i in 1:N ) {
           wh <- sample( 1:n, n, replace=TRUE )
           boot.samples[i,] <- Deming( x[wh], y[wh], vr=var.ratio, boot=FALSE )$coef
        } 
        vcov=cov( boot.samples )
        dimnames(vcov)=list(pn, pn)
        out$vcov=vcov
        alfa=0.05
        boot.percentiles= t( apply( boot.samples, 2, quantile, probs=c(0.5,alfa/2,1-alfa/2 ), na.rm=T ) ) 
        out$boot.percentiles=boot.percentiles
        if(keep.boot) out$boot.samples=boot.samples
    }
    
    class(out)="Deming"
    out
}

summary.Deming=function  (object, ...) {
        getFixedEf(object)
}
coef.Deming = function  (object, ...) {
    object$coef
}
vcov.Deming = function  (object, robust, ...) {
    object$vcov
}
getFixedEf.Deming = function (object, ...) {
    res <- cbind( object$coef,
                   se <- sqrt( diag(object$vcov) ),
                   object$boot.percentiles[,2:3])
    colnames( res )<- c("Estimate", "S.e.(boot)", colnames(res)[4:5] )
    res = cbind(res, "p value"= 2*pnorm(abs(res[,1])/res[,2], lower.tail=FALSE) )
    colnames(res)=c("est", "se(est)", "(lower", "upper)", "p.value")
    res
}
predict.Deming = function (object, newdata = NULL, ...) {
    out=object$coef["Intercept"] + object$coef["Slope"]*c(newdata)
    out
}
