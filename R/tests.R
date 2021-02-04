signtest = function (x) {
    binom.test( sum( x > 0 , na.rm=TRUE) , sum( x != 0 , na.rm=TRUE), p = 0.5,
           alternative = c("two.sided", "less", "greater"),
           conf.level = 0.95)   
}


quick.t.test = function (x, y, var.equal=FALSE) {
    mean.x = mean(x)
    mean.y = mean(y)
    m=length(x)
    n = length(y)
    if (var.equal) {
        (mean.x-mean.y)/sqrt( ( sum((x-mean.x)**2) + sum((y-mean.y)**2) ) * (1/m + 1/n) /(m+n-2) )
    }    else {
        (mean.x-mean.y)/sqrt( sum((x-mean.x)**2)/(m-1)/m + sum((y-mean.y)**2)/(n-1)/n )         
    }
}

vector.t.test = function (mean.x, mean.y, var.x, var.y, n) {
    new.var = (var.x + var.y) /n
    t.stat = abs(mean.x-mean.y)/sqrt(new.var)
    names(t.stat)=NULL
    t.stat
}

tukey.mtest = function (mu, ms, n) {
    #mu = c(45, 58, 46, 45, 56 )
    #ms = 5
    #n=3
    m=length(mu)
    cutoff = qtukey(p=.01, nmeans=m, df=(n-1)*(m-1), nranges = 1, lower.tail = FALSE, log.p = FALSE)/sqrt(2)
    
    t=matrix(0, m,m)
    for (i in 1:m) {
        for (j in i:m) {
            t[i,j] =  ( (mu[i]-mu[j])/ sqrt( 2/n * ms ) )
        }
    }
    
    cat ("The t statistics for Tukey method are calculated below:\n")
    print (signif(t,3))
    cat ("\n")
    
    cat ("By comparing the t values with ", signif (cutoff,3), ", Tukey method declares that the following pairs are significantly different: ")
    for (i in 1:m) {
        for (j in i:m) {
            if (abs(t[i,j]) > cutoff) {
                if (t[i,j]>0) cat (i, "&", j)
                else cat (j, "&", i)
                if (i==m & j==n) cat (". ")    
                else cat (", ")   
            }
        }
    }
    cat ("In other words, the following pairs are not significantly different: ")
    for (i in 1:m) {
        for (j in i:m) {
            if (abs(t[i,j]) <= cutoff & i!=j) {
                cat (i, "&", j)
                if (i==m & j==n) cat (". ")    
                else cat (", ")   
            }
        }
    }
    cat("\n")
    
}


#Hosmer-Lemeshow goodness of fit test
hosmerlem = function(y, yhat, g=10) {
  cutyhat = cut(yhat,
     breaks = quantile(yhat, probs=seq(0,
       1, 1/g)), include.lowest=TRUE)
  obs = xtabs(cbind(1 - y, y) ~ cutyhat)
  expect = xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
  chisq = sum((obs - expect)^2/expect)
  P = 1 - pchisq(chisq, g - 2)
  return(list(chisq=chisq,p.value=P))
}

# x is a matrix or dataframe
mycor.test=function(x, method = c("pearson", "kendall", "spearman"), idx=NULL) {
    p=ncol(x)
    out=matrix(NA,p,p)
    if (is.null(idx)) {
        for (j in 1:(p-1)) 
            for (i in (j+1):p) 
                out[i,j] <- out[j,i] <- suppressWarnings(cor.test(x[,i], x[,j], method=method))$p.value
                
    } else {
        stopifnot(length(idx)==2)
        for (j in idx[[2]]) 
            for (i in idx[[1]]) 
                out[i,j] <- suppressWarnings(cor.test(x[,i], x[,j], method=method))$p.value
        
    }
    out
}

# pop1 is the number of cases and controls in the first population
# pop2 is the number of cases and controls in the second population
mybarnard.test=function(pop1, pop2, method="z-pooled"){
    # each row of tab is a population
    tab=rbind(pop1, pop2)
    test=Exact::exact.test(tab, npNumbers = 100, method = method, to.plot=F)
    test    
}
