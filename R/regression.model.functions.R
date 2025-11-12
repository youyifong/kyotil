formatPvalues=function(pvalues, p.digits=3) {
  ifelse(pvalues<10^(-p.digits), paste0("<0.",concatList(rep("0",p.digits-1)),"1"), 
         format(round(pvalues, p.digits), nsmall=p.digits, scientific=FALSE) )
}

# type 3 and 7 do not give the right output for glm fits
# robust can be passed as part of .... Sometimes robust=T generates errors
# trim: get rid of white space in confidence intervals for alignment
getFormattedSummary=function(fits, type=12, est.digits=2, se.digits=2, robust, random=FALSE, VE=FALSE, 
                             to.trim=FALSE, rows=NULL, coef.direct=FALSE, trunc.large.est=TRUE, 
                             scale.factor=1, p.digits=3, remove.leading0=FALSE, p.adj.method="fdr", sig.level=0.05, ...){
# type=12; est.digits=2; se.digits=2; robust; random=FALSE; VE=FALSE; to.trim=FALSE; rows=NULL; coef.direct=FALSE; trunc.large.est=TRUE; scale.factor=1; p.digits=3; remove.leading0=FALSE; p.adj.method="fdr"    
    if(is.null(names(fits))) names(fits)=seq_along(fits)
    idxes=seq_along(fits); names(idxes)=names(fits)
    
    if (type==11) {
      if (is.null(rows)) stop("for type==11, the rows argument needs to be provided")
      if(length(rows)!=1) stop("type==11 requires only one coef selected")     
    } else if (type==13) {
      if (is.null(rows)) stop("for type==13, the rows argument needs to be provided")
    }
    
    if(length(robust)==1) {
      robust=as.list(rep(robust,length(fits))) 
    } else if (length(robust)!=length(fits)) {
      stop("length of robust needs to match length of fits")    
    }
    
    stopifnot(is.list(robust))
    
    res = sapply(idxes, simplify="array", function (fit.idx) {
        
        # "est","se","(lower","upper)","p.value"
        
        fit=fits[[fit.idx]]
        if(coef.direct) {
            tmp=fit
        } else {
            if (random) {
                tmp = getVarComponent (fit)
                if (inherits(fit, "mer") & type!=1) {
                    warning ("only point estimate is available for variance components from lmer fit, forcing type to 1")
                    type=1
                }
            } else {
                tmp = getFixedEf (fit, robust=robust[[fit.idx]], scale.factor=scale.factor, ...)
            }
            
            if (VE) {
                tmp[,1]=1-tmp[,1]
                # reverse lb and ub
                tmpv=tmp[,4]
                tmp[,4]=1-tmp[,3]
                tmp[,3]=1-tmpv
            }
        }
        
        
        # take only some rows
        if (!is.null(rows)) tmp=tmp[intersect(rows,1:nrow(tmp)),,drop=F]
    
        est.=as.matrix(tmp[,1,drop=FALSE], ncol=1) # if do not cast into matrix, it will be a data frame, which leads to problems next
        too.big=trunc.large.est & est.>100
        
        # replace large values in est.
        # find a round that takes multipled digits!
        # trim is necessary here
        est.=ifelse(too.big,">100",trim(formatDouble(est., est.digits, remove.leading0=remove.leading0)))
        
        p.val.col=which(startsWith(tolower(colnames(tmp)),"p"))
        
        if(ncol(tmp)>1) {
            lb=tmp[,3,drop=FALSE]
            lb[too.big]=NA
            lb=ifelse(trunc.large.est & lb>100, ">100", formatDouble(lb, est.digits, remove.leading0=remove.leading0) )
            if(to.trim) lb=trim(lb)
            
            ub=tmp[,4,drop=FALSE]
            ub[too.big]=NA
            ub=ifelse(trunc.large.est & ub>100, ">100", formatDouble(ub, est.digits, remove.leading0=remove.leading0) )
            if(to.trim) ub=trim(ub)
        }
                
        # str(lb); str(ub); str(est.) # make sure they are not data frames 
        
        if (type==0)
            # return tmp
            out=drop(tmp)
        else if (type==1)
            # est only
            out=drop(est. )
        else if (type==2)
            # est (se)
            out=est. %.% " (" %.% formatDouble(tmp[,2,drop=FALSE], se.digits, remove.leading0=remove.leading0) %.% ")" %.% ifelse (round(tmp[,p.val.col],3)<=sig.level, ifelse (tmp[,p.val.col]<0.01,"**","*"),"")
        else if (type==3) 
            # est (lb, ub)
            out=est. %.% " (" %.% lb %.% "," %.% ub %.% ")" 
        else if (type==7)
          # (lb, ub)
          out=ifelse(drop(too.big), rep("",nrow(lb)), " (" %.% lb %.% ", " %.% ub %.% ")")
        else if (type==4)
            # a space is inserted between est and se, they could be post-processed in swp
            out=est. %.% " " %.% formatDouble(tmp[,2,drop=FALSE], est.digits, remove.leading0=remove.leading0)
        else if (type==5)
            # est **
            out=est. %.%
                ifelse (round(tmp[,p.val.col],3)<=sig.level,ifelse (tmp[,p.val.col]<0.01,"**","*"),"")
        else if (type==6)
            # est (pval)*
            out=est. %.% " (" %.% formatDouble(tmp[,p.val.col,drop=FALSE], 3, remove.leading0=remove.leading0) %.% ")" %.% ifelse (round(tmp[,p.val.col],3)<=sig.level,ifelse (tmp[,p.val.col]<0.01,"**","*"),"")
        else if (type==8)
            # est (p value #)
            out=est. %.% " (p value " %.% 
                formatDouble(tmp[,p.val.col,drop=FALSE], 3, remove.leading0=remove.leading0) %.% ")" 
        else if (type==9)
            # est (pval)
            out=est. %.% " (" %.% formatDouble(tmp[,p.val.col,drop=FALSE], p.digits, remove.leading0=remove.leading0) %.% ")" 
        else if (type==10) {
            # pval
            tmp.out=tmp[,p.val.col,drop=TRUE]
            # e.g., transform 0.000 to <0.001
            out = ifelse(tmp.out<10^(-p.digits), paste0("<0.",concatList(rep("0",p.digits-1)),"1"), 
                         format(round(tmp.out, p.digits), nsmall=p.digits, scientific=FALSE) )
        }
        else if (type==11) {
            # adj pval will be done later because it needs info from all models
            out=tmp[,p.val.col]

        } else if (type==12) {
            # est (lb, up, pval *)
            out=est. %.% 
                " (CI=" %.% lb %.% "," %.% ub %.% 
                ", p=" %.% formatDouble(tmp[,p.val.col,drop=FALSE], 3, remove.leading0=remove.leading0) %.% ")" %.%
                ifelse (round(tmp[,p.val.col],3)<=sig.level,ifelse (tmp[,p.val.col]<0.01,"**","*"),"") 
  
        } else if (type==13) {
            # generalized Wald test p value
            fit=fits[[fit.idx]]
            tmp.out = if (length(fit)==1) NA else {
              stat=coef(fit)[rows] %*% solve(vcov(fit,robust=robust)[rows,rows]) %*% coef(fit)[rows]
              pchisq(stat, length(rows), lower.tail = FALSE)
            }
            out = ifelse(tmp.out<10^(-p.digits), paste0("<0.",concatList(rep("0",p.digits-1)),"1"), 
                         format(round(tmp.out, p.digits), nsmall=p.digits, scientific=FALSE) )
            
        
        } else 
            stop ("getFormattedSummary(). type not supported: "%.%type)
        
        # 13: a single output from multiple terms
        if(!type %in% c(0,13)) {
            names(out)=rownames(tmp)
            out=gsub("NA","",out)
            out=gsub("\\( +\\)","",out)
        }
        
        out
    })
        
    if (is.list(res)) {
    # if the fits have different number of parameters, we need this
        res=cbinduneven(li=lapply(res, function (x) as.matrix(x, ncol=1)))
        colnames(res)=names(fits)
    } else if (!is.matrix(res)) {
    # if there is only one coefficient, we need this
        res=matrix(res, nrow=1)
        colnames(res)=names(fits)
    }
    
    # the following processing requires info across fits
    if(type==11) {
      res=format(round(p.adjust(res, method=p.adj.method), p.digits), nsmall=3, scientific=FALSE) 
      names(res)=names(fits)
    }
    
    res
}


# return a matrix, columns are: est, p value, low bound, upper bound
getFixedEf <- function(object, ...) UseMethod("getFixedEf") 
# variance component
getVarComponent <- function(object, ...) UseMethod("getVarComponent")

## need to export in namespace
# getFixedEf.try-error <- function (object, ...) NA
  

# if there is missing data and robust=T, some errors will happen. it is a useful error to have.
# if ret.robcov TRUE, then returns robust variance-covariance matrix
# infjack.glm defined later in this file
getFixedEf.glm = function (object, exp=FALSE, robust=TRUE, ret.robcov=FALSE, scale.factor=1, ...) {
  
    x=summary(object)
    out=x$coef
    if (robust | ret.robcov) {
        V=infjack.glm(object, 1:nrow(object$model)) # object$data may have NAs
        if (ret.robcov) return (V)
        out[,2]=sqrt(diag(V))
        out[,3]=out[,1]/out[,2]
        out[,4]=2 * pnorm(-abs(out[,3]))
    }
    
    out[,1]=out[,1]*scale.factor
    out[,2]=out[,2]*scale.factor
    
    # to get est, p value, low bound, upper bound
    out=cbind(out[,1],out[,2],out[,1]-1.96*out[,2],out[,1]+1.96*out[,2],out[,4])
    colnames(out)=c("est","se","(lower","upper)","p.value")
    if(exp) {
        out[,c(1,3,4)]=exp(out[,c(1,3,4)])
        colnames(out)[1]="OR"
    }
    
    # sometimes some coef are null, we may want to preserve, modified from print.summary.glm
    if(!is.null(aliased <- x$aliased) && any(aliased)) {
        cn <- names(aliased)
        coefs <- matrix(NA, length(aliased), ncol(out),
                        dimnames=list(cn, colnames(out)))
        coefs[!aliased, ] <- out
        out = coefs
    }
    
    out
}

# to use finite sample correction, object has to come from geepack::geeglm
# robust: T, F, or a string
getFixedEf.gee = function (object, robust=TRUE, exp=FALSE, scale.factor=1, ...) {
  
  # whether this is model-based or sandwich depends on how the model is fit
  out=as.matrix(summary(object)$coef)
  
  if (robust==FALSE) {
    # do nothing
    
  } else {
    # get robust se
    res = ssase(object)

    if (robust==TRUE) {
      se = res$NO
    } else {
      se = res[[robust]]
      if (is.null(se)) stop (robust %.% "not in available robust options MBN, MD, KC, mFG")
    }
    
    se = sqrt(diag(se))

    p.value = 2 * pnorm(abs(object$coefficients / se), lower.tail = FALSE)   # two-sided p-value
    
    out = cbind(object$coefficients, se, 1, p.value)
  }
  

  out[,1]=out[,1]*scale.factor
  out[,2]=out[,2]*scale.factor

  out=cbind(out[,1:2], 
            out[,1]-1.96*out[,2], 
            out[,1]+1.96*out[,2], 
            out[,4,drop=FALSE])
  colnames(out)=c("est", "se",
                  "(lower",
                  "upper)",
                  "p.value")
  
  if(exp) {
      out[,c(1,3,4)]=exp(out[,c(1,3,4)])
      colnames(out)[1]="RR"
  }
  out
}


getFixedEf.MIresult=function(object,exp=FALSE, scale.factor=1, ...) {
    capture.output({
        tmp=summary(object)
    }, type="output") # type = message captures stderr, type=output is for stdout
    
    tmp=tmp[,names(tmp)!="missInfo"] # MIresult has this string column  #tmp=subset(tmp, select=-missInfo) fails check
    tmp=as.matrix(tmp)# data frame fails in format and formatDouble
    
    out = cbind(tmp, "p.val"=2*pt(abs(tmp[,1]/tmp[,2]), df=object$df, lower.tail = FALSE))
    
    out[,1:4]=out[,1:4]*scale.factor

    if(exp) {
        out[,c(1,3,4)]=exp(out[,c(1,3,4)])
    }
    
    out
}

#get estimates, variances, sd from lmer fit
getFixedEf.mer = function (object, scale.factor=1, ...) {
    Vcov <- lme4::VarCorr(object, useScale = FALSE) 
    betas <- lme4::fixef(object) 
    se <- sqrt(diag(Vcov)) 
    zval <- betas / se 
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE) 
    
    out = cbind("Estimate"=betas, se, zval, pval) 
    
    out[,1:4]=out[,1:4]*scale.factor

    out
}
getVarComponent.mer = function (object, ...) {
    tmp=lme4::VarCorr(object)
    mysapply(tmp, function (comp) attr(comp, "stddev") )
}

#get estimates, variances, sd from lmer fit
getFixedEf.glmerMod = function (object, exp=FALSE, scale.factor=1, ...) {
    out=summary(object)$coefficients
    out=cbind(est=out[,1], se=out[,2], lb=out[,1]-1.96*out[,2], ub=out[,1]+1.96*out[,2],"p"=out[,4])
    
    out[,1:4]=out[,1:4]*scale.factor

    if(exp) out[,c(1,3,4)]=exp(out[,c(1,3,4)])
    out #est, se, lb, ub, pvalue 
}
getFixedEf.merMod=getFixedEf.glmerMod

getVarComponent.glmerMod = function (object, ...) {
    tmp=lme4::VarCorr(object)
    mysapply(tmp, function (comp) attr(comp, "stddev") )
}
getVarComponent.merMod=getVarComponent.glmerMod


#get estimates, variances, sd from lme fit
getFixedEf.lme = function (object, scale.factor=1, ...) {
    betas <- object$coef$fixed
    se <- sqrt (diag (object$varFix))
    zval <- betas / se 
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE) 
    out = cbind(betas, se, zval, pval) 
    
    out[,1:4]=out[,1:4]*scale.factor
    out
}
# getVarComponent.lme = function (object, ...) {
#     nlme::VarCorr(object)
# }

getFixedEf.lmerMod = function (object, exp=F, exp10=F, scale.factor=1, ...) {
    betas <- nlme::fixef(object)
    se <- sqrt (diag (getVarComponent(object)))
    zval <- betas / se 
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE) 
    out=cbind(betas, 
              se, 
              "(lower"=betas-qnorm(0.975)*se, 
              "upper)"=betas+qnorm(0.975)*se, 
              zval, 
              pval) 
    out[,1:4]=out[,1:4]*scale.factor
    if(exp) out[,c(1,3,4)]=exp(out[,c(1,3,4)])
    if(exp10) out[,c(1,3,4)]=10**(out[,c(1,3,4)])
    out
}
getVarComponent.lmerMod = function (object, ...) {
    as.matrix(vcov(object)) # otherwise will complain about S4 class convertibility problem
}

getFixedEf.geese = function (object, robust=TRUE, scale.factor=1, ...) {
    tmp=summary(object)$mean
    # tmp is: estimate     san.se         wald            p
    # tmp should be: est, se, lb, ub, pvalue 
    out = as.matrix(cbind(tmp[,c(1, ifelse(robust,2,3))], confint(object), tmp[,4,drop=F])) # if leave as data.frame, formatDouble fails
    
    out[,1:4]=out[,1:4]*scale.factor

    out
}
    
getFixedEf.logistf = function (object, exp=FALSE, scale.factor=1, ...) {
    temp = summary(object)
    out = cbind (coef=temp$coef, se=NA, p.value=temp$prob)
    out[,1]=out[,1]*scale.factor
    out[,2]=out[,2]*scale.factor
    if(exp) out[,1]=exp(out[,1])
    out
}
vcov.logistf=function(object, ...){
    object$var
}

getFixedEf.gam = function (object, scale.factor=1, ...) {
    temp = summary(object)
    out = cbind (temp$p.coef, temp$se[1:length(temp$p.coef)])
    out[,1]=out[,1]*scale.factor
    out[,2]=out[,2]*scale.factor
    out
}

getFixedEf.lm = function (object, exp=F, scale.factor=1, robust=FALSE, ...) {
  
  if (robust) stop("robust has not been implemented")
  
  out=summary(object)$coef
  ci=confint(object)    
  out=cbind(out[,1:2,drop=F], ci, out[,4,drop=F])
  colnames(out)[5]="p-val"
  out[,1:4]=out[,1:4]*scale.factor
  if(exp) {
      out[,c(1,3,4)]=exp(out[,c(1,3,4)])
  }
  out
}

getFixedEf.inla = function (object, ...) {
    tmp = summary(object)$fixed
    n=nrow(tmp)
    tmp.name = row.names(tmp)[n]
    # move intercept to the top
    if (tmp.name=="intercept") {
        tmp = rbind (tmp[n,],tmp)[1:n,,drop=FALSE]
        dimnames (tmp)[[1]][1] = tmp.name
    }
    # rename first column
    dimnames (tmp)[[2]][1] = "Estimate"
    tmp
}
# return the mean, sd, CI of the transformed variable
inla.getMeanSd=function (marginal, f="identity") {
    
    # interpolations suggested by Havard: do it on the original scale
    logtao=log(marginal[,1]); p.logtao=marginal[,2]*marginal[,1]
    fun = splinefun(logtao, log(p.logtao)) 
    h=0.001
    x = seq(min(logtao),max(logtao),by=h) 
    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
#    x = seq(min(logtao)-sd(logtao)/2,max(logtao)+sd(logtao)/2,by=h) 
#    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
#    x = seq(min(logtao)-sd(logtao),max(logtao)+sd(logtao),by=h) 
#    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
#    x = seq(min(logtao)-sd(logtao)*2,max(logtao)+sd(logtao)*2,by=h) 
#    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
    
    lower.boundary = rle(cumsum(pmf)>.025)$lengths[1]+1
    upper.boundary = rle(cumsum(pmf)>.975)$lengths[1]+1
#    if (pmf[lower.boundary]>.04) {
#        #stop ("getMeanSd(): pmf too large at lower boundary: "%.%pmf[lower.boundary])
#        return (rep(NA, 4))
#    }
#    if (pmf[upper.boundary]>.04) {
#        stop ("getMeanSd(): pmf too large at upper boundary"%.%pmf[upper.boundary])
#        return (rep(NA, 4))
#    }
    
    if (f=="identity") {
        func=function(x) { exp(x) }
    } else if (f=="inverse") {
        func=function(x) { exp(x)**-1 }
    } else if (f=="inversesqrt") {
        func=function(x) { exp(x)**-.5 }
    } else 
        stop ("getMeanSd(): function not supported "%.%f)
    
    mu = sum( pmf * func(x))
    stdev = (sum( pmf * func(x)**2 ) - mu**2) ** .5
    out = c("mean"=mu, "stddev"=stdev, 0,0 ) # we may run into problems with the boundary
    #out = c("mean"=mu, "stddev"=stdev, sort(func(x)[c(lower.boundary, upper.boundary)]) ); 
    names(out)[3]="2.5%"; names(out)[4]="97.5%";
    out
}
# returns estimate of standard deviation and the estimated sd of that estimate
getVarComponent.hyperpar.inla = function (object, transformation=NULL, ...) {
    marginals = object$marginals
    out = mysapply(1:length(marginals), function (i) {  
        # this is a little precarious, but hey
        if (startsWith(names(marginals)[i],"Prec")) {
            if (is.null (transformation)) {      
                inla.getMeanSd(marginals[[i]],"inversesqrt")
            } else {
                inla.getMeanSd(marginals[[i]],transformation)
            }
        } else if (startsWith(names(marginals)[i],"Rho")) {
            object$summary[i, c(1,2,3,5)]
        } else {
            stop ("don't know what to do with this names(marginals)[i]: "%.% names(marginals)[i] )
        }
    })
    dimnames (out)[[1]]="sigma.or.rho."%.%dimnames (out)[[1]]
    out
}

getFixedEf.svycoxph=function (object, exp=FALSE, robust=TRUE, scale.factor=1, ...){
    if (!robust) warning("svycoxph variance estimate is always design-based, which is close to robust variance estimate. No model-based variance estimate is available")
    hr=object$coefficients
    se=sqrt(diag(object$var))
    pval=pnorm(abs(hr/se), lower.tail=F)*2
    
    hr=hr*scale.factor
    se=se*scale.factor
    
    out=cbind(HR=hr,
            "se"=se,
            "(lower"=hr-qnorm(0.975)*se, 
            "upper)"=hr+qnorm(0.975)*se, 
            "p.value"=pval
        )
    if (exp) out[,c(1,3,4)]=exp(out[,c(1,3,4)])
    out
#    round(sqrt(diag(attr(object$var,"phases")$phase1)),3)
#    round(sqrt(diag(attr(object$var,"phases")$phase2)),3)    
}


getFixedEf.svyglm=function (object, exp=FALSE, robust=TRUE, scale.factor=1, ...){
    if (!robust) warning("svyglm variance estimate is always design-based, which is close to robust variance estimate. No model-based variance estimate is available")
    tmp=summary(object)$coefficients
    or=tmp[,"Estimate"]
    se=tmp[,"Std. Error"]
    pval=tmp[,"Pr(>|t|)"]
    out=cbind(OR=or,
            "se"=se,
            "(lower"=or-qnorm(0.975)*se, 
            "upper)"=or+qnorm(0.975)*se, 
            "p.value"=pval
        )
    out[,1:4]=out[,1:4]*scale.factor
    if (exp) out[,c(1,3,4)]=exp(out[,c(1,3,4)])
    out
#    round(sqrt(diag(attr(object$var,"phases")$phase1)),3)
#    round(sqrt(diag(attr(object$var,"phases")$phase2)),3)    
}


getFixedEf.svy_vglm=function (object, exp=FALSE, robust=TRUE, scale.factor=1, ...){
    if (!robust) warning("svy_vglm variance estimate is always design-based, which is close to robust variance estimate. No model-based variance estimate is available")
    tmp=summary(object)$coeftable
    or=tmp[,"Coef"]
    se=tmp[,"SE"]
    pval=tmp[,"p"]
    out=cbind(HR=or,
            "se"=se,
            "(lower"=or-qnorm(0.975)*se, 
            "upper)"=or+qnorm(0.975)*se, 
            "p.value"=pval
        )
    out[,1:4]=out[,1:4]*scale.factor
    if (exp) out[,c(1,3,4)]=exp(out[,c(1,3,4)])
    out
}


getFixedEf.coxph=function (object, exp=FALSE, robust=FALSE, scale.factor=1, ...){
    #capture.output()# summary.svycoxph prints some stuff, capture.output absorbs it
    sum.fit<-summary(object)
    se.idx=ifelse(!robust,"se(coef)","robust se")
    if(robust) {
        pvals=sum.fit$coef[,"Pr(>|z|)"]
    } else {
        pvals=pnorm(abs(sum.fit$coef[,1]/sum.fit$coef[,"se(coef)"]), lower.tail=F)*2
    }
    
    out=cbind(HR=sum.fit$coef[,1], 
              "se"=sum.fit$coef[,se.idx],
              "(lower"=(sum.fit$coef[,1]-qnorm(0.975)*sum.fit$coef[,se.idx]), 
              "upper)"=(sum.fit$coef[,1]+qnorm(0.975)*sum.fit$coef[,se.idx]), 
              "p.value"=pvals
    )
    out[,1:4]=out[,1:4]*scale.factor
    if (exp) out[,c(1,3,4)]=exp(out[,c(1,3,4)])
    out
}

getFixedEf.matrix = function (object, ...) {
    t(apply(object, 2, function (x) c("Estimate"=mean(x), "sd"=sd(x), "2.5%"=quantile(x,.025), "97.5"=quantile(x,.975))))
}
getVarComponent.matrix = function (object, ...) {
    t(apply(object, 2, function (x) c("Estimate"=mean(x), "sd"=sd(x), "2.5%"=quantile(x,.025), "97.5"=quantile(x,.975))))
}


#################################################################################################

coef.geese = function  (object, ...) {
    tmp = summary(object)$mean[,1]
    names (tmp)=names (object$beta)
    tmp
}
vcov.geese = function  (object, ...) {
    tmp = object$vbeta
    dimnames (tmp)=list (names (object$beta), names (object$beta))
    tmp
}
residuals.geese = function (object, y, x, ...) {
    y - x %*% object$beta   
}
predict.geese = function (object, x, ...) {
    x %*% object$beta 
}



#################################################################################################

coef.tps = function  (object, ...) {
    object$coef
}
vcov.tps = function  (object, robust, ...) {
    if(robust) object$cove else object$covm
}
getFixedEf.tps = function (object, exp=FALSE, robust=TRUE, scale.factor=1, ...) {
    res = summary(object)$coef
    idx=ifelse(robust,"Emp ","Mod ")
    res = cbind(res[,c("Value",idx%.%"SE")], "lower bound"=res[,"Value"]-1.96*res[,idx%.%"SE"], "upper bound"=res[,1]+1.96*res[,idx%.%"SE"], "p value"=res[,idx%.%"p"])
    res[,1:4]=res[,1:4]*scale.factor
    
    if (exp) res[,c(1,3,4)] = exp(res[,c(1,3,4)])
    colnames(res)=c(ifelse(exp,"OR","est"), "se(est)", "(lower", "upper)", "p.value")
    res
}
predict.tps = function (object, newdata = NULL, type = c("link", "response"), ...) {
    type=match.arg(type)
    out=as.matrix(newdata) %*% coef.tps(object) 
    if(type=="response") out=expit(out)
    out
}





# make two tables, rotating v1 and v2
# fit is an object that needs to have coef and vcov
# list the effect of one variable at three levels of another variable, if data is in fit, use quantiles, otherwise -1,0,1
# v1.type and v2.type: continuous, binary, etc
interaction.table=function(fit, v1, v2, v1.type="continuous", v2.type="continuous", logistic.regression=TRUE){
    
    coef.=coef(fit) 
    cov.=vcov(fit)    
    var.names=names(coef.)
    v1.ind = match(v1, var.names)
    v2.ind = match(v2, var.names)
    itxn.ind = match(v1%.%":"%.%v2, var.names)
    if(is.na(itxn.ind)) itxn.ind = match(v2%.%":"%.%v1, var.names)
    if (any(is.na(c(v1.ind, v2.ind, itxn.ind)))) {
        stop("v1, v2, or interaction not found in var.names")
    }
    
    ret=list()
    for (i in 1:2) {
        
        if (i==1) {
            ind.1=v1.ind; type.1=v1.type
            ind.2=v2.ind; type.2=v2.type
        } else {
            ind.1=v2.ind; type.1=v2.type
            ind.2=v1.ind; type.2=v1.type
        }
        
        # three levels of v.2
        if(is.null(fit$data)) {
            lev=c(-1,0,1)
            names(lev)=lev
        } else {
            # get quantiles
            lev=quantile(fit$data[[var.names[ind.2]]], c(.25,.5,.75), na.rm=TRUE)
        }
        
        # increase ind.1 by 1 at lev
        lin.combs=NULL
        if(type.2=="binary") {
            lin.comb.1=rep(0, length(coef.)); lin.comb.1[ind.1]=1; lin.comb.1[itxn.ind]=0;  lin.combs=rbind(lin.combs, "0"=lin.comb.1)
            lin.comb.1=rep(0, length(coef.)); lin.comb.1[ind.1]=1; lin.comb.1[itxn.ind]=1;  lin.combs=rbind(lin.combs, "1"=lin.comb.1)
        } else {
            lin.comb.1=rep(0, length(coef.)); lin.comb.1[ind.1]=1; lin.comb.1[itxn.ind]=lev[1]; lin.combs=rbind(lin.combs, lin.comb.1)
            lin.comb.1=rep(0, length(coef.)); lin.comb.1[ind.1]=1; lin.comb.1[itxn.ind]=lev[2]; lin.combs=rbind(lin.combs, lin.comb.1)
            lin.comb.1=rep(0, length(coef.)); lin.comb.1[ind.1]=1; lin.comb.1[itxn.ind]=lev[3]; lin.combs=rbind(lin.combs, lin.comb.1)             
            rownames(lin.combs)=names(lev)         
        }
    
        effect = lin.combs%*%coef.
        sd. = sqrt(diag(lin.combs%*%cov.%*%t(lin.combs)))
        p.val = pnorm(abs(effect)/sd., lower.tail=FALSE)
        lci=effect-1.96*sd.
        uci=effect+1.96*sd.
        
        res = cbind(effect, lci, uci, p.val)
        colnames(res)=c("coef","(lower","upper)","p value")
        if (logistic.regression) {
            res[,1:3]=exp(res[,1:3])
            colnames(res)[1]="OR"
        }
        
        ret[[i]]=res 
    }    
    
    names(ret) = "Effect of increasing "%.%c(v1,v2) %.% " by 1 at selected values of " %.% c(v2,v1)
    ret
       
}

# these functions are needed by some of the functions in this file
# returns sandwich estimator of variance matrix
# from Thomas Lumley
infjack.glm<-function(glm.obj,groups){
    umat<-estfun.glm(glm.obj)
    usum<-rowsum(umat,groups,reorder=FALSE)
    modelv<-summary(glm.obj)$cov.unscaled
    modelv%*%(t(usum)%*%usum)%*%modelv
}
jack.glm<-function(glm.obj,groups){
    umat<-jackvalues(glm.obj)
    usum<-rowsum(umat,groups,reorder=FALSE)
    t(usum)%*%usum*(nrow(umat)-1)/nrow(umat)
}
jackvalues<-function(glm.obj){
    db<-lm.influence(glm.obj)$coef
    t(t(db)-apply(db,2,mean))
}   
estfun.glm<-function(glm.obj){
    if (is.matrix(glm.obj$x)) 
        xmat<-glm.obj$x
    else {
        mf<-model.frame(glm.obj)
        xmat<-model.matrix(terms(glm.obj),mf)       
    }
    residuals(glm.obj,"working")*glm.obj$weights*xmat[,rownames(summary(glm.obj)$cov.unscaled)] # the column selection is added to fix a bug caused by weird formula
}


# goodness of fit
#  ngroups=10; main=""; add=FALSE; show.emp.risk=TRUE; lcol=1; ylim=NULL; weights=NULL
risk.cal=function(risk, binary.outcome, weights=NULL, ngroups=NULL, cuts=NULL, main="", add=FALSE, show.emp.risk=TRUE, lcol=2, ylim=NULL, scale=c("logit","risk")){
    
    scale=match.arg(scale)
    
    stopifnot(length(risk)==length(binary.outcome))
    if(!is.null(weights)) stopifnot(length(risk)==length(weights)) else weights=rep(1,length(risk))
    
    if(is.null(cuts)) {
        if(length(table(risk))<ngroups) ngroups=length(table(risk))
        risk.cat=as.numeric(Hmisc::cut2(risk,g=ngroups))
    } else {
        risk.cat=as.numeric(Hmisc::cut2(risk,cuts=cuts))
        print(table(risk.cat))
        ngroups=length(table(risk.cat))
    } 
    
    emp.risk=numeric(ngroups)
    for(i in 1:ngroups){
        emp.risk[i]=sum(binary.outcome[risk.cat==i])/sum(weights[risk.cat==i])
    }
    print(emp.risk)
    
    xx=cbind(risk, "emp.risk"=emp.risk[risk.cat])
    xx=xx[order(xx[,"risk"]),]
    
    if (is.null(ylim)) {
        if(show.emp.risk) all.=c(emp.risk,risk) else all.=risk
        if(scale=="logit") all.=setdiff(all.,0)
        ylim=range(all.,na.rm=T) 
    }
    ylab="risk"
    if (scale=="logit") {
        ylim=logit(ylim)
        xx=logit(xx)
        ylab="log odd"
    }

    plot(xx[,"risk"], type="l", xlab="", xaxt="n", main=main, col=1, ylim=ylim, ylab=ylab)
    if(show.emp.risk) lines(xx[,"emp.risk"], col=lcol, lty=2)
    
    print(table(xx[,"emp.risk"], useNA="ifany"))
        
    mylegend(col=1:lcol, lty=1:2, x=1, legend=c("model fit","empirical"))
    
}
