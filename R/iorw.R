# @ Youyi Fong, Fred Hutchinson Cancer Research Center
# Modified from code written by Ben Cowling and Wey Wen Lim (March 30, 2018) downloaded from https://datadryad.org/stash/dataset/doi:10.5061/dryad.cv37539

iorw=function(formula.effect, formula.mediators, data, family=NULL, nboot=10000, numCores=1, save.steps=FALSE, verbose=FALSE) {
    
    # only keep complete.case
    nomissing=complete.cases(model.frame(formula.effect, data, na.action=na.pass)) & complete.cases(model.frame(formula.mediators, data, na.action=na.pass))
    data=data[nomissing,]
    n <- nrow(data)
    
    if (is.null(family)) {
        if (length(terms(formula.effect)[[2]]==3)) {
            # Surv(,)
            coxreg=TRUE
        } else stop("Please specify a family")
    } else coxreg=FALSE
    
    total.regression=function(dat) {
        if (coxreg) {
            # bootstrap datasets sometime will results in problems that only shows up as warnings
            res=keepWarnings(survival::coxph(formula.effect, data=dat))
            out=res$value
            attr(out,"exception.thrown") = length(res$warning)>0
            out
        } else {
            out=glm(formula.effect, data=dat, family=family)     
            attr(out,"exception.thrown") = FALSE
            out
        }        
    }
    
    # make formula.mediators
    tmp=strsplit(as.character(formula.effect)[3], "\\+")[[1]]
    trt=tmp[1] # first covariate on the right hand side is assumed to be treatment
    formula.mediators=update(formula.mediators, as.formula(paste0(trt, "~.+", concatList(tmp[-1],sep="+"))))
    
    if (verbose) {
        print(formula.effect)
        print(formula.mediators)
    }
    
    do.est=function(dat, save=FALSE) {
        # Step 1: regression to estimate total effect 
        total.model <- total.regression(dat)
        if(attr(total.model, "exception.thrown")) return (c(total=NA, ve=NA, direct=NA, indirect=NA, prop=NA))
        
        # Step 2: regress treatment ~ mediators + covariates, to generate weights for direct effect estimation
        # weight is inversely associated with odds of being in trt=1
        temp.lreg <- glm(formula.mediators, data=dat, family="binomial") 
        iorw.weight=NULL # this is to trick R CMD check so that we don't get: Error in eval(extras, data, env) : object 'iorw.weight' not found
        dat$iorw.weight <- (1/exp(predict(temp.lreg, data=dat, type="link")))^dat[[trim(trt)]] 
#        dat$iorw.weight = dat$iorw.weight^0.3
        
                
        # Step 3: weighted regression to estimate direct effect
        if (coxreg) {
            direct.model <- survival::coxph(formula.effect, data=dat, weights=iorw.weight) 
        } else {
            direct.model <- suppressWarnings(glm(formula.effect, data=dat, family=family, weights=iorw.weight) )
            # warnings: In eval(family$initialize) : non-integer #successes in a binomial glm!
        }
        
        # the hazard ratio for treatment without weights
        total= unname(exp(coef(total.model)))[1+!coxreg] #cox reg no intercept, so it is the first
        # Vaccine efficacy = 1-total effect
        ve=1-total
        # the hazard ratio for treatment with weights
        direct=unname(exp(coef(direct.model)))[1+!coxreg]
        # indirect effect is the ratio of total effect and direct effect
        indirect=total/direct
        # proportion of vaccine effect explained by marker
        prop=log(indirect) / log(total)
        
        out=c(total=total, ve=ve, direct=direct, indirect=indirect, prop=prop)     
        if(save) {
            attr(out,"step1fit")=total.model
            attr(out,"step2fit")=temp.lreg
            attr(out,"step3fit")=direct.model
            attr(out,"step2w")=dat$iorw.weight
        }        
        out                   
    }
    
    out=list()    
    out$formula.effect=formula.effect
    out$formula.mediators=formula.mediators
    tmp=do.est(data, save.steps)
    out$step1fit=attr(tmp,"step1fit")
    out$step2fit=attr(tmp,"step2fit")
    out$step3fit=attr(tmp,"step3fit")
    out$step2weights=attr(tmp,"step2w")
    attributes(tmp)=attributes(tmp)[1]
    out$est=tmp
    
    #### bootstrap ####
    if (nboot>0) {
        # save rng state before set.seed in order to restore before exiting this function
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      
        boot.out=parallel::mclapply(1:nboot, mc.cores = numCores, FUN=function(seed) {    
            set.seed(seed)
            new.data <- data[sample.int(n, replace=TRUE),] 
            do.est(new.data)
        })
        boot.out=do.call(cbind, boot.out) 
        # restore rng state 
        assign(".Random.seed", save.seed, .GlobalEnv)     
    
        out$boot.perc=apply(boot.out, 1, function(x) quantile(x, c(0.025, 0.975), na.rm=TRUE)) 
        attr(out$boot.perc,"nboot") = nboot
        attr(out$boot.perc,"na.cnt") = sum(is.na(boot.out[1,]))
    }
    
    class(out)=c("iorw", class(out))
    
    out
    
}

print.iorw=function(x, ...) {
    for (a in names(x)) {
        if (a=="step2weights") {
            cat("\n#### ", a, ": \n", sep="")
            print(summary(x[[a]]))
        } else {
            cat("\n#### ", a, ": \n", sep="")
            print(x[[a]])
        }
    }
}
