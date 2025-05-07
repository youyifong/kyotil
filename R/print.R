mykable <- function(x, ...) {
  # If input is a matrix, convert to data frame first
  if (is.matrix(x)) {
    rown <- rownames(x)
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    rownames(x) <- rown
  }
  
  # Save row names
  rown <- rownames(x)
  
  # Replace NA and NaN with empty strings
  x_clean <- data.frame(lapply(x, function(col) {
    ifelse(is.na(col), "", col)
  }), stringsAsFactors = FALSE)
  
  # Restore row names
  rownames(x_clean) <- rown
  
  # Pass to kable
  kable(x_clean, ...)
}


mytex=function(dat=NULL, file.name="temp", 
    digits=NULL, display=NULL, align="r", 
    include.rownames=TRUE, include.colnames=TRUE,
    col.headers=NULL,
    comment=FALSE, floating=FALSE, 
    lines=TRUE, hline.after=NULL, 
    add.to.row=NULL, 
    sanitize.text.function = NULL, #function(x) x,
    append=FALSE, preamble="", input.foldername=NULL, save2input.only=NULL,
    caption=NULL, label=paste("tab",mylast(strsplit(file.name, "/")[[1]]),sep=" "), table.placement="h!",
    add.clear.page.between.tables=FALSE,
    longtable=FALSE,
    verbose=FALSE,
	silent=TRUE,
...) {
# file.name="temp"; digits=NULL; display=NULL; align="r"; include.rownames=TRUE; include.colnames=TRUE;    col.headers=NULL;    comment=FALSE; floating=FALSE;     lines=TRUE; hline.after=NULL;add.to.row=NULL; sanitize.text.function = NULL;append=FALSE; preamble="";  save2input.only=NULL;    caption=NULL; table.placement="h!";    add.clear.page.between.tables=FALSE;    longtable=FALSE;    verbose=FALSE;	silent=TRUE; save2input.only=T; label=""
  # input.foldername=NULL;

#    if(exists("tablePath") && file.exists(tablePath)) {
#        file.name=tablePath%.%"/"%.%file.name
#    } else {
#        file.name=file.name
#    }    
        
    if(is.null(save2input.only)) save2input.only = !is.null(input.foldername)
    if (endsWith(file.name,".tex")) file.name=substr(file.name, 1, nchar(file.name)-4)
    tmp=strsplit(file.name, split="/")[[1]] 
    if(is.null(input.foldername)) {
      if (length(tmp)>1) input.foldername=concatList(tmp[-length(tmp)], "/")%.%"/input" else input.foldername="input"
    }
    dir.exists(input.foldername) || dir.create(input.foldername, recursive = TRUE)
    file.name.2=input.foldername%.%"/"%.%tmp[length(tmp)]
    
    if(is.data.frame(dat)) dat=list(dat)
    if (!is.list(dat)) dat=list(dat)
    
    if (!append) { #start a new file
        #document tag, preamble etc
        if(!save2input.only) mytex.begin(file.name%.%".tex", preamble)
        #empty file
        cat ("", file=file.name.2%.%".tex", append=FALSE)
    } 
    
    add.to.row.0=add.to.row
    include.colnames.0=include.colnames
    include.rownames.0=include.rownames
    align.0=align

    if (length(dat)==0) stop("length of dat is 0")
    names(dat)=gsub("_"," ",names(dat))
    for (i in 1:length(dat)) {
    
        include.rownames=include.rownames.0
        align=align.0
        dat1 = dat[[i]]   
        if (is.null(dat1)) warning("some element of dat list is null")   
         
        # character only 
        if (!is.matrix(dat1) & is.character(dat1)) {
            if(!save2input.only) cat (dat1%.%"\n\n\n", file=file.name%.%".tex", append=TRUE)
            cat (dat1%.%"\n\n\n", file=file.name.2%.%".tex", append=TRUE)
            next
        }
            
        # convert vector to matrix
        if (length(dim(dat1))==1) dat1=matrix(c(dat1),nrow=1, dimnames=list(NULL,names(dat1)))
        #if (is.vector(dat1)) dat1=as.matrix(dat1)
    
        dimnam=names(dimnames(dat1))
        
        .ncol=ncol(dat1)
        if(verbose) myprint(align)
        if (length(align)==1) {
            align=rep(align,.ncol+1); align[1]="l"
        } else if (length(align)==.ncol) {
            align=c("l",align) #align may not include alignment for the rownames, just pad it
        } else if (length(align)!=.ncol+1) {
            str(align); str(dat1); myprint(.ncol); stop("length of align incorrect")
        }
        if(verbose) myprint(align)
        
        # add rownames as the first column if necessary
        if(include.rownames & anyDuplicated(rownames(dat1))) {
            tmp=suppressWarnings(data.frame(rownames(dat1),data.frame(dat1))) # may generate warnings about duplicate row names
            names(tmp)[1]=""
            if (!is.null(colnames(dat1))) colnames(tmp)[2:ncol(tmp)]=colnames(dat1)
            if (!is.null(dimnam)) if (is.na(dimnam[2])) colnames(tmp)[1]=dimnam[1] else if (dimnam[2]=="") colnames(tmp)[1]=dimnam[1] 
            dat1=tmp
            include.rownames=FALSE
            if (length(align)==ncol(dat1)) align=c("l",align) # need to extend align on the left
            .ncol=.ncol+1
            #str(align)
        } 

        if (is.null(digits)) if (is.integer(dat1)) digits=0
        
        if(!is.null(dimnam) & is.null(col.headers)) {
            if(!is.na(dimnam[2]))
              if(trim(dimnam[2])!="")
                col.headers="\\hline\n  "%.%ifelse(include.rownames,dimnam[1]%.%"&","")%.%"  \\multicolumn{"%.% ncol(dat1) %.%"}{c}{"%.%dimnam[2]%.%"}   \\\\  \n"
        }
        if (!is.null(col.headers)) top=col.headers else top="\\hline  "
        
        if(include.colnames & !is.null(colnames(dat1))) {
            # to make border in the column names, but centrally aligned
            coln=if(!is.null(sanitize.text.function)) sanitize.text.function(colnames(dat1)) else mysanitize.text(mysanitize.numbers(colnames(dat1)))
            align.1=align[-1]
            align.1=gsub("[lr]","c",align.1)
            # multicolumn env makes sanitize.text result not compilable
            top.1=concatList(coln, sep="& ") %.% "\\\\ \n"%.% # center aligned column titles
#            top.1=concatList(" \\multicolumn{1}{"%.%align.1%.%"}{"%.%coln%.%"} ", sep="&") %.% "\\\\ \n"%.% # center aligned column titles
                "\\hline\n" # insert at the beginning of table, "\n" is added so that there is no need to keep it in col.title
            #print(coln);print(top.1);print(align.1);print(align)
            
            # add a column for rownames, which may include names of rownames
            if(include.rownames) {
                top.1="&" %.% top.1
                if (!is.null(dimnam)) if (is.na(dimnam[2])) top.1=dimnam[1] %.% top.1 else if (trim(dimnam[2])=="") top.1=dimnam[1] %.% top.1
            }
            top=top%.%top.1
            #print(coln);print(top.1);print(align.1)
        }
            
        if (include.colnames & (is.null(hline.after) | length(dat)>1) ) hline.after=c(nrow(dat1)) # cannot use default due to add.to.row    
        include.colnames=FALSE        
        
        if (is.null(add.to.row)) {
            add.to.row=list(list(0), top)
        } else {
            add.to.row=list(c(list(0), add.to.row[[1]]), c(top, add.to.row[[2]]))
        }
        #print(add.to.row)
 
        if (length(dat)>1) {
            if(!save2input.only) cat (ifelse(add.clear.page.between.tables, names(dat)[i]%.%"\n\n", "\\vspace{20pt}"%.%names(dat)[i]%.%"\n\n"), file=file.name%.%".tex", append=TRUE)
            cat (ifelse(add.clear.page.between.tables, names(dat)[i]%.%"\n\n", "\\vspace{20pt}"%.%names(dat)[i]%.%"\n\n"), file=file.name.2%.%".tex", append=TRUE)
        }
        #if (!is.null(attr(dat1,"caption"))) caption=attr(dat1,"caption") else caption=NULL
        
        if (is.null(hline.after)) {
            if (lines) hline.after=c(-1,0,nrow(dat1)) else hline.after=c(nrow(dat1))
            if (!include.colnames) hline.after=hline.after[-(1:2)]
        }
        #print(hline.after)
        if (verbose) print(file.name%.%".tex")
        
        if(!include.rownames) rownames(dat1)=1:nrow(dat1)# otherwise there will be a warning from xtable
        if(!save2input.only) print(..., xtable::xtable(dat1, 
                digits=(if(is.null(digits)) rep(3, .ncol+1) else digits), # cannot use ifelse here!!!
                display=(if(is.null(display)) rep("f", .ncol+1) else display), # or here
                align=align, caption=caption, label=label, ...), # for caption to work, floating needs to be T
            hline.after=hline.after, type = "latex", file = file.name%.%".tex", append = TRUE, floating = floating, table.placement=table.placement, 
            include.rownames=include.rownames, include.colnames=include.colnames, comment=comment, 
            add.to.row=add.to.row, sanitize.text.function =sanitize.text.function )

        print(..., xtable::xtable(dat1, 
                digits=(if(is.null(digits)) rep(3, .ncol+1) else digits), # cannot use ifelse here!!!
                display=(if(is.null(display)) rep("f", .ncol+1) else display), # or here
                align=align, caption=caption, label=label, ...), 
            hline.after=hline.after, type = "latex", file = file.name.2%.%".tex", append = TRUE, floating = floating, table.placement=table.placement, 
            include.rownames=include.rownames, include.colnames=include.colnames, comment=comment, 
            add.to.row=add.to.row, sanitize.text.function =sanitize.text.function, tabular.environment=ifelse(longtable, "longtable","tabular"))
        
        if(i!=length(dat) & add.clear.page.between.tables) {
          cat ("\\clearpage\n", file=file.name.2%.%".tex", append=TRUE)
        }

        if(!save2input.only) cat ("\n", file=file.name%.%".tex", append=TRUE)
        #cat ("\n", file=file.name.2%.%".tex", append=TRUE) # don't add this line since extra lines at the end will prevent two tabular from being put on the same line
        # restore some variables that have changed in this function
        add.to.row=add.to.row.0
        include.colnames=include.colnames.0
    }
    
    if(!append) {
        if(!save2input.only) mytex.end(file.name%.%".tex")
    }
    if(!save2input.only & !silent) cat ("Writing table to "%.%getwd()%.%"/"%.%file.name%.%"\n")
}
#x=matrix(0,2,2,dimnames=list(a=1:2, b=1:2));  mytex(x)
#x=matrix(0,2,2,dimnames=list(a=1:2, 1:2));  mytex(x)

# print a matrix/table or a list of them to a latex file as xtable
# note file.name can not have space in it
# e.g. mytex(matrix(0,2,2));
# e.g. mytex(matrix(0,2,2), digits=4);
# e.g. mytex(list(matrix(0,2,2), c(1,1))); 
# default arguments: file.name="temp"; digits=NULL; display=NULL; align="r"; append=FALSE; preamble=""; keep.row.names=TRUE
mytex.begin=function(file.name,preamble=""){
#    if(exists("tablePath") && file.exists(tablePath)) {
#        file.name=tablePath%.%"/"%.%file.name
#    } else {
#        file.name=file.name
#    }    
    if (!endsWith(file.name,".tex")) file.name=file.name%.%".tex"
    cat ("\\documentclass{article}\n", file=file.name, append=FALSE)
    cat (preamble, file=file.name, append=TRUE)
    cat("\n\\usepackage{geometry}\n", file=file.name, append=TRUE)    
    cat("\n\\begin{document}\n", file=file.name, append=TRUE)    
}
mytex.end=function(file.name){
#    if(exists("tablePath") && file.exists(tablePath)) {
#        file.name=tablePath%.%"/"%.%file.name
#    } else {
#        file.name=file.name
#    }    
    if (!endsWith(file.name,".tex")) file.name=file.name%.%".tex"
    cat ("\n\\end{document}", file=file.name, append=TRUE)
}

# adapted from print.xtable.R
mysanitize.text <- function(str) {
    result <- str
    result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
    result <- gsub("$", "\\$", result, fixed = TRUE)
    result <- gsub(">", "$>$", result, fixed = TRUE)
    result <- gsub("<", "$<$", result, fixed = TRUE)
    result <- gsub("|", "$|$", result, fixed = TRUE)
    result <- gsub("{", "\\{", result, fixed = TRUE)
    result <- gsub("}", "\\}", result, fixed = TRUE)
    result <- gsub("%", "\\%", result, fixed = TRUE)
    result <- gsub("&", "\\&", result, fixed = TRUE)
    result <- gsub("_", "\\_", result, fixed = TRUE)
    result <- gsub("#", "\\#", result, fixed = TRUE)
    result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
    result <- gsub("~", "$\\sim$", result, fixed = TRUE) # this is changed by Y.F. 
    result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$",
                   result, fixed = TRUE)
    return(result)
}
mysanitize.numbers <- function(x) {
    result <- x
#    if ( math.style.negative ) {
        ## Jake Bowers <jwbowers@illinois.edu> in e-mail
        ## from 2008-08-20 suggested disabling this feature to avoid
        ## problems with LaTeX's dcolumn package.
        ## by Florian Wickelmaier <florian.wickelmaier@uni-tuebingen.de>
        ## in e-mail from 2008-10-03 requested the ability to use the
        ## old behavior.
        for(i in 1:length(x)) {
            result[i] <- gsub("-", "$-$", result[i], fixed = TRUE)
#        }
    }
    return(result)
}

        
# write a table that contains mean and sd to temp.tex in the current working directory, getwd()
# models can be a list of models, or a single model
make.latex.coef.table = function (models, model.names=NULL, row.major=FALSE, round.digits=NULL) {
# e.g.: models=list(gam1, gam2); round.digits= c(3,3,3,3,3); model.names=c("gam1", "gam2");  row.major=TRUE   
    if (! ("list" %in% class (models) ) ) {models=list(models)}
    
    numParams = nrow (getFixedEf(models[[1]]))
    numModels = length (models)
    
    if (is.null (model.names)) {model.names=rep("",numModels)}
    if (is.null(round.digits)) round.digits=rep(3,numParams)    
    
    coef.table = mysapply (1:numModels, function (i.model) {
        temp = getFixedEf(models[[i.model]]) [,1:2,drop=FALSE]
        for (i.param in 1:numParams) {
            temp[i.param,] = round (temp[i.param,], round.digits[i.param])
        }
        temp2 = paste (format(temp[,1]), "(", format(temp[,2]), ")")
        names (temp2) = dimnames(temp)[[1]]
        temp2
    })
    dimnames (coef.table)[[1]] = model.names
    
    if (row.major) mytex ( coef.table, align="r" ) 
    else mytex (t(coef.table), align="r") 
}


roundup=function (value, digits, na.to.empty=TRUE, remove.leading0=FALSE) {
    if (length(digits)==1) {
        out=format(round(value, digits), nsmall=digits, scientific=FALSE) 
    } else {
        if (length(digits)!=length(value)) stop("length of value and length of values different")
        out = sapply(1:length(digits), function (i) roundup (value[i], digits[i], na.to.empty))
    }
    if(remove.leading0) out=sub("^0\\.","\\.",out)
    if(na.to.empty) sub("NA|NaN","",out) else out
}
formatInt=function (x, digits, fill="0", ...) {
    formatC(x, format="d", flag=fill, width=digits) 
}
formatDouble=roundup
# test
#formatDouble(c(1,2,3.12344), 1:3)


prettyprint=function (value, digit=2) {
  if (value>=1e4) {
    gsub("e\\+0*", "x10^", format(value, digit=digit, scientific=T) )
  } else if (value>=100) {
      round(value)
  } else if (value>=1) {
      round(value,1)
  } else {
      signif(value, digit)
  }
}


# don't have to transpose x
mywrite=function(x, ...){
    if (is.list(x)) x=fill.jagged.array(x)
    if (is.null(ncol(x))) i=length(x)
    else i=ncol(x)
    write (t(x), ncolumns=i, ...)
}

# default row.names to FALSE
# file name needs no file extension
mywrite.csv = function(x, file="tmp", row.names=FALSE, digits=NULL, silent=TRUE, ...) {  
    if (!is.null(digits)) {
        if(length(digits)==1) {
            x=round(x,digits)
        } else {
            for (i in 1:ncol(x)) {
                x[,i]=round(x[,i], digits[i])
            }                
        }
    }
    x[is.na(x)]=""
    write.csv(x, file=file%.%".csv", row.names=row.names, ...)
    if(!silent) cat("Writing table to "%.%getwd()%.%"/"%.%file%.%".csv\n")
}


myprint <- function(object, ...) UseMethod("myprint") 

# this function is placed at the bottom of the file because it contains "\""), which makes all the following line be miss-interpreted as being in quotes
myprint.default = function (..., newline=TRUE, digits=3, print.name=TRUE) {   
    digits.save=getOption("digits")
    options(digits=digits)
    object <- as.list(substitute(list(...)))[-1]
    x=list(...)
    for (i in 1:length(x)) {
        if (inherits(x[[i]],"formula")) {cat(as.character(x[[i]]), "; "); next}
        tmpname <- deparse(object[[i]])[1]
        #str(tmpname)
        #str(gsub("\\\\","\\",gsub("\"", "", tmpname)))
        #str(x[[i]])
        #if (gsub("\\\\","\\",gsub("\"", "", tmpname))!=x[[i]]) {
        
#        # the following line fails when I redefined contain using grepl
#        if (contain(tmpname, "\"") | contain(tmpname, "\\")) {
#            for (a in x[[i]]) cat(a)
#        } else {
            if(print.name) cat (tmpname %.% " = ")
            for (a in x[[i]]) cat(a,"") # by putting "" there, a space is introduced b/c cat prints a sep
            if (i!=length(x)) cat ("; ")
#        }
    }
    if (newline)  cat("\n")
    options(digits=digits.save)
}
#a="hello"; b="you"; myprint (a, b); myprint ("test"); myprint.default ("\t")

myprint.matrix=function(object, ...){
    tmpname <- deparse(substitute(object))
    cat(tmpname, "\n")
    for (i in 1:nrow(object)) myprint(object[i,], print.name=FALSE, ...)
}
