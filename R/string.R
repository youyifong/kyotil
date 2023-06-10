# paste two strings together
# e.g. "a" %.% "b"

#"%+%" <- function (a, b) {
#    .Deprecated("%.%")
#    out=paste(a,b,sep="")
#    out
#}
"%.%" <- function (a, b) out=paste(a,b,sep="")




#### file name manipulation

getFileStem=function(file.name){
    substr(file.name, 1, lastIndex(file.name, ".")-1 )
}
getStem=getFileStem
fileStem=getFileStem

getExt = function(file.name){
    substr(file.name, lastIndex(file.name, ".")+1, nchar(file.name) )
}

# s="prefix_name" is changed to "name"
remove.prefix=function(s,sep="_"){
    tmp=strsplit(s,sep)    
    sapply(tmp, function (x) {
        if (length(x)==1) return (x)
        concatList(x[-1],sep)
    })
}


remove.postfix=function(s,sep="_"){
    tmp=strsplit(s,sep)    
    sapply(tmp, function (x) {
        if (length(x)==1) return (x)
        concatList(x[-length(x)],sep)
    })
}



#### misc

escapeUnderline=function (name) {
    gsub("_", "\\\\_", name)
}



#### string search

firstIndex =function (s1, s2) {
    k=nchar (s2)
    ret=-1;
    for (i in 1:(nchar(s1)-k+1) ) {
        if (substr(s1, i, i+k-1)==s2) {
            ret=i;
            break;
        }
    }
    ret
}

lastIndex =function (s1, s2) {
    k=nchar (s2)
    ret=-1;
    for (i in 1:(nchar(s1)-k+1) ) {
        if (substr(s1, i, i+k-1)==s2) {
            ret=i;
        }
    }
    ret
}

# return TRUE if s1 contains s2
contain =function (s1, s2) {
    grepl(s2, s1)
#    sapply (s1, function (s) {
#        if (is.na(s)) return (NA)
#        k=nchar (s2)
#        matched=0
#        for (i in 1:(nchar(s)-k+1) ) {
#            if (substr(s, i, i+k-1)==s2) {
#                matched=i
#                break
#            }
#        }
#        matched!=0
#    })
}
## test
#contain("abc","bc")
#grepl("bc",c("abc","abc"))

# copied from http://stackoverflow.com/questions/2261079/how-to-trim-leading-and-trailing-whitespace-in-r
trim.leading <- function (x)  sub("^\\s+", "", x)
trim.trailing <- function (x) sub("\\s+$", "", x)
trim <- function (x, trim.trailing=TRUE, trim.leading=TRUE)  {
    if(trim.trailing & trim.trailing) gsub("^\\s+|\\s+$", "", x) else if (trim.trailing) trim.trailing(x) else if (trim.trailing) trim.leading(x) else x
}
