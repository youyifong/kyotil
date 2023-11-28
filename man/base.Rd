\name{base.functions}
\alias{cbinduneven}
\alias{binary}
\alias{binary2}
\alias{f2c}
\alias{ftoi}
\alias{keepWarnings}
\alias{meanmed}
\alias{methods4}
\alias{myaggregate}
\alias{myreshapewide}
\alias{mysapply}
\alias{myscale}
\alias{mysystem}
\alias{mytapply}
\alias{read.sv}
\alias{read.tsv}
\alias{table.prop}
\alias{unix}
\alias{table.cases}
\alias{table.cases.3}
\alias{multi.outer}
\alias{myreshapelong}
\alias{mycor}

\title{
Some Base Functions
}
\description{
cbinduneven binds together a list of matrixes/dataframes of different lengths, rows are matched by names
binary returns binary representation of an integer. 
binary2 returns binary representatin of an integer with leading 0, the length of string is n.
mysystem can call any exe file that is in the PATH
f2c convert temperature from f to c/



}
\usage{
cbinduneven(li)
binary(i)

multi.outer (f, ... ) 

myreshapelong(dat, cols.to.be.stacked, label.cols.to.be.stacked, new.col.name)

binary2(i, n)

f2c(f)

ftoi(f)

keepWarnings(expr)

meanmed(x, na.rm = FALSE)

methods4(classes, super = FALSE, ANY = FALSE)

myaggregate(x, by, FUN, new.col.name = "aggregate.value", ...)

myreshapewide(formula, dat, idvar, keep.extra.col=FALSE)

mysapply(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, ret.mat = TRUE)

myscale(x)

mysystem(cmd, ...)

mytapply(X, INDEX, FUN = NULL, ..., simplify = TRUE)

read.sv(file, header = TRUE, ...)

read.tsv(file, header = TRUE, sep = "\t", ...)

table.prop(x,y=NULL,digit=1,style=2,whole.table.add.to.1=FALSE,useNA="ifany",
    add.perc=FALSE, add.total.column = FALSE)

table.cases  (case,group,include.all=TRUE,desc="cases")
table.cases.3(case,group1,group2)

unix()

mycor (x, use = "everything", method = c("pearson", "kendall", "spearman"), 
    alternative = c("two.sided", "less", "greater"), exact = NULL, 
    conf.level = 0.95, continuity = FALSE, 
    digits.coef=2, digits.pval=3,
    ...) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{add.total.column}{
tbdi
}
  \item{use}{
tbdi
}
  \item{method}{
tbdi
}
  \item{alternative}{
tbdi
}
  \item{exact}{
tbdi
}
  \item{conf.level}{
tbdi
}
  \item{continuity}{
tbdi
}
  \item{digits.coef}{
tbdi
}
  \item{digits.pval}{
tbdi
}


  \item{cols.to.be.stacked}{
tbdi
}
  \item{label.cols.to.be.stacked}{
tbdi
}


  \item{li}{a list}  
  \item{i}{
tbdi
}

  \item{n}{
tbdn
}

  \item{f}{In multi.out, f is a function.}
  \item{case}{vector of 0/1}
  \item{group}{vector of multi-group indicators}
  \item{formula}{
a formula object.}


  \item{expr}{
tbdexpr
}

  \item{x}{
tbdx
}
  \item{na.rm}{
tbdna.rm
}

  \item{classes}{
tbdclasses
}
  \item{super}{
tbdsuper
}
  \item{ANY}{
tbdANY
}

  \item{desc}{
tbdby
}
  \item{by}{
tbdby
}
  \item{whole.table.add.to.1}{
Boolean
}
  \item{new.col.name}{
tbdnew.col.name
}
  \item{\dots}{
tbd\dots
}

  \item{dat}{
tbddat
}
  \item{idvar}{
tbdidvar
}

  \item{X}{
tbdX
}
  \item{simplify}{
tbdsimplify
}
  \item{USE.NAMES}{
tbdUSE.NAMES
}
  \item{ret.mat}{
tbdret.mat
}


  \item{cmd}{
tbdcmd
}

  \item{INDEX}{
tbdINDEX
}

  \item{file}{
tbdfile
}
  \item{header}{
tbdheader
}

  \item{sep}{
tbdsep
}

  \item{y}{
tbdy
}
  \item{digit}{
tbddigit
}
  \item{style}{
tbdstyle
}

  \item{FUN}{
tbdFUN
}
  \item{keep.extra.col}{
tbdFUN
}
  \item{useNA}{
tbdFUN
}
  \item{add.perc}{
tbdFUN
}
  \item{include.all}{
tbdFUN
}
  \item{group1}{
tbdFUN
}
  \item{group2}{
tbdFUN
}

}
\examples{

binary(5) ###  101
binary2(5, 4)

a=data.frame("x"=1:2)
b=data.frame("y"=3:5);#rownames(b)[3]=""
cbinduneven(list(a,b))

\dontrun{
# the formula in myreshapewide can only have one variable in the right hand side
    myreshapewide(fi~week, dat, c("ptid","stim"))

    myreshapelong(dat.201.neut, cols.to.be.stacked=c("MN.3","SF162","SVA.MLV"), 
        label.cols.to.be.stacked="antigen", new.col.name="y")

    myaggregate(subset(dat.poc, select=c(HIV, trt)), list(dat.poc$f), function(x) 
      with(x, c(fisher.test(HIV, trt)$estimate, fisher.test(HIV, trt)$p.value)))


}

}
