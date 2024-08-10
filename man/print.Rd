\name{print.functions}
\alias{myprint}
\alias{myprint.default}
\alias{myprint.matrix}
\alias{formatInt}
\alias{prettyprint}
\alias{formatDouble}
\alias{make.latex.coef.table}
\alias{mytex}
\alias{mytex.begin}
\alias{mytex.end}
\alias{mywrite}
\alias{mywrite.csv}
\alias{roundup}
\alias{mysanitize.numbers}
\alias{mysanitize.text}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print Functions
}
\description{
roundup prints a specified number of digits after decimal point even if 0s are needed at the end. 
formatInt prints a specified number of digits before decimal point even if 0s are needed at the beginning.
}
\usage{

myprint(object, ...)

\method{myprint}{default}(..., newline = TRUE, digits = 3, print.name=TRUE)

\method{myprint}{matrix}(object, ...)

formatInt(x, digits, fill = "0", ...)

prettyprint (value, digit=2) 

make.latex.coef.table(models, model.names = NULL, row.major = FALSE, round.digits = NULL)

mysanitize.text(str) 
mysanitize.numbers(x) 

mytex(dat = NULL, file.name = "temp", digits = NULL, display
                 = NULL, align = "r", include.rownames = TRUE,
                 include.colnames = TRUE, col.headers = NULL, comment =
                 FALSE, floating = FALSE, lines = TRUE, hline.after =
                 NULL, add.to.row = NULL, sanitize.text.function =
                 NULL, append = FALSE, preamble = "", input.foldername
                 = NULL, save2input.only = NULL, caption = NULL, label
                 = paste("tab", last(strsplit(file.name, "/")[[1]]),
                 sep = " "), table.placement = "h!",
                 add.clear.page.between.tables = FALSE, longtable =
                 FALSE, verbose = FALSE, silent = TRUE, ...)

mytex.begin(file.name, preamble = "")

mytex.end(file.name)

mywrite(x, ...)

mywrite.csv(x, file = "tmp", row.names = FALSE, digits = NULL,
                 silent = TRUE, ...)

roundup (value, digits, na.to.empty=TRUE, remove.leading0=TRUE) 

formatDouble(value, digits, na.to.empty=TRUE, remove.leading0=TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{digit}{
tbddigit
}
  \item{silent}{
tbdnewline
}
  \item{input.foldername}{
tbdnewline
}
  \item{object}{
tbdnewline
}
  \item{newline}{
tbdnewline
}
  \item{print.name}{
tbddigits
}


  \item{save2input.only}{Boolean}
  \item{include.colnames}{Boolean}
  \item{col.headers}{string. Column headers}
  \item{comment}{Boolean, whether to include the version and timestamp comment}
  \item{hline.after}{vector}
  \item{add.to.row}{a list}
  \item{sanitize.text.function}{a function}

  \item{str}{
tbdvalue
}
  \item{remove.leading0}{
tbdvalue
}
  \item{caption}{
tbdvalue
}
  \item{longtable}{
tbdvalue
}
  \item{label}{ default to be the same as file.name stem
}
  \item{table.placement}{
tbdvalue
}
  \item{na.to.empty}{
tbdvalue
}
  \item{value}{
tbdvalue
}
  \item{digits}{
tbddigits
}
  \item{fill}{
tbdfill
}

  \item{models}{
tbdmodels
}
  \item{model.names}{
tbdmodel.names
}
  \item{row.major}{
tbdrow.major
}
  \item{round.digits}{
tbdround.digits
}

  \item{dat}{
tbddat
}
  \item{file.name}{
tbdfile.name
}
  \item{display}{
tbddisplay
}
  \item{align}{
tbdalign
}
  \item{append}{
tbdappend
}
  \item{preamble}{
tbdpreamble
}
  \item{include.rownames}{
tbdinclude.rownames
}
  \item{floating}{
tbdfloating
}
  \item{lines}{
tbdlines
}
  \item{\dots}{
tbd\dots
}
  \item{verbose}{
tbd\dots
}



  \item{x}{
tbdx
}

  \item{file}{
tbdfile
}
  \item{row.names}{
tbdrow.names
}

  \item{add.clear.page.between.tables}{
tbdrow.names
}
}
\examples{

roundup (3.1, 2) # 3.10

formatInt(3, 2) # 03


\dontrun{ 

# demo of dimnames
tab=diag(1:4); rownames(tab)<-colnames(tab)<-1:4; names(dimnames(tab))=c("age","height")
# for greek letter in the labels, we need sanitize.text.function=identity
rownames(tab)[1]="$\\alpha$"
# note that to use caption, floating needs to be TRUE
mytex (tab, file="tmp1", sanitize.text.function=identity, 
    caption="This is a caption .........................", caption.placement="top",
    floating=TRUE)

# col.headers has to have the RIGHT number of columns 
# but align is more flexible, may not need to include the rownames col
tab=diag(1:4); rownames(tab)<-colnames(tab)<-1:4
mytex (tab, file="tmp", include.rownames = TRUE, 
    align=c("c","c","c|","c","c"), col.headers=
    "\\hline\n & \\multicolumn{2}{c|}{Vaccine} & \\multicolumn{2}{c}{Control} \\\\ \n")
# not include rownames
mytex (tab, file="tmp", include.rownames = FALSE, 
    align=c("c","c","c|","c","c"), col.headers=
    "\\hline\n     \\multicolumn{2}{c|}{Vaccine} & \\multicolumn{2}{c}{Control} \\\\ \n") 
# It should work even if some rownames are duplicated
tab=diag(1:4); rownames(tab)=rep(1,4); colnames(tab)<-1:4
mytex (tab, file="tmp", include.rownames = TRUE, 
    align=c("c","c|","c","c"), col.headers=
    "\\hline\n & \\multicolumn{2}{c|}{Vaccine} & \\multicolumn{2}{c}{Control} \\\\ \n") 


# add.to.rows
tab=diag(1:4); rownames(tab)<-1:4; colnames(tab)<-c("a","b","c","d")
mytex (tab, file="tmp",
    add.to.row=list( list(0,2),
        c("          \\multicolumn{5}{l}{Heading 1} \\\\ \n",
          "\\hline\n \\multicolumn{5}{l}{Heading 2}\\\\ \n"
    ))
)


}



}
