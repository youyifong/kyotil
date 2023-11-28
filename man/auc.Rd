\name{auc}
\alias{auc}
\alias{ratio}
\alias{trainauc}
\alias{coef.auc}
\alias{predict.auc}
\alias{print.auc}
\alias{summary.auc}
\alias{ratio.glm}
\alias{ratio.auc}
\alias{trainauc.glm}
\alias{trainauc.auc}

\title{AUC}
\description{
AUC methods.
}
\usage{
\method{coef}{auc}(object, \dots) 
\method{predict}{auc}(object, newdata, case.percentage = NULL, \dots) 
\method{print}{auc}(x, \dots) 
\method{summary}{auc}(object, \dots) 
\S3method{trainauc}{auc}(fit, training.data = NULL, \dots)
\S3method{ratio}{auc}(fit)

\S3method{trainauc}{glm}(fit, \dots)
\S3method{ratio}{glm}(fit)

}
\arguments{
  \item{fit}{an object that inherits from class 'auc' such as 'rauc' or 'sauc'}
  \item{object}{an object that inherits from class 'auc' such as 'rauc' or 'sauc'}
  \item{x}{an object that inherits from class 'auc' such as \code{rauc}, \code{sauc} or \code{sauc.dca}.}
  \item{newdata}{data at which to predict}
  \item{case.percentage}{used for class prediction, defaults to NULL}
  \item{training.data}{data frame used to compute auc based on a fit obtained by a call to  \code{rauc}, \code{sauc} or \code{sauc.dca}}
  \item{\dots}{arguments passed to or from methods}
}

%\value{}


\author{
  Youyi Fong \email{youyifong@gmail.com}\cr
  Krisztian Sebestyen \email{}\cr
}
