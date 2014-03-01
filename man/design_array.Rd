\name{design_array}
\alias{design_array}
\title{
Computes the design socioarray of covariate values
}
\description{
Computes the design socioarray of covariate values for an AME fit
}
\usage{
design_array(Xrow=NULL,Xcol=NULL,Xdyad=NULL,intercept=TRUE,n)
}
\arguments{
\item{Xrow}{ an n x pr matrix of row covariates }
\item{Xcol}{ an n x pc matrix of column covariates  }
\item{Xdyad}{ an n x n x pd array of dyadic covariates }
\item{intercept}{ logical }
\item{n}{ number of rows/columns }
}
\value{
an n x n x (pr+pc+pd+intercept) 3-way array 
}
\author{
Peter Hoff
}
\examples{
## The function is currently defined as 
design_array(intercept=TRUE,n=10) 
}
