\name{probit_start}
\alias{probit_start}
\title{
Starting values for probit SRM
}
\description{
Starting values for probit social relations model
}
\usage{
probit_start(Y, X)
}
\arguments{
  \item{Y}{
n x n binary relational matrix 
}
  \item{X}{
n x n x p covariate array
}
}
\value{
\item{beta}{probit estimate of beta}
\item{a}{probit estimate of row effects}
\item{b}{probit estimate of column effects}
\item{rho}{probit estimate of dyadic correlation}
}
\author{
Peter Hoff
}

\examples{
## The function is currently defined as
probit_start <-
function(Y,X)
{
  fit<-glm(c(1*(Y>0))~ -1 + apply(X,3,c),family=binomial(link=probit))
  E<-matrix(0,nrow(Y),ncol(Y)) ; E[!is.na(Y)]<-fit$res

  a<-apply(E,1,mean) ; b<-apply(E,2,mean)
  E<-E - outer(a,b,"+")
  rho<-cor(cbind(E[upper.tri(E)],t(E)[upper.tri(E)]))[1,2]
  beta<-fit$coef
  list(beta=beta,a=a,b=b,rho=rho)
}
}
