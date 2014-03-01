\name{amen-package}
\alias{amen-package}
\alias{amen}
\docType{package}
\title{
Additive and multiplicative effects modeling of networks and relational data
}
\description{
This package computes Bayesian estimates for additive and multiplicative effects (AME) models for networks and relational data. The basic model includes regression terms, the covariance structure of the social relations model (Warner, Kenny and Stoto (1979), Wong (1982)), and multiplicative factor effects (Hoff(2009)). Four different link functions accommodate different relational data structures, including binary/network data (bin), normal relational data (nrm), ordinal relational data (ord) and data from fixed-rank nomination schemes (frn).  Several of these link functions are discussed in Hoff, Fosdick, Volfovsky and Stovel (2013).Development of this software was supported in part by NICHD grant 1R01HD067509-01A1.
}
\details{
\tabular{ll}{
Package: \tab amen\cr
Type: \tab Package\cr
Version: \tab 0.999 \cr
Date: \tab 2014-3-1\cr
License: \tab GPL-3 \cr
}
}
\author{
Peter Hoff,
Bailey Fosdick, Alex Volfovsky, Kate Stovel

Maintainer: Peter Hoff <pdhoff@uw.edu>
}
\keyword{ package }
\examples{

data(YX_frn)
fit<-ame(YX_frn$Y,YX_frn$X,burn=5,nscan=5,odens=1,model="frn")

summary(fit)

plot(fit) 

}
