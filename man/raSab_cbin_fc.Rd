\name{raSab_cbin_fc}
\alias{raSab_cbin_fc}
\title{
Simulate a and Sab from full conditional distributions under the cbin likelihood
}
\description{
Simulate a and Sab from full conditional distributions under the cbin likelihood
}
\usage{
raSab_cbin_fc(Z, Y, a, b, Sab, odobs, odmax, SS = round(sqrt(nrow(Z))))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
a square matrix, the current value of Z 
}
  \item{Y}{
square matrix of ranked nomination data
}
  \item{a}{
current value of row effects
}
  \item{b}{
current value of column effects
}
  \item{Sab}{
current value of Cov(a,b)
}
  \item{odobs}{
observed outdegree
}
  \item{odmax}{
a scalar or vector giving the maximum number of nominations for each individual
}
  \item{SS}{
number of iterations
}
}
\value{
\item{Z}{new value of Z}
\item{Sab}{new value of Sab}
\item{a}{new value of a}
}
\author{
Peter Hoff
}

\examples{
## The function is currently defined as
function(Z,Y,a,b,Sab,odobs,odmax,SS=round(sqrt(nrow(Z))))
{
  E<-Z-a%*%t(rep(1,nrow(Z))) 
  MEL<-MEU<- -E
  MEL[!is.na(Y) & Y==0]<- -Inf
  MEU[!is.na(Y) & Y==1]<-Inf 
  MEL[is.na(Y)]<- -Inf ; MEU[is.na(Y)]<- Inf
  lba<-apply(MEL,1,max) 
  lba[is.na(lba)]<- -Inf
  uba<-apply(MEU,1,min) 
  uba[is.na(uba)]<- Inf

  uba[ odobs==odmax ]<- Inf  # to account for censoring

  for(ss in 1:SS)
  {
    ea<-b*Sab[1,2]/Sab[2,2]
    sa<-sqrt(Sab[1,1]-Sab[1,2]^2/Sab[2,2])
    a<-ea+sa*qnorm(runif(nrow(Z),pnorm((lba-ea)/sa),pnorm((uba-ea)/sa)))
    Sab<-solve(rwish(solve(diag(2)+crossprod(cbind(a,b))),3+nrow(Z)))
  }
  list(Z=E+a \%*\% t(rep(1,nrow(Z))),a=a,Sab=Sab)
}



}

