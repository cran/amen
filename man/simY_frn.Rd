\name{simY_frn}
\alias{simY_frn}
\title{
Simulate an relational matrix based on a fixed rank nomination scheme
}
\description{
Simulate an relational matrix based on a fixed rank nomination scheme
}
\usage{
simY_frn(EZ, rho, odmax, YO)
}
\arguments{
  \item{EZ}{
a square matrix giving the expected value of the latent Z matrix
}
  \item{rho}{
dyadic correlation
}
  \item{odmax}{
a scalar or vector giving the maximum number of nominations for each node
}
\item{YO}{
a square matrix identifying where missing values should be maintained
}
}
\value{
a square matrix, where higher values represent stronger relationshis
}
\author{
Peter Hoff
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(EZ,rho,odmax,YO=NULL)
{ 
  if(length(odmax)==1) { odmax<-rep(odmax,nrow(EZ)) }
  ZS<-simZ(EZ,rho)  
  
  diag(ZS)<- -Inf
  if(!is.null(YO)) { ZS[is.na(YO)]<- -Inf } 

  YS<-ZS*0 
  for(i in 1:nrow(EZ))
  {
    rs<-rank(ZS[i,])  -  (nrow(EZ)-odmax[i]) 
    YS[i,]<-rs*(rs>0)*(ZS[i,]>0) 
    YS[i,YS[i,]>0 ] <- match( YS[i,YS[i,]>0 ] ,sort(unique(YS[i,YS[i,]>0 ]))) 
  }
  diag(YS)<-NA 
  if(!is.null(YO)) { YS[is.na(YO)]<- NA } 

YS
}


}
