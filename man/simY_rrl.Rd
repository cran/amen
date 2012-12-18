\name{simY_rrl}
\alias{simY_rrl}
\title{
Simulate an relational matrix based on a relative rank nomination scheme
}
\description{
Simulate an relational matrix based on a relative rank nomination scheme
}
\usage{
simY_rrl(EZ, rho, odobs, YO)
}
\arguments{
  \item{EZ}{
a square matrix giving the expected value of the latent Z matrix
}
  \item{rho}{
dyadic correlation
}
  \item{odobs}{
a scalar or vector giving the observed number of nominations for each node
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
function(EZ,rho,odobs,YO=NULL)
{
  ZS<-simZ(EZ,rho)
  diag(ZS)<- -Inf
  if(!is.null(YO)) { ZS[is.na(YO)]<- -Inf }

  YS<-ZS*0
  for(i in 1:nrow(EZ))
  {
    ri<-order( -ZS[i,] )[seq(1,odobs[i],length=odobs[i]) ]
    YS[i,ri]<- seq(odobs[i],1,length=odobs[i])
  }
  diag(YS)<-NA
  if(!is.null(YO)) { YS[is.na(YO)]<- NA }
YS
}
}
