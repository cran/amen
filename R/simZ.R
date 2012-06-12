simZ <-
function(EZ,rho,s2=1)
{
  w1<-sqrt((1 + sqrt(1-rho^2))/2)
  w2<-sign(rho)*sqrt(1-w1^2)
  EC<-matrix(rnorm(length(EZ)),nrow(EZ),nrow(EZ))
  EC<- sqrt(s2)*( w1*EC + w2*t(EC) )
  EZ+EC    
}
