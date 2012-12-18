Xbeta <-
function(X,beta)
{
  XB<-matrix(0,nrow=dim(X)[1],ncol=dim(X)[2] )
  for(k in seq(1,length(beta),length=length(beta))){XB<-XB + beta[k]*X[,,k]}
  XB
}
