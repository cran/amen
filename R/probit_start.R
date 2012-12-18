probit_start <-
function(Y,X,rm.int=FALSE)
{
  # starting values for the probit srm
  X<-X[,,which(apply(X,3,function(x){var(c(x))})!=0),drop=FALSE]
  fit<-glm(c(1*(Y>0))~ apply(X,3,c),family=binomial(link=probit))
  E<-matrix(0,nrow(Y),ncol(Y)) ; E[!is.na(Y)]<-fit$res

  a<-apply(E,1,mean) ; b<-apply(E,2,mean)
  E<-E - outer(a,b,"+")
  rho<-cor(cbind(E[upper.tri(E)],t(E)[upper.tri(E)]))[1,2]
  beta<-fit$coef; if(rm.int){beta<-beta[-1] }
  list(beta=beta,a=a,b=b,rho=rho)
}
