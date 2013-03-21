probit_start<-
function(Y,X)
{ 
  if(dim(X)[3]>0) 
  {
  fit<-glm(c(1*(Y>0))~ -1 + apply(X,3,c),family=binomial(link=probit)) 
  beta<-fit$coef
  res<-fit$res 
  }

  if(dim(X)[3]==0)
  { 
  beta<-numeric(0) 
  z<- rank(Y,na.last="keep")[!is.na(Y)]  
  res<-qnorm( z/(length(z)+1) ) 
  }


  E <- matrix(NA, nrow(Y), ncol(Y))
  E[!is.na(Y)] <- res
  a <- apply(E, 1, mean,na.rm=TRUE)
  b <- apply(E, 2, mean,na.rm=TRUE)

  E<-E - outer(a,b,"+")
  rho<-cor(cbind(E[upper.tri(E)],t(E)[upper.tri(E)]),use="complete.obs")[1,2]
  list(beta=beta,a=a,b=b,rho=rho)
}
