rwish <-
function(S0,nu=dim(S0)[1]+1)
{
  # sample from a Wishart distribution 
  # with expected value nu*S0 
  sS0<-chol(S0)
  Z<-matrix(rnorm(nu*dim(S0)[1]),nu,dim(S0)[1])%*%sS0
  t(Z)%*%Z
}
