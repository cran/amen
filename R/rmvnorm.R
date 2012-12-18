rmvnorm <-
function(n,mu,Sigma,Sigma.chol=chol(Sigma))
{
  # sample from a matrix normal distribution
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%Sigma.chol) +c(mu))
}
