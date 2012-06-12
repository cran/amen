rZ_bin_fc <-
function(Z,EZ,rho,Y)
{
  sz<-sqrt(1-rho^2)
  ut<-upper.tri(EZ)
  lt<-lower.tri(EZ)
  
  for(y in 0:1)
  {
    lb<-c(-Inf,0)[y+1] ; ub<-c(0,Inf)[y+1]  

    up<- ut & Y==y
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    Z[up]<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))

    up<- lt & Y==y
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    Z[up]<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))
  }

  diag(Z)<-rnorm(nrow(Z),diag(EZ),1)
  Z
}
