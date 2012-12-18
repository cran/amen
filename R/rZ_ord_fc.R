rZ_ord_fc <-
function(Z,EZ,rho,W,FY)
{
  # simulates Z under the contraints
  # (1)  Y[i,j]>Y[k,l] => Z[i,j]>Z[k,l]

  sz<-sqrt(1-rho^2)
  ut<-upper.tri(Z)
  lt<-lower.tri(Z)

  W[is.na(W)]<- -1

  for(w in sample(c(-1,1:length(FY)) ))
  {
    lb<-suppressWarnings(max(Z[!is.na(W) & W==w-1],na.rm=TRUE))
    ub<-suppressWarnings(min(Z[!is.na(W) & W==w+1],na.rm=TRUE))

    up<- ut & W==w
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    Z[up]<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))

    up<- lt & W==w
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    Z[up]<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))
  }
  diag(Z)<-rnorm(nrow(Z),diag(EZ),1)
  Z
}
