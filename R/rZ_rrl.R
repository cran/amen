#####
rZ_rrl_fc<-function(Z,EZ,rho,Y,YL)
{
  # simulates Z under the contraints
  # (1)  Y[i,j]>Y[i,k] => Z[i,j]>Z[i,k]

  sz<-sqrt(1-rho^2)
  ut<-upper.tri(Z)
  lt<-lower.tri(Z)
  rws<-outer(1:nrow(Z),rep(1,nrow(Z)))

  Y[is.na(Y)]<- -1

  for(y in c((-1):ncol(YL)) )
  {
    if(y<2)
    {
      if(y<=0){lbm<- rep(-Inf,nrow(Z))}
      if(y==1){lbm<- apply(Z - (Y!=0)*(Inf^(Y!=0)),1,max,na.rm=TRUE)    }
    }
    if(y>=2) {lbm<-Z[cbind(1:nrow(Z),YL[,y-1])] }

    if(y== -1) { ubm<-rep(Inf,nrow(Z))}
    if(y<ncol(YL) & y>=0)  { ubm<- Z[ cbind(1:nrow(Z), YL[,y+1] )]  }
    if(y==ncol(YL)) { ubm<- rep(Inf,nrow(Z)) }
    ubm[is.na(ubm)]<-Inf ; lbm[is.na(lbm)]<- -Inf

    up<- ut & Y==y
    rwb<-rws[up]
    lb<-lbm[rwb] ; ub<-ubm[rwb]
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    Z[up]<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))

    up<- lt & Y==y
    rwb<-rws[up]
    lb<-lbm[rwb] ; ub<-ubm[rwb]
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    Z[up]<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))
  }

  diag(Z)<-rnorm(nrow(Z),diag(EZ),1)
  Z
}
#####




