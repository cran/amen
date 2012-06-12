rZ_frn_fc <-
function(Z,EZ,rho,Y,YL,odmax,odobs)
{
  # simulates Z under the contraints
  # (1)  Y[i,j]>Y[i,k]                => Z[i,j]>Z[i,k]  (same as rank)
  # (2)  Y[i,j]>0                     => Z[i,j]>0    
  # (3)  Y[i,j]=0 & odobs[i]<odmax[i] => Z[i,j]<0

  sz<-sqrt(1-rho^2)
  ut<-upper.tri(Z)
  lt<-lower.tri(Z)
  rws<-outer(1:nrow(Z),rep(1,nrow(Z)))

  for(y in sample(c(0:ncol(YL))) )
  {
    if(y<2)
    {
      if(y==0){lbm<- rep(-Inf,nrow(Z))}
      if(y==1){lbm<-pmax(0,apply(Z - (Y!=0)*(Inf^(Y!=0)),1,max,na.rm=TRUE)) }
    }
    if(y>=2) {lbm<-Z[cbind(1:nrow(Z),YL[,y-1])] }
   
    if(y<ncol(YL))  { ubm<- Z[ cbind(1:nrow(Z), YL[,y+1] )]  }
    if(y==0) { ubm[odobs<odmax]<- 0 }
    if(y==ncol(YL)) { ubm<- rep(Inf,nrow(Z)) }
    ubm[is.na(ubm)]<-Inf ; lbm[is.na(lbm)]<- -Inf

    for(k in sample(1:2)) 
    { 
    if(k==1) { 
    up<- ut & Y==y
    rwb<-rws[up]
    lb<-lbm[rwb] ; ub<-ubm[rwb]
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    Z[up]<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))
              }
    if(k==2)  { 
    up<- lt & Y==y
    rwb<-rws[up]
    lb<-lbm[rwb] ; ub<-ubm[rwb]
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    Z[up]<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz))) 
              } 
    } 
  }

  diag(Z)<-rnorm(nrow(Z),diag(EZ),1)
  Z
}
