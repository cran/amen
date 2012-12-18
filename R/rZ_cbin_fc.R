rZ_cbin_fc <-
function(Z,EZ,rho,Y,odmax,odobs)
{
  # simulates Z under the contraints 
  # (1)  Y[i,j] > Y[i,k]              => Z[i,j]>Z[i,k]  
  # (2)  Y[i,j]>0                     => Z[i,j]>0    
  # (3)  Y[i,j]=0 & odobs[i]<odmax[i] => Z[i,j]<0

  sz<-sqrt(1-rho^2)
  ut<-upper.tri(EZ)
  lt<-lower.tri(EZ)

  for(y in sample(0:1))
  {
    if(y==1)
    { 
      ub<- Inf
      lbm<-matrix(pmax(apply(Z-(Y!=0)*(Inf^(Y!=0)),1,max,na.rm=TRUE),0),
                  nrow(Z),nrow(Z))
    }

    if(y==0)
    { 
      lb<- -Inf
      ubm<-matrix(apply(Z+(Y!=1)*(Inf^(Y!=1)),1,min,na.rm=TRUE),nrow(Z), 
                  nrow(Z))
      ubm[ odobs<odmax ] <- 0 
    }

    up<- ut & Y==y 
    if(y==0)  { ub<-ubm[up] }  
    if(y==1)  { lb<-lbm[up] }
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    Z[up]<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))

    up<- lt & Y==y
    if(y==0)  { ub<-ubm[up] }
    if(y==1)  { lb<-lbm[up] }
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    Z[up]<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))
  }

  #diag(Z)<-rnorm(nrow(Z),diag(EZ),1)
  Z[is.na(Y)]<- rnorm(sum(is.na(Y)),EZ[is.na(Y)],1)
  Z
}

