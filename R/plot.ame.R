plot.ame <-
function(x, ...)
{  
  fit<-x 
  require(amen) 

  gof<-1*(nrow(fit$GOF)>1) 
  par(mfrow=c(1+2*gof,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))

  mSABR<-apply(fit$SABR,2,median)
  matplot(fit$SABR,type="l",lty=1,ylab="SABR")
  abline(h=mSABR,col=1:length(mSABR) )

  if(ncol(fit$BETA)>0) 
  { 
    mBETA<-apply(fit$BETA,2,median)
    matplot(fit$BETA,type="l",lty=1,col=1:length(mBETA),ylab="BETA")
    abline(h=mBETA,col=1:length(mBETA) )
    abline(h=0,col="gray")
  }

  if(gof)
  {
    for(k in 1:4)
    {
      hist(fit$GOF[-1,k],xlim=range(fit$GOF[,k]),main="",prob=TRUE,
           xlab=colnames(fit$GOF)[k],col="lightblue",ylab="",yaxt="n")
           abline(v=fit$GOF[1,k],col="red")
    }
  }
}

