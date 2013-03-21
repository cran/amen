lfcplot<-function(U,V,Y=NULL,ecol="lightblue",sncol="darkgreen",rncol="darkred",snlab=NULL,rnlab=NULL)
{
  uhat<-U[,1:2] ; vhat<-V[,1:2]

  vd<-vhat/sqrt(apply(vhat^2,1,sum))
  vm<-sqrt(sqrt(apply(vhat^2,1,sum)) )
  vm<-1.2*vm/max(vm)

  ud<-uhat/sqrt(apply(uhat^2,1,sum))
  um<-sqrt(sqrt(apply(uhat^2,1,sum)))
  um<-1.2*um/max(um)

  par(mar= .5+c(0,0,0,0),mgp=c(1.75,.75,0))
  plot(ud,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")

  edges<-which(Y!=0,arr.ind=TRUE)
  for(k in 1:nrow(edges))
  {
    i<-edges[k,1] ; j<-edges[k,2]
    segments(ud[i,1],ud[i,2],.8*vd[j,1],.8*vd[j,2],col=ecol)
  }

  if(is.null(snlab)){snlab<-as.character(1:nrow(Y))}
  if(is.null(rnlab)){rnlab<-as.character(1:nrow(Y))}

  if( is.character(snlab)){ text(ud[,1],ud[,2],snlab,cex=um,col=sncol) }
  if( is.character(rnlab)){ text(vd[,1]*.8,vd[,2]*.8,rnlab,cex=vm,col=rncol) }

  if(!is.character(snlab)) { points(ud[,1],ud[,2],pch=snlab,cex=um,col=sncol) }
  if(!is.character(rnlab)) { points(vd[,1]*.8,vd[,2]*.8,pch=rnlab,cex=vm,col=rncol) }
}




