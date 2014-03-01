gofstats<-function(Y)
{
  sd.rowmean<-sd(rowMeans(Y,na.rm=TRUE) ,na.rm=TRUE) 
  sd.colmean<-sd(colMeans(Y,na.rm=TRUE) ,na.rm=TRUE)
  
  dyad.dep<- cor( c(Y),c(t(Y)) , use="complete.obs")
 
  E<-Y-mean(Y,na.rm=TRUE) ;  D<-1*(!is.na(E)) ; E[is.na(E)]<-0
  triad.dep<- sum(diag(E%*%E%*%E))/( sum(diag(D%*%D%*%D)) * sd(c(Y),na.rm=TRUE)^3)

  gof<-c(sd.rowmean,sd.colmean, dyad.dep , triad.dep )
  names(gof)<-c("sd.rowmean","sd.colmean","dyad.dep","triad.dep")
  gof
}



