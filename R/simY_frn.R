simY_frn <-
function(EZ,rho,odmax,YO=NULL)
{ 
  if(length(odmax)==1) { odmax<-rep(odmax,nrow(EZ)) }
  ZS<-simZ(EZ,rho)  
  
  diag(ZS)<- -Inf
  if(!is.null(YO)) { ZS[is.na(YO)]<- -Inf } 

  YS<-ZS*0 
  for(i in 1:nrow(EZ))
  {
    rs<-rank(ZS[i,])  -  (nrow(EZ)-odmax[i]) 
    YS[i,]<-rs*(rs>0)*(ZS[i,]>0) 
    YS[i,YS[i,]>0 ] <- match( YS[i,YS[i,]>0 ] ,sort(unique(YS[i,YS[i,]>0 ]))) 
  }
  diag(YS)<-NA 
  if(!is.null(YO)) { YS[is.na(YO)]<- NA } 

YS
}