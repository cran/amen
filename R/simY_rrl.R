simY_rrl <-
function(EZ,rho,odobs,YO=NULL)
{ 
  ZS<-simZ(EZ,rho)  
  diag(ZS)<- -Inf
  if(!is.null(YO)) { ZS[is.na(YO)]<- -Inf } 

  YS<-ZS*0  
  for(i in 1:nrow(EZ))
  {
    ri<-order( -ZS[i,] )[seq(1,odobs[i],length=odobs[i]) ]
    YS[i,ri]<- seq(odobs[i],1,length=odobs[i])
  }
  diag(YS)<-NA 
  if(!is.null(YO)) { YS[is.na(YO)]<- NA }
YS
}



