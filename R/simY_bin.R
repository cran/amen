simY_bin <-
function(EZ,rho)
{
  ZS<-simZ(EZ,rho) 
  YS<-1*(ZS>0) ; diag(YS)<-NA
YS
}
