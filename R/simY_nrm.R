simY_nrm <-
function(EY,rho,s2) 
{
  YS<-simZ(EY,rho,s2) 
  diag(YS)<-NA 
  YS
}
