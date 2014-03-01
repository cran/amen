rZ_nrm_fc<-function(Z,EZ,rho,s2,Y)
{
  ZS<-simY_nrm(EZ,rho,s2)
  diag(ZS)<-rnorm(nrow(Y),diag(EZ),sqrt(s2))
  Z[is.na(Y)]<-ZS[is.na(Y)]
  Z
}

