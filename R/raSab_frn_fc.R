raSab_frn_fc <-
function(Z,Y,YL,a,b,Sab,odobs,odmax,SS=round(sqrt(nrow(Z))))
{
  E<-Z-a%*%t(rep(1,nrow(Z)))
  lba<- -E[ cbind(1:nrow(Z), YL[,1]) ] ; lba[is.na(lba)]<- -Inf 
  uba<-  -apply(E - (Y!=0)*(Inf^(Y!=0)),1,max,na.rm=TRUE)
  uba[odobs==odmax]<- Inf 

  for(ss in 1:SS)
  {
    ea<-b*Sab[1,2]/Sab[2,2] 
    sa<-sqrt(Sab[1,1]-Sab[1,2]^2/Sab[2,2])
    a<-ea+sa*qnorm(runif(nrow(Z),pnorm((lba-ea)/sa),pnorm((uba-ea)/sa))) 
    Sab<-solve(rwish(solve(diag(2)+crossprod(cbind(a,b))),3+nrow(Z)))
  }
  list(Z=E+a%*%t(rep(1,nrow(Z))),a=a,Sab=Sab)
}
