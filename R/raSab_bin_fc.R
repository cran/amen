raSab_bin_fc <-
function(Z,Y,a,b,Sab,SS=round(sqrt(nrow(Z))))
{
  E<-Z-a%*%t(rep(1,nrow(Z))) 
  MEL<-MEU<- -E;MEL[Y==0]<- -Inf;MEU[Y==1]<-Inf;diag(MEU)<-Inf;diag(MEL)<- -Inf
  lba<-apply(MEL,1,max) 
  lba[is.na(lba)]<- -Inf
  uba<-apply(MEU,1,min) 
  uba[is.na(uba)]<- Inf

  for(ss in 1:SS)
  {
    ea<-b*Sab[1,2]/Sab[2,2]
    sa<-sqrt(Sab[1,1]-Sab[1,2]^2/Sab[2,2])
    a<-ea+sa*qnorm(runif(nrow(Z),pnorm((lba-ea)/sa),pnorm((uba-ea)/sa)))
    Sab<-solve(rwish(solve(diag(2)+crossprod(cbind(a,b))),3+nrow(Z)))
  }
  list(Z=E+a%*%t(rep(1,nrow(Z))),a=a,Sab=Sab)
}
