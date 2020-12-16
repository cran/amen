## -----------------------------------------------------------------------------
library(amen)

Y<-sheep$dom 
age<-sheep$age 

Y[1:8,1:8] 

age[1:8]

## ----echo=c(-1)---------------------------------------------------------------
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
boxplot( apply(Y,1,mean,na.rm=TRUE) ~ age )  
boxplot( apply(Y,2,mean,na.rm=TRUE) ~ age )  

## -----------------------------------------------------------------------------
x<-sheep$age - mean(sheep$age)

X<-design_array(Xrow=x,Xcol=x,Xdyad=outer(x,x)) 

## -----------------------------------------------------------------------------
Z<-log(Y+1) ; diag(Z)<-0 

Sab<-diag(2) 

R<-1 ; Suv<-diag(2*R) ; U<-V<-matrix(0,nrow(Y),R) 

s2<-1 ; rho<-0

## -----------------------------------------------------------------------------
for(s in 1:10)
{ 
  # update beta, a and b
  tmp<-rbeta_ab_fc(Z,Sab,rho,X,s2,offset=U%*%t(V))
  beta<-tmp$beta ; a<-tmp$a ; b<-tmp$b

  # update UV
  tmp<-rUV_fc(Z,U,V,Suv,rho,offset=Xbeta(X,beta)+outer(a,b,"+"))
  U<-tmp$U ; V<-tmp$V

  # update Sab
  Sab<-rSab_fc(a,b)

  # update Suv
  Suv<-rSuv_fc(U,V)

  # update s2
  s2<-rs2_fc(Z,rho,offset=Xbeta(X,beta)+outer(a,b,"+")+U%*%t(V))

  # update rho
  rho<-rrho_mh(Z,rho,s2,offset=Xbeta(X,beta)+outer(a,b,"+")+U%*%t(V))
}

## -----------------------------------------------------------------------------
# propose candidate Z
Zp<-Z+matrix(rnorm(nrow(Y)^2),nrow(Y),nrow(Y))*sqrt(s2)

# compute acceptance ratio
EZ<-Xbeta(X,beta)+outer(a,b,"+")+U%*%t(V)
lr<-ldZgbme(Zp,Y,function(y,z){ dpois(y,exp(z),log=TRUE) },EZ,rho,s2) -
    ldZgbme(Z, Y,function(y,z){ dpois(y,exp(z),log=TRUE) },EZ,rho,s2)

# simulate symmetric matrix of (log) uniform rvs
lh<-matrix(log(runif(nrow(Y)^2)),nrow(Y),nrow(Y))
lh[lower.tri(lh)]<-t(lh)[lower.tri(lh)]

# update dyads for which lr>lh, and keep track of acceptances
Z[lr>lh]<-Zp[lr>lh]

## ----results="hide",cache=TRUE,fig.keep='last'--------------------------------
#### Starting values
Z<-log(Y+1) ; diag(Z)<-0

Sab<-diag(2)

R<-1 ; Suv<-diag(2*R) ; U<-V<-matrix(0,nrow(Y),R) 

s2<-1 ; rho<-0

#### Parameter values to be saved
BETA<-VE<-LL<-NULL
ACC<-matrix(0,nrow(Y),nrow(Y))

#### MCMC 
set.seed(1)
for(s in 1:10000)
{ 
  ## Update AMEN parameters

  # update beta, a and b
  tmp<-rbeta_ab_fc(Z,Sab,rho,X,s2,offset=U%*%t(V))
  beta<-tmp$beta ; a<-tmp$a ; b<-tmp$b

  # update UV
  tmp<-rUV_fc(Z,U,V,Suv,rho,offset=Xbeta(X,beta)+outer(a,b,"+"))
  U<-tmp$U ; V<-tmp$V

  # update Sab
  Sab<-rSab_fc(a,b)

  # update Suv
  Suv<-rSuv_fc(U,V)

  # update s2
  s2<-rs2_fc(Z,rho,offset=Xbeta(X,beta)+outer(a,b,"+")+U%*%t(V))

  # update rho
  rho<-rrho_mh(Z,rho,s2,offset=Xbeta(X,beta)+outer(a,b,"+")+U%*%t(V))



  ## Update Z 
 
  # propose candidate Z 
  Zp<-Z+matrix(rnorm(nrow(Y)^2),nrow(Y),nrow(Y))*sqrt(s2)   

  # compute acceptance ratio 
  EZ<-Xbeta(X,beta)+outer(a,b,"+")+U%*%t(V)
  lr<-ldZgbme(Zp,Y,function(y,z){ dpois(y,exp(z),log=TRUE) },EZ,rho,s2) -
      ldZgbme(Z, Y,function(y,z){ dpois(y,exp(z),log=TRUE) },EZ,rho,s2) 

  # simulate symmetric matrix of (log) uniform rvs
  lh<-matrix(log(runif(nrow(Y)^2)),nrow(Y),nrow(Y))
  lh[lower.tri(lh)]<-t(lh)[lower.tri(lh)]  

  # update dyads for which lr>lh, and keep track of acceptances
  Z[lr>lh]<-Zp[lr>lh] 
  ACC[lr>lh]<-ACC[lr>lh]+1  



  ## Output
  if(s%%25==0) 
  {  
    cat(s,range(ACC/s),"\n") 
    BETA<-rbind(BETA,beta) 
    VE<-rbind(VE,c(s2,rho))    

    par(mfrow=c(1,3),mar=c(3,3,1,1),mgp=c(1.75,.75,0)) 
    matplot(BETA,type="l") 
    LL<-c(LL,sum(dpois(Y,exp(Z),log=TRUE),na.rm=TRUE)) ; plot(LL,type="l")
    hist(ACC/s) 
  }
}

## -----------------------------------------------------------------------------
apply(BETA,2,mean)

apply(BETA,2,sd)

## -----------------------------------------------------------------------------
apply(VE,2,mean) 

apply(VE,2,sd) 

