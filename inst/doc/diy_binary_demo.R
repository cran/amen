## ----message=FALSE------------------------------------------------------------
library(amen)

Y<-1*(IR90s$dyadvars[,,"conflicts"] >0 )

netplot(Y,plot.iso=FALSE,plotnames=TRUE)

## -----------------------------------------------------------------------------
Xn<-log(IR90s$nodevars[,2])
Xd<-log(IR90s$dyadvars[,,3]+1) 

## -----------------------------------------------------------------------------
X<-design_array(Xrow=Xn,Xcol=Xn,Xdyad=Xd,n=nrow(Y)) 

## -----------------------------------------------------------------------------
dim(X) 

dimnames(X)[[3]] 

X[1:3,1:3,1]

X[1:3,1:3,2]

X[1:3,1:3,3]

X[1:3,1:3,4]

## -----------------------------------------------------------------------------
Sab<-diag(2) 

rho<-0 

R<-2 

Suv<-diag(2*R)

U<-V<-matrix(0,nrow(Y),R) 

# just need Z[i,j]<0 if Y[i,j]=0, Z[i,j]>0 if Y[i,j]=1  
Z<-matrix(zscores(Y,ties.method="random"),nrow(Y),nrow(Y))
Z<-Z-max(Z[!is.na(Y) & Y==0 ] ) 
diag(Z)<-0

## -----------------------------------------------------------------------------
ybar<-mean(Y,na.rm=TRUE) ; mu<-qnorm(ybar)
E<- (Y - ybar)/dnorm(qnorm(ybar)) ; diag(E)<-0
a<-rowMeans(E,na.rm=TRUE)  ; a[is.na(a)]<-0
b<-colMeans(E,na.rm=TRUE)  ; b[is.na(b)]<-0

vscale<-mean(diag(cov(cbind(a,b))))
PHAT<-pnorm(mu+outer(a,b,"+"))
vdfmlt<-.25/mean(PHAT*(1-PHAT))

Sab0<-diag(2)*vscale 
eta0<-round(4*vdfmlt)       

Suv0<-diag(2*R)*vscale 
kappa0<-round((2*R+2)*vdfmlt)  

## ----results="hide",cache=TRUE,fig.keep='last'--------------------------------
BETA<-NULL ; RHO<-NULL ; UV<-matrix(0,nrow(Y),nrow(Y))

for(s in 1:10000)
{ 
  # update beta, a and b  
  tmp<-rbeta_ab_fc(Z,Sab,rho,X,offset=U%*%t(V))
  beta<-tmp$beta ; a<-tmp$a ; b<-tmp$b

  # update UV 
  tmp<-rUV_fc(Z,U,V,Suv,rho,offset=Xbeta(X,beta) + outer(a,b,"+")) 
  U<-tmp$U ; V<-tmp$V

  # update Suv
  Suv<-rSuv_fc(U,V,Suv0=Suv0,kappa0=kappa0) 

  # update Sab
  Sab<-rSab_fc(a,b,Sab0=Sab0,eta0=eta0)

  # update rho 
  rho<-rrho_mh(Z,rho,offset=Xbeta(X,beta)+outer(a,b,"+")+U%*%t(V)) 

  # update Z 
  Z<-rZ_bin_fc(Z, Xbeta(X,beta)+outer(a,b,"+")+U%*%t(V), rho, Y)  

  # periodically save some output   
  if(s%%25==0 & s>500) 
  { 
    cat(round(100*s/10000),"% complete\n") 
    BETA<-rbind(BETA,beta) 
    RHO<-c(RHO,rho) 
    UV<-UV+U%*%t(V)  
    matplot(BETA,type="l") 
  }
}

## ----results='hide',fig.keep='last',cache=TRUE--------------------------------
fit<-ame(Y,Xrow=Xn,Xcol=Xn,Xdyad=Xd,R=2,family="bin")

## ----echo=c(-1)---------------------------------------------------------------
par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
for(j in 1:4)
{
  plot(density(BETA[,j] ),main="",xlab="")
  lines(density(fit$BETA[,j] ),col="green") 
}

