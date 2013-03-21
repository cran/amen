ame_rrl <-
function(Y,X=NULL,cvar=TRUE,dcor=TRUE,R=0,
             seed=1,nscan=5e4,burn=5e2,odens=25,plot=TRUE,print=TRUE)
{
set.seed(seed)


diag(Y)<-NA

if(is.null(X)) {  X<-array(dim=c(nrow(Y),nrow(Y),0)) }
if( length(dim(X))==2 ) { X<-array(X,dim=c(dim(X),1)) }
if(  any(apply( apply(X,c(1,3),var),2,sum )==0 ) )
{cat("WARNING: row effects are not estimable using this procedure ","\n") }

## marginal means and regression sums of squares
Xr<-apply(X,c(1,3),sum)            # row sum
Xc<-apply(X,c(2,3),sum)            # col sum
mX<- apply(X,3,c)                  # design matrix
mXt<-apply(aperm(X,c(2,1,3)),3,c)  # dyad-transposed design matrix
XX<-t(mX)%*%mX                     # regression sums of squares
XXt<-t(mX)%*%mXt                   # crossproduct sums of squares

## list of ordered friends: higher the number -> better the friend
ymx<-max(apply(1*(Y>0),1,sum,na.rm=TRUE))
YL<-NULL
warn<-FALSE
for(i in 1:nrow(Y))
{
  yi<-Y[i,] ; rnkd<-which( !is.na(yi)&yi>0 )
  if(length(yi[rnkd])>length(unique(yi[rnkd]))){warn<-TRUE}
  yi[rnkd]<-rank(yi[rnkd],ties.method="random")
  Y[i,]<-yi
  YL<-rbind(YL, match(1:ymx,yi))
}

if(warn){cat("WARNING: Random reordering applied to break ties in ranks\n")}


#################### model and starting values 
fit<-probit_start(Y,X)
beta<-fit$beta ; rho<- .99*fit$rho*dcor ; 
Sab<-matrix(c(0,0,0,(1+var(fit$b))*cvar),2,2)
Z<-Y
for(i in 1:nrow(Y)) { Z[i,-i]<-qnorm(rank(Y[i,-i])/nrow(Y)) }
diag(Z)<-apply(Z,1,max,na.rm=TRUE)
U<-V<-matrix(0,nrow(Y),R)
a<-rep(0,nrow(Y))
####################

## MCMC setup
tr.obs<-t_recip(1*(Y>0))   # reciprocation
tt.obs<-t_trans(1*(Y>0))   # transitivity
td.obs<-t_degree(1*(Y>0))  # degree distributions
odobs<-apply(Y>0,1,sum,na.rm=TRUE) # obs outdegrees
idobs<-apply(Y>0,2,sum,na.rm=TRUE) # obs indegrees

TT<-TR<-TID<-TOD<-SABR<-NULL
BETA<-matrix(nrow=0,ncol=dim(X)[3]) ; colnames(BETA)<-dimnames(X)[[3]]
UVPS<-U%*%t(V)*0
APS<-BPS<-rep(0,nrow(Y))
have.coda<-suppressWarnings(try(require(coda,quietly=TRUE),silent=TRUE))

## MCMC
for(s in 1:(nscan+burn)) 
{

  ## update beta,a,b
  tmp<-rbeta_ab_fc(Z-U%*%t(V),Sab,rho,X,mX,mXt,XX,XXt,Xr,Xc)
  beta<-tmp$beta ; b<-tmp$b*cvar
  ##

  ## update Z
  Z<-rZ_rrl_fc(Z,Xbeta(X,beta) + outer(a,b,"+")+U%*%t(V),rho,Y,YL)
  ##

  ## update covariance model
  if(cvar){ Sab[2,2]<-1/rgamma(1,(1+nrow(Y))/2,(1+sum(b^2))/2) }
  if(dcor){rho<-rrho_mh(Z-(Xbeta(X,beta)+outer(a,b,"+")+U%*%t(V)),rho)}
  ##

  ## update multiplicative effects 
  if(R>0){UV<-rUV_fc(Z-(Xbeta(X,beta)+outer(a,b,"+")),U,V,rho);U<-UV$U;V<-UV$V}
  ##

  ## output
  if(s%%odens==0 & s<=burn){cat(round(100*s/burn,2)," pct burnin complete \n")}
  if(s%%odens==0 & s>burn)
  { 

    ## save current parameter values
    BETA<-rbind(BETA, beta)
    SABR<-rbind(SABR,c( Sab[upper.tri(Sab,diag=T)],rho))
    UVPS<-UVPS+U%*%t(V)
    APS<-APS+a ; BPS<-BPS+b 
 
    ## simulate gof stats
    YS<-simY_rrl(Xbeta(X,beta)+outer(a,b,"+")+U%*%t(V),rho,odobs,YO=Y)       
    TR<-c(TR,t_recip(1*(YS>0)))
    TT<-c(TT,t_trans(1*(YS>0)))
    td<-t_degree(1*(YS>0)) ; TOD<- rbind(TOD,td$od) ; TID<- rbind(TID,td$id) 

    if(print)
    {
      cat(s,round(apply(BETA,2,mean),2),":",round(apply(SABR,2,mean),2),"\n")
      if(have.coda & length(TR)>3 & length(beta)>0) {cat(round(effectiveSize(BETA)),"\n") } 
    } 

    ## plots
    if(plot)
    {
      par(mfrow=c(4,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
      hist(a,main="",col="lightblue",prob=TRUE)
      hist(b,main="",col="lightblue",prob=TRUE)

      
      mSABR<-apply(SABR,2,median)
      matplot(SABR,type="l",lty=1)
      abline(h=mSABR,col=1:length(mSABR) )

      if(length(beta)>0)
      {
      mBETA<-apply(BETA,2,median)
      matplot(BETA,type="l",lty=1,col=1:length(mBETA))
      abline(h=mBETA,col=1:length(mBETA) )
      abline(h=0,col="gray")
      }
 
      mod<-max(odobs) ; mid<-max(idobs)
      qod<-apply(TOD,2,quantile,prob=c(.975,.5,.25))
      qid<-apply(TID,2,quantile,prob=c(.975,.5,.25))
      odd<-rbind(qod,td.obs$od)[,0:(mod+1)]
      plot(c(0,mod),range(odd),type="n")
      for(k in 1:4) { lines(0:mod,odd[k,],col=c("gray","black")[1+(k==4)]) }
      idd<-rbind(qid,td.obs$id)[,0:(mid+1)]
      plot(c(0,mid),range(idd),type="n")
      for(k in 1:4) { lines(0:mid,idd[k,],col=c("gray","black")[1+(k==4)]) }

      plot(TR,ylim=range(c(TR,tr.obs)),type="l") ; abline(h=tr.obs)
      plot(TT,ylim=range(c(TT,tt.obs)),type="l") ; abline(h=tt.obs)
    }   
  }
}


colnames(SABR)<-c("va","cab","vb","rho") 

###
UVPM<-UVPS/length(TT)
UDV<-svd(UVPM)
U<-UDV$u[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)]) 
V<-UDV$v[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)])
A<-APS/length(TT)
B<-BPS/length(TT) 
rownames(A)<-rownames(B)<-rownames(U)<-rownames(V)<-rownames(Y)
EZ<-Xbeta(X,apply(BETA,2,mean)) + outer(A,B,"+")+U%*%t(V)
dimnames(EZ)<-dimnames(Y) 
###

fit<-list(BETA=BETA,SABR=SABR,A=A,B=B,U=U,V=V,EZ=EZ,
          TT=TT,TR=TR,TID=TID,TOD=TOD, 
          tt=tt.obs,tr=tr.obs,td=td.obs) 
class(fit)<-"ame"
fit
}
