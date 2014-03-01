###
ame<-
function (Y,Xdyad=NULL, Xrow=NULL, Xcol=NULL, 
    rvar = !(model=="rrl") , cvar = TRUE, dcor = TRUE, R = 0,
    model="nrm",
    intercept=!is.element(model,c("rrl","ord")),
    odmax=rep(max(apply(Y>0,1,sum,na.rm=TRUE)),nrow(Y)),
    seed = 1, nscan = 50000, burn = 500, odens = 25,
    plot=TRUE, print = TRUE, gof=TRUE)
{ 
 
  ## 
  set.seed(seed)

  diag(Y) <- NA
  if(is.element(model,c("bin","cbin"))) { Y<-1*(Y>0) } 

  if(is.element(model,c("cbin","frn","rrl"))){odobs<-apply(Y>0,1,sum,na.rm=TRUE)}

  X<-design_array(Xrow,Xcol,Xdyad,intercept,nrow(Y)) 
  if(is.element(model,c("rrl","ord"))&any(apply(apply(X,c(1,3),var),2,sum)==0))
  {
    cat("WARNING: row effects are not estimable using this procedure ","\n")
  }
  ##

  ##
  if(length(odmax)==1) { odmax<-rep(odmax,nrow(Y)) }
  ## 
 
  ##
  if(is.element(model,c("frn","rrl")))
  {
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
  }
  ##

  ## starting Z values
  if(model=="nrm") { Z<-Y }
  if(model=="ord") { Z<-matrix(zscores(Y),nrow(Y),ncol(Y)) } 
  if(model=="rrl") { Z<-matrix(t(apply(Y,1,zscores)),nrow(Y),ncol(Y)) }  
  if(model=="bin")
  { 
    Z<-matrix(zscores(Y),nrow(Y),nrow(Y)) 
    z01<- .5* ( max(Z[Y==0],na.rm=TRUE) + min(Z[Y==1],na.rm=TRUE) ) 
    Z<-Z - z01
  } 

  if(is.element(model,c("cbin","frn")))
  {
    Z<-Y
    for(i in 1:nrow(Y))
    {
      yi<-Y[i,]
      zi<-zscores(yi)
      rnkd<-which( !is.na(yi) & yi>0 ) 
      if(length(rnkd)>0 && min(zi[rnkd])<0) { zi[rnkd]<-zi[rnkd] - min(zi[rnkd]) + 1e-3 }

      if(length(rnkd)<odmax[i]) 
      {
        urnkd<-which( !is.na(yi) & yi==0 ) 
        if(max(zi[urnkd])>0) { zi[urnkd]<-zi[urnkd] - max(zi[urnkd]) -1e-3 }
      }

    Z[i,]<-zi
    } 
  }

  Z0<-Z #"observed" value of Z 
  ZA<-outer(rowMeans(Z,na.rm=TRUE),colMeans(Z,na.rm=TRUE),"+")/2
  Z[is.na(Z)]<-ZA[is.na(Z)] 
  ##


  ## zeros for other starting values
  beta<-rep(0,dim(X)[3]) 
  a <-b<-rep(0,nrow(Y)) 
  s2 <-1 
  rho<-0
  Sab<-diag(c(rvar,cvar))
  U<-V<-matrix(0, nrow(Y), R)  ##


  ##  output items
  BETA <- matrix(nrow = 0, ncol = dim(X)[3])
  SABR<-matrix(nrow=0,ncol=5) 
  UVPS <- U %*% t(V) * 0 
  APS<-BPS<- rep(0,nrow(Y))  
  YPS<-matrix(0,nrow(Y),ncol(Y)) ; dimnames(YPS)<-dimnames(Y)
  GOF<-matrix(gofstats(Y),1,4)  
  rownames(GOF)<-"obs"
  colnames(GOF)<- c("sd.rowmean","sd.colmean","dyad.dep","triad.dep")

  colnames(BETA) <- dimnames(X)[[3]] 
  colnames(SABR) <- c("va", "cab", "vb", "rho", "ve")
  names(APS)<-names(BPS)<- rownames(U)<-rownames(V)<-rownames(Y) 

  have_coda<-suppressWarnings(try(require(coda,quietly = TRUE),silent=TRUE))

  ## marginal means and regression sums of squares
  Xr<-apply(X,c(1,3),sum)            # row sum
  Xc<-apply(X,c(2,3),sum)            # col sum
  mX<- apply(X,3,c)                  # design matrix
  mXt<-apply(aperm(X,c(2,1,3)),3,c)  # dyad-transposed design matrix
  XX<-t(mX)%*%mX                     # regression sums of squares
  XXt<-t(mX)%*%mXt                   # crossproduct sums of squares


  ## MCMC
  for (s in 1:(nscan + burn)) 
  { 

    ## update Z 
    EZ<-Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V)
    if(model=="nrm"){ Z<-rZ_nrm_fc(Z,EZ,rho,s2,Y) ; s2<-rs2_fc(Z-EZ,rho) }
    if(model=="bin"){ Z<-rZ_bin_fc(Z,EZ,rho,Y) }
    if(model=="ord"){ Z<-rZ_ord_fc(Z,EZ,rho,Y) }
    if(model=="cbin"){Z<-rZ_cbin_fc(Z,EZ,rho,Y,odmax,odobs)}
    if(model=="frn"){ Z<-rZ_frn_fc(Z,EZ,rho,Y,YL,odmax,odobs)}
    if(model=="rrl"){ Z<-rZ_rrl_fc(Z,EZ,rho,Y,YL)} 

    ## update beta, a b
    tmp <- rbeta_ab_fc(Z-U%*%t(V), Sab, rho, X, mX, mXt,XX, XXt, Xr, Xc, s2)
    beta <- tmp$beta
    a <- tmp$a * rvar
    b <- tmp$b * cvar 

    ## update Sab 
    if(rvar & cvar)
    {
      if(is.element(model,c("nrm","ord")))
      { 
        Sab<-solve(rwish(solve(diag(2)+crossprod(cbind(a, b))),3+nrow(Z))) 
      }

      if(model=="bin")
      {
        tmp<-raSab_bin_fc(Z,Y,a,b,Sab) ; Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      }

      if(model=="cbin")
      {
        tmp<-raSab_cbin_fc(Z,Y,a,b,Sab,odmax,odobs)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      }
      
      if(model=="frn")
      { 
        tmp<-raSab_frn_fc(Z,Y,YL,a,b,Sab,odmax,odobs) 
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a 
      }
    }
    if (rvar & !cvar) 
    {
      Sab[1, 1] <- 1/rgamma(1, (1 + nrow(Y))/2, (1 + sum(a^2))/2)
    }

    if (!rvar & cvar) 
    {
      Sab[2, 2] <- 1/rgamma(1, (1 + nrow(Y))/2, (1 + sum(b^2))/2)
    }

    ## update rho
    if (dcor) 
    {
      rho<-rrho_mh(Z-(Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V)), rho,s2)
    }

    ## update U,V
    if (R > 0) 
    {
      UV <- rUV_fc(Z - (Xbeta(X, beta) + outer(a, b, "+")), U, V, rho, s2)
      U <- UV$U
      V <- UV$V
    } 

    ## output
    if(s%%odens==0&s<=burn){cat(round(100*s/burn,2)," pct burnin complete \n")}

    if(s%%odens==0 & s>burn) 
    { 
      BETA <- rbind(BETA, beta)
      SABR <- rbind(SABR, c(Sab[upper.tri(Sab, diag = T)], rho,s2))
      UVPS <- UVPS + U %*% t(V)
      APS <- APS + a
      BPS <- BPS + b 

      # simulate from posterior predictive 
      EZ<-Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V)
      if(model=="bin") { Ys<-simY_bin(EZ,rho) }
      if(model=="cbin"){ Ys<-1*(simY_frn(EZ,rho,odmax,YO=Y)>0) }
      if(model=="frn") { Ys<-simY_frn(EZ,rho,odmax,YO=Y) }
      if(model=="rrl") { Ys<-simY_rrl(EZ,rho,odobs,YO=Y ) }
      if(model=="nrm") { Ys<-simY_nrm(EZ,rho,s2) }
      if(model=="ord") { Ys<-simY_ord(EZ,rho,Y) }
      YPS<-YPS+Ys

      if(gof){ Ys[is.na(Y)]<-NA ; GOF<-rbind(GOF,gofstats(Ys)) }

      if (print) 
      {
        cat(s,round(apply(BETA,2,mean),2),":",round(apply(SABR,2,mean),2),"\n")
        if (have_coda & nrow(SABR) > 3 & length(beta)>0) 
        {
          cat(round(effectiveSize(BETA)), "\n")
        }
      }
      if(plot) 
      {
        par(mfrow=c(1+2*gof,2),mar=c(3,3,1,1),mgp=c(1.75,0.75,0))
        mSABR <- apply(SABR, 2, median)
        matplot(SABR, type = "l", lty = 1)
        abline(h = mSABR, col = 1:length(mSABR)) 

        if(length(beta)>0) 
        {
          mBETA <- apply(BETA, 2, median)
          matplot(BETA, type = "l", lty = 1, col = 1:length(mBETA))
          abline(h = mBETA, col = 1:length(mBETA))
          abline(h = 0, col = "gray") 
        } 
        if(gof)
        {
          for(k in 1:4)
          {
            hist(GOF[-1,k],xlim=range(GOF[,k]),main="",prob=TRUE,
                 xlab=colnames(GOF)[k],col="lightblue",ylab="",yaxt="n")  
            abline(v=GOF[1,k],col="red") 
          }
        } 
      }
    }
  }   

  ## output
  UVPM<-UVPS/nrow(SABR)
  UDV<-svd(UVPM)
  U<-UDV$u[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
  V<-UDV$v[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
  APM<-APS/nrow(SABR)
  BPM<-BPS/nrow(SABR) 
  names(APM)<-names(BPM)<-rownames(U)<-rownames(V)<-rownames(Y) 
  EZ<-Xbeta(X,apply(BETA,2,mean)) + outer(APM,BPM,"+")+UVPM
  dimnames(EZ)<-dimnames(Y) 
  YPM<-YPS/nrow(SABR)

  fit <- list(BETA=BETA,SABR=SABR,APM=APM,BPM=BPM,U=U,V=V,UVPM=UVPM,EZ=EZ,
             YPM=YPM,GOF=GOF)
  class(fit) <- "ame"
  fit
}



