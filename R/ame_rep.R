#' AME model fitting routine for replicated relational data
#' 
#' An MCMC routine providing a fit to an additive and multiplicative effects
#' (AME) regression model to replicated relational data of
#' various types. 
#' 
#' This command provides posterior inference for parameters in AME models of
#' independent replicated relational data, assuming one of six possible data
#' types/models:
#' 
#' "nrm": A normal AME model.
#' 
#' "bin": A binary probit AME model.
#' 
#' "ord": An ordinal probit AME model. An intercept is not identifiable in this
#' model.
#' 
#' "cbin": An AME model for censored binary data.  The value of 'odmax'
#' specifies the maximum number of links each row may have.
#' 
#' "frn": An AME model for fixed rank nomination networks. A higher value of
#' the rank indicates a stronger relationship. The value of 'odmax' specifies
#' the maximum number of links each row may have.
#' 
#' "rrl": An AME model based on the row ranks. This is appropriate if the
#' relationships across rows are not directly comparable in terms of scale. An
#' intercept, row random effects and row regression effects are not estimable
#' for this model.
#' 
#' @usage ame_rep(Y,Xdyad=NULL, Xrow=NULL, Xcol=NULL, rvar = !(model=="rrl")
#' , cvar = TRUE, dcor = TRUE, R = 0, model="nrm",
#' intercept=!is.element(model,c("rrl","ord")),
#' odmax=rep(max(apply(Y>0,c(1,3),sum,na.rm=TRUE)),nrow(Y[,,1])), seed = 1,
#' nscan = 50000, burn = 500, odens = 25, plot=TRUE, print = TRUE, gof=TRUE)
#' @param Y an n x n x T array of relational matrix, where the third dimension correponds to replicates (over time, for example). See
#' model below for various data types.
#' @param Xdyad an n x n x pd x T array of covariates
#' @param Xrow an n x pr x T array of nodal row covariates
#' @param Xcol an n x pc x T array of nodal column covariates
#' @param rvar logical: fit row random effects?
#' @param cvar logical: fit column random effects?
#' @param dcor logical: fit a dyadic correlation?
#' @param R integer: dimension of the multiplicative effects (can be zero)
#' @param model character: one of "nrm","bin","ord","cbin","frn","rrl" - see
#' the details below
#' @param intercept logical: fit model with an intercept?
#' @param odmax a scalar integer or vector of length n giving the maximum
#' number of nominations that each node may make - used for "frn" and "cbin"
#' models
#' @param seed random seed
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param burn burn in for the Markov chain
#' @param odens output density for the Markov chain
#' @param plot logical: plot results while running?
#' @param print logical: print results while running?
#' @param gof logical: calculate goodness of fit statistics?
#' @return \item{BETA}{posterior samples of regression coefficients}
#' \item{SABR}{posterior samples of Cov(a,b) and the dyadic correlation}
#' 
#' \item{APM}{posterior mean of additive row effects a} \item{BPM}{posterior
#' mean of additive column effects b} \item{U}{posterior mean of multiplicative
#' row effects u} \item{V}{posterior mean of multiplicative column effects v}
#' \item{UVPM}{posterior mean of UV} \item{EZ}{estimate of expectation of Z
#' matrix} \item{YPM}{posterior mean of Y (for imputing missing values)}
#' \item{GOF}{observed (first row) and posterior predictive (remaining rows)
#' values of four goodness-of-fit statistics}
#' @author Peter Hoff, Yanjun He
#' @examples
#' 
#' data(YX_bin_long) 
#' fit<-ame_rep(YX_bin_long$Y,YX_bin_long$X,burn=5,nscan=5,odens=1,model="bin")
#' # you should run the Markov chain much longer than this
#' 
#' @export ame_rep 
ame_rep<-
  function(Y,Xdyad=NULL, Xrow=NULL, Xcol=NULL, 
            rvar = !(model=="rrl") , cvar = TRUE, dcor = TRUE, R = 0,
            model="nrm",
            intercept=!is.element(model,c("rrl","ord")),
            odmax=rep(max(apply(Y>0,c(1,3),sum,na.rm=TRUE)),nrow(Y[,,1])),
            seed = 1, nscan = 50000, burn = 500, odens = 25,
            plot=TRUE, print = TRUE, gof=TRUE)
  { 
    
    ## 
    set.seed(seed)
    
    N<-dim(Y)[3]
    
    for (t in 1:N) diag(Y[,,t]) <- NA
    if(is.element(model,c("bin","cbin"))) { Y<-1*(Y>0) } 
    
    if(is.element(model,c("cbin","frn","rrl"))){odobs<-apply(Y>0,c(1,3),sum,na.rm=TRUE)}
    
    n<-nrow(Y[,,1])
    pr<-length(Xrow[,,,1])/n
    pc<-length(Xcol[,,,1])/n
    pd<-length(Xdyad[,,,1])/n^2
    X<-array(dim=c(n,n,pr+pc+pd+intercept,N))
    for (t in 1:N){
      X[,,,t]<-design_array(Xrow[,,,t],Xcol[,,,t],Xdyad[,,,t],intercept,n) 
    }
    
    if(is.element(model,c("rrl","ord"))&any(apply(apply(X,c(1,3,4),var),2,sum)==0))
    {
      cat("WARNING: row effects are not estimable using this procedure ","\n")
    }
    ##
    
    ##
    if(length(odmax)==1) { odmax<-rep(odmax,nrow(Y[,,1])) }
    ## 
    
    ##
    if(is.element(model,c("frn","rrl")))
    {
      ymx<-max(apply(1*(Y>0),c(1,3),sum,na.rm=TRUE))
      YL<-list()
      for (t in 1:N){
        YL.t<-NULL
        warn<-FALSE
        for(i in 1:nrow(Y[,,1]))
        {
          yi<-Y[i,,t] ; rnkd<-which( !is.na(yi)&yi>0 )
          if(length(yi[rnkd])>length(unique(yi[rnkd]))){warn<-TRUE}
          yi[rnkd]<-rank(yi[rnkd],ties.method="random")
          Y[i,,t]<-yi
          YL.t<-rbind(YL.t, match(1:ymx,yi))
        }
        YL[[t]]<-YL.t
        if(warn){cat("WARNING: Random reordering applied to break ties in ranks\n")}
      }
    }
    ##
    
    ## starting Z values
    Z<-array(dim=dim(Y))
    for (t in 1:N){
      if(model=="nrm") { Z[,,t]<-Y[,,t] }
      if(model=="ord") { Z[,,t]<-matrix(zscores(Y[,,t]),nrow(Y[,,t]),ncol(Y[,,t])) } 
      if(model=="rrl") { Z[,,t]<-matrix(t(apply(Y[,,t],1,zscores)),nrow(Y[,,t]),ncol(Y[,,t])) }  
      if(model=="bin")
      { 
        Z[,,t]<-matrix(zscores(Y[,,t]),nrow(Y[,,t]),nrow(Y[,,t])) 
        z01<- .5* ( max(Z[,,t][Y[,,t]==0],na.rm=TRUE) + min(Z[,,t][Y[,,t]==1],na.rm=TRUE) ) 
        Z[,,t]<-Z[,,t] - z01
      } 
      
      if(is.element(model,c("cbin","frn")))
      {
        Z[,,t]<-Y[,,t]
        for(i in 1:nrow(Y[,,t]))
        {
          yi<-Y[i,,t]
          zi<-zscores(yi)
          rnkd<-which( !is.na(yi) & yi>0 ) 
          if(length(rnkd)>0 && min(zi[rnkd])<0) { zi[rnkd]<-zi[rnkd] - min(zi[rnkd]) + 1e-3 }
          
          if(length(rnkd)<odmax[i]) 
          {
            urnkd<-which( !is.na(yi) & yi==0 ) 
            if(max(zi[urnkd])>0) { zi[urnkd]<-zi[urnkd] - max(zi[urnkd]) -1e-3 }
          }
          
          Z[i,,t]<-zi
        } 
      }
    }
    
    
    Z0<-Z #"observed" value of Z
    ZA<-Z
    for (t in 1:N){
      ZA[,,t]<-outer(rowMeans(Z[,,t],na.rm=TRUE),colMeans(Z[,,t],na.rm=TRUE),"+")/2
    }
    Z[is.na(Z)]<-ZA[is.na(Z)] 
    ##
    
    
    ## zeros for other starting values
    beta<-rep(0,dim(X)[3]) 
    a <-b<-rep(0,nrow(Y[,,1])) 
    s2 <-1 
    rho<-0
    Sab<-diag(c(rvar,cvar))
    U<-V<-matrix(0, nrow(Y[,,1]), R)  ##
    
    
    ##  output items
    BETA <- matrix(nrow = 0, ncol = dim(X)[3])
    SABR<-matrix(nrow=0,ncol=5) 
    UVPS <- U %*% t(V) * 0 
    APS<-BPS<- rep(0,nrow(Y[,,1]))  
    YPS<-array(0,dim=dim(Y)) ; dimnames(YPS)<-dimnames(Y)
    GOF<-matrix(gofstats(Y[,,1]),1,4)  
    rownames(GOF)<-"obs"
    colnames(GOF)<- c("sd.rowmean","sd.colmean","dyad.dep","triad.dep")
    
    colnames(BETA) <- dimnames(X)[[3]] 
    colnames(SABR) <- c("va", "cab", "vb", "rho", "ve")
    names(APS)<-names(BPS)<- rownames(U)<-rownames(V)<-rownames(Y[,,1]) 
    
    have_coda<-suppressWarnings(try(requireNamespace("coda",quietly = TRUE),silent=TRUE))
    

    ## MCMC
    for (s in 1:(nscan + burn)) 
    { 
      
      ## update Z
      E.nrm<-array(dim=dim(Z))
      for (t in 1:N){
        EZ<-Xbeta(X[,,,t], beta) + outer(a, b, "+") + U %*% t(V)
        if(model=="nrm"){ Z[,,t]<-rZ_nrm_fc(Z[,,t],EZ,rho,s2,Y[,,t]) ;E.nrm[,,t]<-Z[,,t]-EZ}
        if(model=="bin"){ Z[,,t]<-rZ_bin_fc(Z[,,t],EZ,rho,Y[,,t]) }
        if(model=="ord"){ Z[,,t]<-rZ_ord_fc(Z[,,t],EZ,rho,Y[,,t]) }
        if(model=="cbin"){Z[,,t]<-rZ_cbin_fc(Z[,,t],EZ,rho,Y[,,t],odmax,odobs)}
        if(model=="frn"){ Z[,,t]<-rZ_frn_fc(Z[,,t],EZ,rho,Y[,,t],YL[[t]],odmax,odobs)}
        if(model=="rrl"){ Z[,,t]<-rZ_rrl_fc(Z[,,t],EZ,rho,Y[,,t],YL[[t]])}
      }
      if (model=="nrm") s2<-rs2_rep_fc(E.nrm,rho) 
      
      ## update beta, a b
      UVprod<-U%*%t(V)
      ZmUV<-Z
      for (t in 1:N){
        ZmUV[,,t]<-Z[,,t]-UVprod
      }
      tmp <- rbeta_ab_rep_fc(ZmUV,Sab,rho,X,s2)
      beta <- tmp$beta
      a <- tmp$a * rvar
      b <- tmp$b * cvar 
      
      ## update Sab 
      if(rvar & cvar)
      {
        #if(is.element(model,c("nrm","ord")))
        #{ 
          Sab<-solve(rwish(solve(diag(2)+crossprod(cbind(a, b))),3+nrow(Z[,,1]))) 
        #}
        
        #if(model=="bin")
        #{
        #  tmp<-raSab_bin_fc(Z,Y,a,b,Sab) ; Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
        #}
        
        #if(model=="cbin")
        #{
        #  tmp<-raSab_cbin_fc(Z,Y,a,b,Sab,odmax,odobs)
        #  Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
        #}
        
        #if(model=="frn")
        #{ 
        #  tmp<-raSab_frn_fc(Z,Y,YL,a,b,Sab,odmax,odobs) 
        #  Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a 
        #}
      }
      if (rvar & !cvar) 
      {
        Sab[1, 1] <- 1/rgamma(1, (1 + nrow(Y[,,t]))/2, (1 + sum(a^2))/2)
      }
      
      if (!rvar & cvar) 
      {
        Sab[2, 2] <- 1/rgamma(1, (1 + nrow(Y[,,t]))/2, (1 + sum(b^2))/2)
      }
      
      ## update rho
      if (dcor) 
      {
        E.T<-array(dim=dim(Z))
        for (t in 1:N){
          E.T[,,t]<-Z[,,t]-(Xbeta(X[,,,t], beta) + outer(a, b, "+") + U %*% t(V))
        }
        rho<-rrho_mh_rep(E.T, rho,s2)
      }
      
      ## update U & V
      if (R > 0) 
      {
        E.T2<-array(dim=dim(Z))
        for (t in 1:N){
          E.T2[,,t]<-Z[,,t] - (Xbeta(X[,,,t], beta) + outer(a, b, "+"))
        }
        UV <- rUV_rep_fc(E.T2, U, V, rho, s2)
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
        EZ<-Ys<-array(dim=dim(Z))
        for (t in 1:N){
          EZ[,,t]<-Xbeta(X[,,,t], beta) + outer(a, b, "+") + U %*% t(V)
          if(model=="bin") { Ys[,,t]<-simY_bin(EZ[,,t],rho) }
          if(model=="cbin"){ Ys[,,t]<-1*(simY_frn(EZ[,,t],rho,odmax,YO=Y[,,t])>0) }
          if(model=="frn") { Ys[,,t]<-simY_frn(EZ[,,t],rho,odmax,YO=Y[,,t]) }
          if(model=="rrl") { Ys[,,t]<-simY_rrl(EZ[,,t],rho,odobs,YO=Y[,,t] ) }
          if(model=="nrm") { Ys[,,t]<-simY_nrm(EZ[,,t],rho,s2) }
          if(model=="ord") { Ys[,,t]<-simY_ord(EZ[,,t],rho,Y[,,t]) }
        }
        
        YPS<-YPS+Ys
        
        if(gof)
        { 
          Ys[is.na(Y)]<-NA
          gof.s<-apply(Ys,3,gofstats)
          GOF<-rbind(GOF,rowMeans(gof.s)) }
        
        if (print) 
        {
          cat(s,round(apply(BETA,2,mean),2),":",round(apply(SABR,2,mean),2),"\n")
          if (have_coda & nrow(SABR) > 3 & length(beta)>0) 
          {
            cat(round(coda::effectiveSize(BETA)), "\n")
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
    names(APM)<-names(BPM)<-rownames(U)<-rownames(V)<-rownames(Y[,,1])
    EZ<-YPM<-array(dim=dim(Y))
    for (t in 1:N){
      EZ[,,t]<-Xbeta(X[,,,t],apply(BETA,2,mean)) + outer(APM,BPM,"+")+UVPM
      dimnames(EZ[,,t])<-dimnames(Y[,,t]) 
      YPM[,,t]<-YPS[,,t]/nrow(SABR)
    }
    
    
    fit <- list(BETA=BETA,SABR=SABR,APM=APM,BPM=BPM,U=U,V=V,UVPM=UVPM,EZ=EZ,
                YPM=YPM,GOF=GOF)
    #fit <- list(BETA=BETA,SABR=SABR,APM=APM,BPM=BPM,U=U,V=V,UVPM=UVPM,EZ=EZ,
    #            YPM=YPM)
    class(fit) <- "ame"
    fit
  }



