ame_nrm<-
function (Y, X=NULL, Xrow=NULL, Xcol=NULL, 
    rvar = TRUE, cvar = TRUE, dcor = TRUE, R = 0, 
    seed = 1, nscan = 50000, burn = 500, odens = 25, plot = TRUE, 
    print = TRUE,intercept=TRUE) 
{ 
set.seed(seed)


### combine dyad, row and col covariates
pd<-length(X)/nrow(Y)^2
pr<-length(Xrow)/nrow(Y)
pc<-length(Xcol)/nrow(Y)

if(pd+pr+pc>0)
{
  Xall<-array(dim=c(nrow(Y),nrow(Y),pr+pc+pd))
  dnX<-NULL

  if(pd>0)
  {
    Xall[,,1:pd]<-X
    dnX<-paste0(dimnames(X)[[3]],rep(".dyad",pd))
  }

  if(pr>0)
  {
    if(length(Xrow)==nrow(Y)) { Xrow<-matrix(Xrow,nrow(Y),1)  }
    Xrowa<-array(dim=c(nrow(Y),nrow(Y),ncol(Xrow)))
    for(j in 1:ncol(Xrow)){ Xrowa[,,j]<-matrix( Xrow[,j], nrow(Y),nrow(Y)) }
    Xall[ ,,pd+1:pr]<- Xrowa
    dnX<-c(dnX,paste0(colnames(Xrow),rep(".row" ,pr)))
  }

  if(pc>0)
  {
    if(length(Xcol)==nrow(Y)) { Xcol<-matrix(Xcol,nrow(Y),1)  }
    Xcola<-array(dim=c(nrow(Y),nrow(Y),ncol(Xcol)))
    for(j in 1:ncol(Xcol)){ Xcola[,,j]<-t(matrix( Xcol[,j], nrow(Y),nrow(Y))) }
    Xall[ ,,pd+pr+1:pc]<- Xcola
    dnX<-c(dnX,paste0(colnames(Xcol),rep(".col" ,pc)))
  }

  if(pd+pr+pc>1) { dimnames(Xall)[[3]]<- dnX }
  if(pd+pr+pc==1){ dimnames(Xall)[[3]]<- list(dnX) }
  X<-Xall
}
###





    ## create intercept if X is NULL and intercept is TRUE
    if(is.null(X) & intercept) {  X<-matrix(1,nrow(Y),nrow(Y)) }
    if(is.null(X) & !intercept) {  X<-array(dim=c(nrow(Y),nrow(Y),0)) }
 
    diag(Y) <- NA

    if(length(dim(X))==2 ) { X<-array(X,dim=c(dim(X),1)) }
    if(dim(X)[1]>0 &  !any(apply(X,3,function(x){var(c(x))})==0) )
    {
      if(!intercept){ cat("WARNING: design matrix lacks an intercept","\n")  }
      if(intercept)  # add an intercept
      {
      X1<-array(dim=c(0,0,1)+dim(X))
      X1[,,1]<-1 ; X1[,,-1]<-X
      X<-X1
      }
    }

    Xr <- apply(X, c(1, 3), sum)
    Xc <- apply(X, c(2, 3), sum)
    mX <- apply(X, 3, c)
    mXt <- apply(aperm(X, c(2, 1, 3)), 3, c)
    XX <- t(mX) %*% mX
    XXt <- t(mX) %*% mXt 
 
    if(dim(X)[3]>0)
    { 
    fit <- lm(c(Y) ~ -1 + apply(X, 3, c))   
    beta<-fit$coef
    res<-fit$res 
    }

    if(dim(X)[3]==0)
    { 
    beta<-numeric(0) 
    res<- Y[!is.na(Y)]  
    }

    E <- matrix(NA, nrow(Y), ncol(Y))
    E[!is.na(Y)] <- res
    a <- apply(E, 1, mean, na.rm = TRUE)*rvar
    b <- apply(E, 2, mean, na.rm = TRUE)*cvar
    E <- E - outer(a, b, "+")
    CE <- cov(cbind(E[upper.tri(E)], t(E)[upper.tri(E)]),use="complete.obs")
    s2 <- 0.5 * (CE[1, 1] + CE[2, 2])
    rho <- dcor*(.99)*cov2cor(CE)[1, 2]
    Sab<-cov(cbind(a,b)) + diag(2)*rvar*cvar 
    Z<-Y 
    EZ <- Xbeta(X, beta) + outer(a, b, "+")
    ZS <- simY_nrm(EZ, rho, s2)
    diag(ZS) <- rnorm(nrow(Y), diag(EZ), sqrt(s2)) 
    Z[is.na(Y)]<-ZS[is.na(Y)]

    U <- V <- matrix(0, nrow(Y), R) 

    YT<- 1*(Y>quantile(c(Y),.95,na.rm=TRUE))
    tr.obs <- t_recip(YT)
    tt.obs <- t_trans(YT)
    td.obs <- t_degree(YT)
    odobs <- apply(YT, 1, sum, na.rm = TRUE)
    idobs <- apply(YT, 2, sum, na.rm = TRUE)
    TT <- TR <- TID <- TOD <- SABR <- NULL
    BETA <- matrix(nrow = 0, ncol = dim(X)[3])
    colnames(BETA) <- dimnames(X)[[3]]
    UVPS <- U %*% t(V) * 0
    APS <- BPS <- rep(0, nrow(Y))
    have.coda <- suppressWarnings(try(require(coda, quietly = TRUE), 
        silent = TRUE))
for (s in 1:(nscan + burn)) 
{ 
  ## update beta, a b
  tmp <- rbeta_ab_fc(Z-U%*%t(V), Sab, rho, X, mX, mXt,XX, XXt, Xr, Xc, s2)
  beta <- tmp$beta
  a <- tmp$a * rvar
  b <- tmp$b * cvar 

  ## update Z - impute diagonal and other missing values 
  EZ <- Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V) 
  ZS <- simY_nrm(EZ, rho, s2)  
  diag(ZS) <- rnorm(nrow(Y), diag(EZ), sqrt(s2))
  Z[is.na(Y)] <- ZS[is.na(Y)]  

  ## update Sab 
  if(rvar & cvar)
  {
    Sab<-solve(rwish(solve(diag(2) + crossprod(cbind(a, b))), 3 + nrow(Z)))
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
    rho <- rrho_mh(Z - (Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V)), rho,s2)
  }
  s2 <- rs2_fc(Z - (Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V)), rho)

  ## update U,V
  if (R > 0) 
  {
    UV <- rUV_fc(Z - (Xbeta(X, beta) + outer(a, b, "+")), U, V, rho, s2)
    U <- UV$U
    V <- UV$V
  } 

  ## output
        if (s%%odens == 0 & s <= burn) {
            cat(round(100 * s/burn, 2), " pct burnin complete \n")
        }
        if (s%%odens == 0 & s > burn) { 

            BETA <- rbind(BETA, beta)
            SABR <- rbind(SABR, c(Sab[upper.tri(Sab, diag = T)], 
                rho, s2))
            UVPS <- UVPS + U %*% t(V)
            APS <- APS + a
            BPS <- BPS + b
            YS <- simY_nrm(Xbeta(X, beta) + outer(a, b, "+") + 
                U %*% t(V), rho, s2)
            YS[is.na(Y)] <- NA  
            YTS<- 1*(YS>quantile(c(YS),.95,na.rm=TRUE))

            TR <- c(TR, t_recip(YTS))
            TT <- c(TT, t_trans(YTS))
            td <- t_degree(YTS)
            TOD <- rbind(TOD, td$od)
            TID <- rbind(TID, td$id)
            if (print) {
                cat(s, round(apply(BETA, 2, mean), 2), ":", round(apply(SABR, 
                  2, mean), 2), "\n")
                if (have.coda & length(TR) > 3  & length(beta)>0) {
                  cat(round(effectiveSize(BETA)), "\n")
                }
            }
            if (plot) {
                par(mfrow = c(4, 2), mar = c(3, 3, 1, 1), mgp = c(1.75, 
                  0.75, 0)) 
                hist(a, main = "", col = "lightblue", prob = TRUE)
                hist(b, main = "", col = "lightblue", prob = TRUE)
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

                mod <- max(odobs)
                mid <- max(idobs)
                qod <- apply(TOD, 2, quantile, prob = c(0.975, 
                  0.5, 0.25))
                qid <- apply(TID, 2, quantile, prob = c(0.975, 
                  0.5, 0.25))
                odd <- rbind(qod, td.obs$od)[, 0:(mod + 1)]
                plot(c(0, mod), range(odd), type = "n")
                for (k in 1:4) {
                  lines(0:mod, odd[k, ], col = c("gray", "black")[1 + 
                    (k == 4)])
                }
                idd <- rbind(qid, td.obs$id)[, 0:(mid + 1)]
                plot(c(0, mid), range(idd), type = "n")
                for (k in 1:4) {
                  lines(0:mid, idd[k, ], col = c("gray", "black")[1 + 
                    (k == 4)])
                }
                plot(TR, ylim = range(c(TR, tr.obs)), type = "l")
                abline(h = tr.obs)
                plot(TT, ylim = range(c(TT, tt.obs)), type = "l")
                abline(h = tt.obs)
            }
        }
    }
    colnames(SABR) <- c("va", "cab", "vb", "rho", "ve")

###
UVPM<-UVPS/length(TT)
UDV<-svd(UVPM)
U<-UDV$u[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
V<-UDV$v[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
A<-APS/length(TT)
B<-BPS/length(TT) 
rownames(A)<-rownames(B)<-rownames(U)<-rownames(V)<-rownames(Y)
EZ<-Xbeta(X,apply(BETA,2,mean)) + outer(A,B,"+")+U%*%t(V)
dimnames(EZ)<-dimnames(Y) 
###


    fit <- list(BETA = BETA, SABR = SABR, A=A,B=B,U=U,V=V,UVPM=UVPM,EZ=EZ,
        TT = TT, 
        TR = TR, TID = TID, TOD = TOD, tt = tt.obs, tr = tr.obs, 
        td = td.obs)
    class(fit) <- "ame"
    fit
}



