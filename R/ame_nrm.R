ame_nrm<-
function (Y, X, rvar = TRUE, cvar = TRUE, dcor = TRUE, R = 0, 
    seed = 1, nscan = 50000, burn = 500, odens = 25, plot = TRUE, 
    print = TRUE) 
{
    diag(Y) <- NA
    if (length(dim(X)) == 2) {
        X <- array(X, dim = c(dim(X), 1))
    }
    if (!any(apply(X, 3, function(x) {
        var(c(x))
    }) == 0)) {
        cat("WARNING: design matrix lacks an intercept", "\n")
    }
    Xr <- apply(X, c(1, 3), sum)
    Xc <- apply(X, c(2, 3), sum)
    mX <- apply(X, 3, c)
    mXt <- apply(aperm(X, c(2, 1, 3)), 3, c)
    XX <- t(mX) %*% mX
    XXt <- t(mX) %*% mXt
    fit <- lm(c(Y) ~ -1 + apply(X, 3, c))
    E <- matrix(0, nrow(Y), ncol(Y))
    E[!is.na(Y)] <- fit$res
    a <- apply(E, 1, mean, na.rm = TRUE)
    b <- apply(E, 2, mean, na.rm = TRUE)
    E <- E - outer(a, b, "+")
    CE <- cov(cbind(E[upper.tri(E)], t(E)[upper.tri(E)]))
    s2 <- 0.5 * (CE[1, 1] + CE[2, 2])
    rho <- cov2cor(CE)[1, 2]
    beta <- fit$coef
    Sab <- cov(cbind(a, b))
    EZ <- Xbeta(X, beta) + outer(a, b, "+")
    Z <- simY_nrm(EZ, rho, s2)
    diag(Z) <- rnorm(nrow(Y), diag(EZ), sqrt(s2))
    U <- V <- matrix(0, nrow(Y), R) 

    YT<- 1*(Y>quantile(c(Y),.95,na.rm=TRUE))
    tr.obs <- t_recip(YT)
    tt.obs <- t_trans(YT)
    td.obs <- t_degree(YT)
    odobs <- apply(YT, 1, sum, na.rm = TRUE)
    idobs <- apply(YT, 2, sum, na.rm = TRUE)
    set.seed(seed)
    TT <- TR <- TID <- TOD <- SABR <- NULL
    BETA <- matrix(nrow = 0, ncol = dim(X)[3])
    colnames(BETA) <- dimnames(X)[[3]]
    UVPS <- U %*% t(V) * 0
    APS <- BPS <- rep(0, nrow(Y))
    have.coda <- suppressWarnings(try(require(coda, quietly = TRUE), 
        silent = TRUE))
    for (s in 1:(nscan + burn)) {
        tmp <- rbeta_ab_fc(Z - U %*% t(V), Sab, rho, X, mX, mXt, 
            XX, XXt, Xr, Xc, s2)
        beta <- tmp$beta
        a <- tmp$a * rvar
        b <- tmp$b * cvar
        EZ <- Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V)
        ZS <- simY_nrm(EZ, rho, s2)
        Z[is.na(Y)] <- ZS[is.na(Y)]
        diag(Z) <- rnorm(nrow(Y), diag(EZ), sqrt(s2))
        if (rvar & cvar) {
            Sab <- solve(rwish(solve(diag(2) + crossprod(cbind(a, 
                b))), 3 + nrow(Z)))
        }
        if (rvar & !cvar) {
            Sab[1, 1] <- 1/rgamma(1, (1 + nrow(Y))/2, (1 + sum(a^2))/2)
        }
        if (!rvar & cvar) {
            Sab[2, 2] <- 1/rgamma(1, (1 + nrow(Y))/2, (1 + sum(b^2))/2)
        }
        if (dcor) {
            rho <- rrho_mh(Z - (Xbeta(X, beta) + outer(a, b, 
                "+") + U %*% t(V)), rho,s2)
        }
        s2 <- rs2_fc(Z - (Xbeta(X, beta) + outer(a, b, "+") + 
            U %*% t(V)), rho)
        if (R > 0) {
            UV <- rUV_fc(Z - (Xbeta(X, beta) + outer(a, b, "+")), 
                U, V, rho, s2)
            U <- UV$U
            V <- UV$V
        }
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
                if (have.coda & length(TR) > 3) {
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
                mBETA <- apply(BETA, 2, median)
                matplot(BETA, type = "l", lty = 1, col = 1:length(mBETA))
                abline(h = mBETA, col = 1:length(mBETA))
                abline(h = 0, col = "gray")
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
    fit <- list(BETA = BETA, SABR = SABR, UVPM = UVPS/length(TT), 
        APM = APS/length(TT), BPM = BPS/length(TT), TT = TT, 
        TR = TR, TID = TID, TOD = TOD, tt = tt.obs, tr = tr.obs, 
        td = td.obs)
    class(fit) <- "ame"
    fit
}



