simY_ord<-
function (EZ, rho,Y) 
{ 
    
    uY<-sort(unique(c(Y))) 
    FY<-table(c(Y)) ; FY<-FY/sum(FY) ; FY<-cumsum(FY) 
    ZS <- simZ(EZ, rho)
    diag(ZS) <- NA
    qZ <- quantile(ZS, FY[-length(FY)], na.rm = TRUE)
    YS <- ZS * 0 + max(uY)
    for (w in 1:(length(uY) - 1)) {
        YS[(YS == max(uY)) & ZS <= qZ[w]] <- uY[w]
    }
    YS
}



