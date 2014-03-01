rZ_ord_fc<-
function (Z, EZ, rho,Y)
{ 
    # this could be done outside this function
    uY<-sort(unique(c(Y)))
    W<-matrix(match(Y,uY),nrow(Y),nrow(Y))
    W[is.na(W)] <- -1

    sz <- sqrt(1 - rho^2)
    ut <- upper.tri(Z)
    lt <- lower.tri(Z)
    for (w in sample(c(-1, 1:length(uY)))) {
        lb <- suppressWarnings(max(Z[!is.na(W) & W == w - 1],
            na.rm = TRUE))
        ub <- suppressWarnings(min(Z[!is.na(W) & W == w + 1],
            na.rm = TRUE))
        up <- ut & W == w
        ez <- EZ[up] + rho * (t(Z)[up] - t(EZ)[up])
        Z[up] <- ez + sz * qnorm(runif(sum(up), pnorm((lb - ez)/sz),
            pnorm((ub - ez)/sz)))
        up <- lt & W == w
        ez <- EZ[up] + rho * (t(Z)[up] - t(EZ)[up])
        Z[up] <- ez + sz * qnorm(runif(sum(up), pnorm((lb - ez)/sz),
            pnorm((ub - ez)/sz)))
    }
    diag(Z) <- rnorm(nrow(Z), diag(EZ), 1)
    Z
}

