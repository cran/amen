\name{rZ_frn_fc}
\alias{rZ_frn_fc}
\title{
Simulate Z given fixed rank nomination data
}
\description{
Simulates a random latent matrix Z given its expectation, dyadic correlation and 
fixed rank nomination data
}
\usage{
rZ_frn_fc(Z, EZ, rho, Y, YL, odmax, odobs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
a square matrix, the current value of Z
}
  \item{EZ}{
expected value of Z
}
  \item{rho}{
dyadic correlation
}
  \item{Y}{
square matrix of ranked nomination data
}
  \item{YL}{
list of ranked individuals, from least to most preferred in each row
}
  \item{odmax}{
a scalar or vector giving the maximum number of nominations for each individual
}
  \item{odobs}{
observed outdegree
}
}
\details{
  simulates Z under the constraints
   (1)  Y[i,j]>Y[i,k]                => Z[i,j]>Z[i,k]  ,  
   (2)  Y[i,j]>0                     => Z[i,j]>0     , 
   (3)  Y[i,j]=0 & odobs[i]<odmax[i] => Z[i,j]<0
}
\value{
a square matrix, the new value of Z
}
\author{Peter Hoff}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## The function is currently defined as
function (Z, EZ, rho, Y, YL, odmax, odobs) 
{
    sz <- sqrt(1 - rho^2)
    ut <- upper.tri(Z)
    lt <- lower.tri(Z)
    rws <- outer(1:nrow(Z), rep(1, nrow(Z)))
    for (y in sample(c(0:ncol(YL)))) {
        if (y < 2) {
            if (y == 0) {
                lbm <- rep(-Inf, nrow(Z))
            }
            if (y == 1) {
                lbm <- pmax(0, apply(Z - (Y != 0) * (Inf^(Y != 
                  0)), 1, max, na.rm = TRUE))
            }
        }
        if (y >= 2) {
            lbm <- Z[cbind(1:nrow(Z), YL[, y - 1])]
        }
        if (y < ncol(YL)) {
            ubm <- Z[cbind(1:nrow(Z), YL[, y + 1])]
        }
        if (y == 0) {
            ubm[odobs < odmax] <- 0
        }
        if (y == ncol(YL)) {
            ubm <- rep(Inf, nrow(Z))
        }
        ubm[is.na(ubm)] <- Inf
        lbm[is.na(lbm)] <- -Inf
        for (k in sample(1:2)) {
            if (k == 1) {
                up <- ut & Y == y
                rwb <- rws[up]
                lb <- lbm[rwb]
                ub <- ubm[rwb]
                ez <- EZ[up] + rho * (t(Z)[up] - t(EZ)[up])
                Z[up] <- ez + sz * qnorm(runif(sum(up), pnorm((lb - 
                  ez)/sz), pnorm((ub - ez)/sz)))
            }
            if (k == 2) {
                up <- lt & Y == y
                rwb <- rws[up]
                lb <- lbm[rwb]
                ub <- ubm[rwb]
                ez <- EZ[up] + rho * (t(Z)[up] - t(EZ)[up])
                Z[up] <- ez + sz * qnorm(runif(sum(up), pnorm((lb - 
                  ez)/sz), pnorm((ub - ez)/sz)))
            }
        }
    }
    diag(Z) <- rnorm(nrow(Z), diag(EZ), 1)
    Z
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
