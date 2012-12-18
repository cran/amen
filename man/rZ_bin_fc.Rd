\name{rZ_bin_fc}
\alias{rZ_bin_fc}
\title{
Simulate Z based on a probit model
}
\description{
Simulates a random latent matrix Z given its expectation, dyadic correlation and a binary relational matrix Y
}
\usage{
rZ_bin_fc(Z, EZ, rho, Y)
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
square binary relational matrix 
}
}
\value{
a square matrix , the new value of Z
}
\author{
Peter Hoff
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Z, EZ, rho, Y) 
{
    sz <- sqrt(1 - rho^2)
    ut <- upper.tri(EZ)
    lt <- lower.tri(EZ)
    for (y in 0:1) {
        lb <- c(-Inf, 0)[y + 1]
        ub <- c(0, Inf)[y + 1]
        up <- ut & Y == y
        ez <- EZ[up] + rho * (t(Z)[up] - t(EZ)[up])
        Z[up] <- ez + sz * qnorm(runif(sum(up), pnorm((lb - ez)/sz), 
            pnorm((ub - ez)/sz)))
        up <- lt & Y == y
        ez <- EZ[up] + rho * (t(Z)[up] - t(EZ)[up])
        Z[up] <- ez + sz * qnorm(runif(sum(up), pnorm((lb - ez)/sz), 
            pnorm((ub - ez)/sz)))
    }
    diag(Z) <- rnorm(nrow(Z), diag(EZ), 1)
    Z
  }
}
