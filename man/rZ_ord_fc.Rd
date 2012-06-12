\name{rZ_ord_fc}
\alias{rZ_ord_fc}
\title{
Simulate Z given the partial ranks
}
\description{
Simulates a random latent matrix Z given its expectation, dyadic correlation and partial rank information provided by W 
}
\usage{
rZ_ord_fc(Z, EZ, rho, W, FY)
}
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
  \item{W}{
matrix of ordinal data
}
  \item{FY}{
cumulative CDF of distribution of W (note this isnt really needed) 
}
}
\value{
a square matrix, the new value of Z
}
\author{
Peter Hoff
}
\examples{
## The function is currently defined as
function (Z, EZ, rho, W, FY) 
{
    sz <- sqrt(1 - rho^2)
    ut <- upper.tri(Z)
    lt <- lower.tri(Z)
    for (w in sample(1:length(FY))) {
        lb <- suppressWarnings(max(Z[W == w - 1], na.rm = TRUE))
        ub <- suppressWarnings(min(Z[W == w + 1], na.rm = TRUE))
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
}
