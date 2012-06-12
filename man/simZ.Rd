\name{simZ}
\alias{simZ}
\title{
Simulate Z given its expectation and covariance
}
\description{
Simulate Z given its expectation and covariance
}
\usage{
simZ(EZ, rho, s2 = 1)
}
\arguments{
  \item{EZ}{
expected value of Z
}
  \item{rho}{
dyadic correlation
}
  \item{s2}{
dyadic variance
}
}
\value{
a simulated value of Z
}
\author{
Peter Hoff
}
\examples{
## The function is currently defined as
function (EZ, rho, s2 = 1) 
{
    w1 <- sqrt((1 + sqrt(1 - rho^2))/2)
    w2 <- sign(rho) * sqrt(1 - w1^2)
    EC <- matrix(rnorm(length(EZ)), nrow(EZ), nrow(EZ))
    EC <- sqrt(s2) * (w1 * EC + w2 * t(EC))
    EZ + EC
  }
}
