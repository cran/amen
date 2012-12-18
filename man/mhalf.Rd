\name{mhalf}
\alias{mhalf}
\title{
Symmetric square root of a matrix
}
\description{
Computes the symmetric square root of a positive definite matrix 
}
\usage{
mhalf(M)
}
\arguments{
  \item{M}{
a positive definite matrix
}
}
\value{
a matrix \code{H} such that \code{H^2} equals \code{M}
}
\author{
Peter Hoff
}
\examples{
## The function is currently defined as
function (M) 
{
    tmp <- eigen(M)
    tmp$vec \%*\% sqrt(diag(tmp$val, nrow = nrow(M))) \%*\% t(tmp$vec)
  }
}
