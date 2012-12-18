\name{rmvnorm}
\alias{rmvnorm}
\title{Simulation from a multivariate normal distribution}
\description{Simulates a matrix where the rows are i.i.d. samples from a multivariate normal distribution }
\usage{
rmvnorm(n, mu, Sigma, Sigma.chol = chol(Sigma))
}
\arguments{
  \item{n}{sample size }

  \item{mu}{multivariate mean  vector
}
  \item{Sigma}{
covariance matrix 
}
  \item{Sigma.chol}{
Cholesky factorization of \code{Sigma}
}
}
\value{
a matrix with \code{n} rows
}
\author{
Peter Hoff
}
\examples{
## The function is currently defined as
function (n, mu, Sigma, Sigma.chol = chol(Sigma)) 
{
    E <- matrix(rnorm(n * length(mu)), n, length(mu))
    t(t(E \%*\% Sigma.chol) + c(mu))
  }
}
