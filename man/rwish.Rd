\name{rwish}
\alias{rwish}
\title{Simulation from a Wishart distribution}
\description{Simulates a random Wishart-distributed matrix}
\usage{rwish(S0, nu = dim(S0)[1] + 1)}
\arguments{
  \item{S0}{a positive definite matrix}
  \item{nu}{a positive integer}
}

\value{a positive definite matrix}

\author{Peter Hoff}

\examples{
## The expectation is S0*nu

S0<-rwish(diag(3)) 

SS<-matrix(0,3,3) 
for(s in 1:1000) { SS<-SS+rwish(S0,5) }

SS/s 

S0*5

}
