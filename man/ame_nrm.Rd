\name{ame_nrm}
\alias{ame_nrm}
\title{
AME fit for normal relational data
}
\description{
An MCMC routine providing a fit to an additive and multiplicative effects
(AME) regression model for normal relational data
}
\usage{
ame_nrm(Y, X, rvar = TRUE, cvar = TRUE, dcor = TRUE, R = 0, seed = 1, nscan = 50000, burn = 500, odens = 25, plot = TRUE, print = TRUE)
}
\arguments{
  \item{Y}{
an n x n square relational matrix 
}
  \item{X}{
an n x n x p array of covariates
}
  \item{rvar}{
logical: fit row random effects?
}
  \item{cvar}{
logical: fit column random effects?
}
  \item{dcor}{
logical: fit a dyadic correlation?
}
  \item{R}{
integer: dimension of the multiplicative effects (can be zero)
}
  \item{seed}{
random seed
}
  \item{nscan}{
number of iterations of the Markov chain (beyond burn-in)
}
  \item{burn}{
burn in for the Markov chain
}
  \item{odens}{
output density for the Markov chain 
}
  \item{plot}{
logical: plot results while running?
}
  \item{print}{
logical: print results while running?
}
}
\value{
\item{BETA}{posterior samples of regression coefficients}
\item{SABR}{posterior samples of Cov(a,b) and the dyadic correlation}
\item{UVPM}{posterior mean of UV'}
\item{APM}{posterior mean of additive row effects a}
\item{BPM}{posterior mean of additive column effects b}
\item{TT}{posterior predictive of number of triangles}
\item{TR}{posterior predictive of fraction of reciprocated dyads}
\item{TID}{posterior predictive of indegree distribution}
\item{TOD}{posterior predictive of outdegree distribution}
\item{tt}{observed triangles}
\item{tr}{observed reciprocity fraction}
\item{td}{observed out and indegree distributions}
}
\author{
Peter Hoff
}

\examples{
data(YX_nrm)
fit<-ame_nrm(YX_nrm$Y,YX_nrm$X,burn=5,nscan=5,odens=1)  
# you should run the Markov chain much longer than this

 }
