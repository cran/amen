% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rZ_cbin_fc.R
\name{rZ_cbin_fc}
\alias{rZ_cbin_fc}
\title{Simulate Z given fixed rank nomination data}
\usage{
rZ_cbin_fc(Z, EZ, rho, Y, odmax, odobs)
}
\arguments{
\item{Z}{a square matrix, the current value of Z}

\item{EZ}{expected value of Z}

\item{rho}{dyadic correlation}

\item{Y}{square matrix of ranked nomination data}

\item{odmax}{a scalar or vector giving the maximum number of nominations for
each individual}

\item{odobs}{observed outdegree}
}
\value{
a square matrix, the new value of Z
}
\description{
Simulates a random latent matrix Z given its expectation, dyadic correlation
and censored binary nomination data
}
\details{
simulates Z under the constraints (1) Y[i,j]=1, Y[i,k]=0 => Z[i,j]>Z[i,k] ,
(2) Y[i,j]=1 => Z[i,j]>0 , (3) Y[i,j]=0 & odobs[i]<odmax[i] => Z[i,j]<0
}
\author{
Peter Hoff
}
