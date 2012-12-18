\name{t_recip}
\alias{t_recip}
\title{Reciprocity statistic}
\description{For square binary matrices, computes the fraction of reciprocated dyads}
\usage{ t_recip(Y) }
\arguments{
  \item{Y}{A square binary matrix with a missing diagonal}
}
\value{The fraction of reciprocated dyads in the network}
\author{Peter Hoff}
\examples{ 
data(YX_bin) 
t_recip(YX_bin$Y)  }
