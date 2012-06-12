\name{t_trans}
\alias{t_trans}
\title{Transitivity statistic}
\description{For square binary matrices, computes the number of triangles}
\usage{ t_trans(Y) }
\arguments{
  \item{Y}{A square binary matrix with a missing diagonal}
}
\value{An integer number of triangles in the network}
\author{Peter Hoff}
\examples{
data(YX_bin) 
t_trans(YX_bin$Y)  }
