\name{t_degree}
\alias{t_degree}
\title{Degree distribution}
\description{For square binary matrices, computes the outdegree and indegree distributions}
\usage{ t_degree(Y) }
\arguments{
  \item{Y}{A square binary matrix with a missing diagonal}
}
\value{A list with the following components:
\item{od}{outdegree distribution}
\item{id}{indegree distribution}
}
\author{Peter Hoff}
\examples{ 
data(YX_bin)
t_degree(YX_bin$Y)  }
