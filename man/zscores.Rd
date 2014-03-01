\name{zscores}
\alias{zscores}
\title{
rank-based z-scores
}
\description{
Computes the normal scores corresponding to the ranks of a data vector
}
\usage{
zscores(y)
}
\arguments{
  \item{y}{
a numeric vector
}
}
\value{ 
a numeric vector 
}
\author{
Peter Hoff
}
\examples{
## The function is currently defined as 
zscores<-function(y)
{
 qnorm( rank(y,na.last="keep",ties.method="average")/(1+sum(!is.na(y))) )
}

}
