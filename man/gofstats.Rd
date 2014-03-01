\name{gofstats}
\alias{gofstats}
\title{
Goodness of fit statistics 
}
\description{
Goodness of fit statistics  evaluating second and third-order dependence patterns
}
\usage{
gofstats(Y) 
}
\arguments{
  \item{Y}{ 
a relational data matrix
}
}
\value{
a vector of gof statistics
}
\author{
Peter Hoff
}
\examples{
data(YX_nrm) 

gofstats(YX_nrm$Y) 

}
