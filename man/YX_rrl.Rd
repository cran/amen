\name{YX_rrl}
\alias{YX_rrl}
\docType{data}
\title{
row-specific 
ordinal relational data and covariates
}
\description{
a synthetic dataset that includes row-specific ordinal relational data as well as information
 on five covariates
}
\usage{data(YX_rrl)}
\format{
  The format is:
List of 2
 $ Y: num [1:100, 1:100] NA 0 3 0 3 1 0 1 1 0 ...
 $ X: num [1:100, 1:100, 1:5] 1 1 1 1 1 1 1 1 1 1 ...
  ..- attr(*, "dimnames")=List of 3
  .. ..$ : NULL
  .. ..$ : NULL
  .. ..$ : chr [1:5] "cgpa" "csmoke" "igrade" "ismoke" ...
}
\examples{
data(YX_rrl)
gofstats(YX_rrl$Y)
}
\keyword{datasets}
