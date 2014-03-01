\name{YX_bin}
\alias{YX_bin}
\docType{data}
\title{
binary relational data and covariates
}
\description{
a synthetic dataset that includes binary relational data as well as information
 on eight covariates
}
\usage{data(YX_bin)}
\format{
  The format is:
List of 2
 $ Y: num [1:100, 1:100] NA 0 0 0 0 0 0 0 0 1 ...
 $ X: num [1:100, 1:100, 1:8] 1 1 1 1 1 1 1 1 1 1 ...
  ..- attr(*, "dimnames")=List of 3
  .. ..$ : NULL
  .. ..$ : NULL
  .. ..$ : chr [1:8] "intercept" "rgpa" "rsmoke" "cgpa" ...
}
\examples{
data(YX_bin)
gofstats(YX_bin$Y) 
}
\keyword{datasets}
