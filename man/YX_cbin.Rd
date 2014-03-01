\name{YX_cbin}
\alias{YX_cbin}
\docType{data}
\title{
Censored binary nomination data and covariates
}
\description{
a synthetic dataset that includes 
relational data where the number of nominations per row is censored at 10, 
along with information on eight covariates
}
\usage{data(YX_cbin)}
\format{
  The format is:
List of 2
 $ Y: num [1:100, 1:100] NA 0 0 0 1 0 0 0 0 3 ...
 $ X: num [1:100, 1:100, 1:8] 1 1 1 1 1 1 1 1 1 1 ...
  ..- attr(*, "dimnames")=List of 3
  .. ..$ : NULL
  .. ..$ : NULL
  .. ..$ : chr [1:8] "intercept" "rgpa" "rsmoke" "cgpa" ...
}
\examples{
data(YX_cbin)
gofstats(YX_cbin$Y) 
}
\keyword{datasets}
