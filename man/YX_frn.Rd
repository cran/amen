\name{YX_frn}
\alias{YX_frn}
\docType{data}
\title{
Fixed rank nomination data and covariates
}
\description{
a synthetic dataset that includes fixed rank nomination data as well as information on eight covariates
}
\usage{data(YX_frn)}
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
data(YX_frn)
t_trans(YX_frn$Y>0 ) 
}
\keyword{datasets}
