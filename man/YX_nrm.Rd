\name{YX_nrm}
\alias{YX_nrm}
\docType{data}
\title{
normal relational data and covariates
}
\description{
a synthetic dataset that includes continuous (normal) relational data as well as information
 on eight covariates
}
\usage{data(YX_nrm)}
\format{
  The format is:
List of 2
 $ Y: num [1:100, 1:100] NA -4.05 -0.181 -3.053 -1.579 ...
 $ X: num [1:100, 1:100, 1:8] 1 1 1 1 1 1 1 1 1 1 ...
  ..- attr(*, "dimnames")=List of 3
  .. ..$ : NULL
  .. ..$ : NULL
  .. ..$ : chr [1:8] "intercept" "rgpa" "rsmoke" "cgpa" ...
}
\examples{
data(YX_nrm)
t_trans(YX_nrm$Y>quantile(YX_nrm$Y,.95,na.rm=TRUE))
}
\keyword{datasets}
