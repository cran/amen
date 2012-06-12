\name{YX_ord}
\alias{YX_ord}
\docType{data}
\title{
ordinal relational data and covariates
}
\description{
a synthetic dataset that includes ordinal relational data as well as information
 on eight covariates
}
\usage{data(YX_ord)}
\format{
  The format is:
List of 2
 $ Y: num [1:100, 1:100] NA 0 3 0 3 1 0 1 1 0 ...
 $ X: num [1:100, 1:100, 1:8] 1 1 1 1 1 1 1 1 1 1 ...
  ..- attr(*, "dimnames")=List of 3
  .. ..$ : NULL
  .. ..$ : NULL
  .. ..$ : chr [1:8] "intercept" "rgpa" "rsmoke" "cgpa" ...
}
\examples{
data(YX_ord)
t_trans(YX_ord$Y>0)
}
\keyword{datasets}
