\name{simY_nrm}
\alias{simY_nrm}
\title{
Simulate a normal relational matrix
}
\description{
Simulates a normal relational matrix
}
\usage{
simY_nrm(EY, rho, s2)
}
\arguments{
  \item{EY}{
square matrix giving the expected value of the relational matrix
}
  \item{rho}{
dyadic correlation
}
  \item{s2}{
dyaidic variance
}
}
\value{
a square matrix
}
\author{
Peter Hoff
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (EY, rho, s2) 
{
    YS <- simZ(EY, rho, s2)
    diag(YS) <- NA
    YS
  }
}
