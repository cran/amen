\name{simY_ord}
\alias{simY_ord}
\title{
Simulate an ordinal relational matrix 
}
\description{
Simulates an ordinal relational matrix having a particular marginal distribution
}
\usage{
simY_ord(EZ, rho, uY, FY)
}
\arguments{
  \item{EZ}{
square matrix giving the expected value of the latent Z matrix
}
  \item{rho}{
dyadic correlation
}
  \item{uY}{
support of the ordinal outcome variable 
}
  \item{FY}{cumulative CDF for \code{uY}
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
function (EZ, rho, uY, FY) 
{
    ZS <- simZ(EZ, rho)
    diag(ZS) <- NA
    qZ <- quantile(ZS, FY[-length(FY)], na.rm = TRUE)
    YS <- ZS * 0 + max(uY)
    for (w in 1:(length(uY) - 1)) {
        YS[(YS == max(uY)) & ZS <= qZ[w]] <- uY[w]
    }
    YS
  }
}


