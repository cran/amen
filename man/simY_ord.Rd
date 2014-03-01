\name{simY_ord}
\alias{simY_ord}
\title{
Simulate an ordinal relational matrix 
}
\description{
Simulates an ordinal relational matrix having a particular marginal distribution
}
\usage{
simY_ord(EZ, rho, Y)
}
\arguments{
  \item{EZ}{
square matrix giving the expected value of the latent Z matrix
}
\item{rho}{ scalar giving the within-dyad correlation }
  \item{Y}{ordinal relational data matrix 
}
}
\value{
a square matrix 
}
\author{
Peter Hoff
}

\examples{
## The function is currently defined as 
simY_ord<-
function (EZ, rho,Y)
{

    uY<-sort(unique(c(Y)))
    FY<-table(c(Y)) ; FY<-FY/sum(FY) ; FY<-cumsum(FY)
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


