\name{raSab_frn_fc}
\alias{raSab_frn_fc}
\title{
Simulate a and Sab from full conditional distributions under frn likelihood
}
\description{
Simulate a and Sab from full conditional distributions under frn likelihood
}
\usage{
raSab_frn_fc(Z, Y, YL, a, b, Sab, odobs, odmax, SS = round(sqrt(nrow(Z))))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
a square matrix, the current value of Z 
}
  \item{Y}{
square matrix of ranked nomination data
}
  \item{YL}{
list of ranked individuals, from least to most preferred in each row
}
  \item{a}{
current value of row effects
}
  \item{b}{
current value of column effects
}
  \item{Sab}{
current value of Cov(a,b)
}
  \item{odobs}{
observed outdegree
}
  \item{odmax}{
a scalar or vector giving the maximum number of nominations for each individual
}
  \item{SS}{
number of iterations
}
}
\value{
\item{Z}{new value of Z}
\item{Sab}{new value of Sab}
\item{a}{new value of a}
}
\author{
Peter Hoff
}

\examples{
## The function is currently defined as
function (Z, Y, YL, a, b, Sab, odobs, odmax, SS = round(sqrt(nrow(Z)))) 
{
    E <- Z - a \%*\% t(rep(1, nrow(Z)))
    lba <- -E[cbind(1:nrow(Z), YL[, 1])]
    lba[is.na(lba)] <- -Inf
    uba <- -apply(E - (Y != 0) * (Inf^(Y != 0)), 1, max, na.rm = TRUE)
    uba[odobs == odmax] <- Inf
    for (ss in 1:SS) {
        ea <- b * Sab[1, 2]/Sab[2, 2]
        sa <- sqrt(Sab[1, 1] - Sab[1, 2]^2/Sab[2, 2])
        a <- ea + sa * qnorm(runif(nrow(Z), pnorm((lba - ea)/sa), 
            pnorm((uba - ea)/sa)))
        Sab <- solve(rwish(solve(diag(2) + crossprod(cbind(a, 
            b))), 3 + nrow(Z)))
    }
    list(Z = E + a \%*\% t(rep(1, nrow(Z))), a = a, Sab = Sab)
  }
}

