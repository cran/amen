\name{plot.ame}
\alias{plot.ame}
\title{
Plot results of an AME object
}
\description{
A set of plots summarizing the MCMC routine for an AME fit, as well as 
some posterior predictive checks. 
}
\usage{
\method{plot}{ame}(x, ...)
}
\arguments{
  \item{x}{the result of fitting an AME model}
  \item{...}{additional parameters (not used)}
}
\value{
a series of plots 
}
\author{
Peter Hoff
}

\examples{
## The function is currently defined as
function (x, ...) 
{
  fit<-x
    par(mfrow = c(3, 2), mar = c(3, 3, 1, 1), mgp = c(1.75, 0.75, 
        0))
    mSABR <- apply(fit$SABR, 2, median)
    matplot(fit$SABR, type = "l", lty = 1, ylab = "SABR")
    abline(h = mSABR, col = 1:length(mSABR))
    mBETA <- apply(fit$BETA, 2, median)
    matplot(fit$BETA, type = "l", lty = 1, col = 1:length(mBETA), 
        ylab = "BETA")
    abline(h = mBETA, col = 1:length(mBETA))
    abline(h = 0, col = "gray")
    plot(fit$TR, ylim = range(c(fit$TR, fit$tr)), type = "l", 
        ylab = "reciprocity", xlab = "")
    abline(h = fit$tr)
    plot(fit$TT, ylim = range(c(fit$TT, fit$tt)), type = "l", 
        ylab = "transitive triples", xlab = "")
    abline(h = fit$tt)
    mod <- max(which(fit$td$od > 0)) - 1
    mid <- max(which(fit$td$id > 0)) - 1
    qod <- apply(fit$TOD, 2, quantile, prob = c(0.975, 0.5, 0.25))
    qid <- apply(fit$TID, 2, quantile, prob = c(0.975, 0.5, 0.25))
    odd <- rbind(qod, fit$td$od)[, 0:(mod + 1)]
    plot(c(0, mod), range(odd), type = "n", xlab = "outdegree", 
        ylab = "count")
    for (k in 1:4) {
        lines(0:mod, odd[k, ], col = c("gray", "black")[1 + (k == 
            4)])
    }
    idd <- rbind(qid, fit$td$id)[, 0:(mid + 1)]
    plot(c(0, mid), range(idd), type = "n", xlab = "indegree", 
        ylab = "count")
    for (k in 1:4) {
        lines(0:mid, idd[k, ], col = c("gray", "black")[1 + (k == 
            4)])
    }
  }
}
