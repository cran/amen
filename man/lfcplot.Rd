\name{lfcplot}

\alias{lfcplot}

\title{
Latent factor circle plot
}

\description{
A plot representing a two-dimensional latent factor fit to a social network
}

\usage{lfcplot(U,V,Y=NULL,ecol="lightblue",sncol="darkgreen",rncol="darkred",snlab=NULL,rnlab=NULL)}

\arguments{
\item{U}{a matrix of sender-specific latent factors with at least two columns}
\item{V}{a matrix of receiver-specific latent factors with at least two columns}
\item{Y}{a sociomatrix}
\item{ecol}{color of edges} 
\item{sncol}{color of nodes as senders (can be node-specific)}
\item{rncol}{color of nodes as receivers (ca be node-specific)}
\item{snlab}{labels of nodes as senders (character or pch values)} 
\item{rnlab}{labels of nodes as receivers (character or pch values)} 
}

\value{
a plot
}

\author{
Peter Hoff
}

\examples{
data(YX_bin)

### You should run your MCMC algorithm much longer than this
fit<-ame_bin(YX_bin$Y,R=2,burn=5,nscan=50,odens=5,plot=FALSE) 
lfcplot(fit$U,fit$V,YX_bin$Y) 
}

