## -----------------------------------------------------------------------------
library(amen)

## -----------------------------------------------------------------------------
data(lazegalaw)

Y<-lazegalaw$Y[,,2]
Xn<-lazegalaw$X[,c(2,4,5,6)]
Xd<-lazegalaw$Y[,,-2]
Xd<-array( c(Xd,outer(Xn[,4],Xn[,4],"==")),dim=dim(Xd)+c(0,0,1))
dimnames(Xd)[[3]]<-c("advice","cowork","samepractice")

dimnames(Xd)[[3]]
dimnames(Xn)[[2]]

## ----results='hide',message=FALSE---------------------------------------------
netplot(lazegalaw$Y[,,2],ncol=Xn[,4])

## ----fig.keep='last',results='hide',cache=TRUE--------------------------------
fitSRRM<-ame(Y, Xd=Xd, Xr=Xn, Xc=Xn, family="bin")

## -----------------------------------------------------------------------------
summary(fitSRRM) 

## ----fig.keep='last',results='hide',cache=TRUE--------------------------------
fitAME<-ame(Y, Xd=Xd, Xr=Xn, Xc=Xn, R=3, family="bin")

## -----------------------------------------------------------------------------
summary(fitAME) 

