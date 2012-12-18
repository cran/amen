summary.ame <-
function(object, ...)
{ 
  fit<-object
  tmp<-cbind(apply(fit$BETA,2,mean), apply(fit$BETA,2,sd) ,
       apply(fit$BETA,2,mean)/apply(fit$BETA,2,sd) , 
       2*(1-pnorm( abs(apply(fit$BETA,2,mean)/apply(fit$BETA,2,sd)))))
  colnames(tmp)<-c("pmean","psd","z-stat","p-val") 
  cat("\nbeta:\n")
  print(round(tmp,3))

  tmp<-apply(fit$SABR,2,mean)
  Sab<-matrix(tmp[c(1,2,2,3)],2,2) ; dimnames(Sab)<-list(c("a","b"),c("a","b"))
  cat("\nSigma_ab pmean:\n") ; print(round(Sab,3))

 cat("\nrho pmean:\n", round(tmp[4],3) ,"\n") 
}
