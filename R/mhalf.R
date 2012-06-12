mhalf <-
function(M) 
{ 
  #symmetric square  root of a pos def matrix
  tmp<-eigen(M)
  tmp$vec%*%sqrt(diag(tmp$val,nrow=nrow(M)))%*%t(tmp$vec)
}
