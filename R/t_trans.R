t_trans <-
function(Y) 
{ 
  YS<-1*(Y+t(Y)>0) 
  sm<-0 
  for(i in 1:nrow(YS))
  {
    ci<-which(YS[i,]>0) 
    sm<-sm+ sum(YS[ci,ci],na.rm=TRUE) 
  }
 sm/6
}
