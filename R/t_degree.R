t_degree <-
function(Y)
{
  list(od=table(c(0:nrow(Y),apply(Y,1,sum,na.rm=TRUE)))-1,
       id=table(c(0:nrow(Y),apply(Y,2,sum,na.rm=TRUE)))-1 )
}
