t_recip <-
function(Y)
{
  sum(Y*t(Y), na.rm=TRUE ) / sum(Y,na.rm=TRUE)
}
