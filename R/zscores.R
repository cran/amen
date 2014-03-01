zscores<-function(y)
{
 qnorm( rank(y,na.last="keep",ties.method="average")/(1+sum(!is.na(y))) )
}

