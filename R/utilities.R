boxcox_trans <- function(x, lambda){
  if(lambda != 0) y <- (x^lambda - 1)/lambda
  else y <- log(x)
  return(y)
}

boxcox_trans_back <- function(x, lambda){
  if(lambda != 0) y <- (x*lambda + 1)^(1/lambda)
  else y <- exp(x)
  return(y)
}
