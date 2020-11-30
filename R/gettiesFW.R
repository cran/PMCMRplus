## for Fligner-Wolfe
gettiesFW <- function(x){
  n <- length(x)
  t <- table(x)
  C <- sum(t * (t^2 - 1))
  return(C)
}
