## are ties present
gettiesDunn <- function(x){
  n <- length(x)
  t <- table(x)
  C <- sum(t^3 - t) / (12 * (n - 1))
  return(C)
}
