CalculPx <-
function(X, beta, intercept=0){
  prob <- 1/(1+exp(-(X%*%beta+intercept)))
  return(prob)
}
