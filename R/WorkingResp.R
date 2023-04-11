WorkingResp <-
function(y, Px, X, beta, intercept=0){
  pred <- X%*%beta+intercept
  temp <- pred+(y-Px)/(Px*(1-Px))
  temp[which(Px==1)]=pred[which(Px==1)]
  temp[which(Px==0)]=pred[which(Px==0)]
  
  return(unlist(temp))
}
