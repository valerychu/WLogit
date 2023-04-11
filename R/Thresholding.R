Thresholding <-
function(X, y, coef, TOP){
  beta_top <- sapply(TOP, top_thresh, vect=coef)
  auc <- c()
  opt_top=max(TOP)
  prob_temp <- CalculPx(X, beta_top[,1])
  ll <- prob_temp^y*(1-prob_temp)^(1-y)
  auc[1] <- -log(prod(ll))
  for (j in 2:length(TOP)) {
    prob_temp <- CalculPx(X, beta_top[,j])
    ll <- prob_temp^y*(1-prob_temp)^(1-y)
    ll <- ifelse(ll<0.000001, 0.1, ll)
    auc[j] <- -log(prod(ll))
    if((auc[j]==Inf) ||(auc[j-1]==Inf)){
      next
    } else{
      ratio_auc <- round(auc[j]/auc[j-1], 6)
      if(ratio_auc >= 0.9999){
        opt_top <- TOP[j]
        break
      }
    }
  }
  
  return(list(opt_top=opt_top, auc=auc))
}
