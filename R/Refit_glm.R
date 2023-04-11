Refit_glm <-
function(X, beta_pred, y){
  X_temp <- X[, which(beta_pred != 0)]
  if(length(which(beta_pred != 0))==0){
    coef_est <- beta_pred
  }else if(is.null(ncol(X_temp))){
    mydata <- data.frame(Y = y, X_temp)
    colnames(mydata) <- c("Y", "X")
    formula <- paste0("Y~-1 +", paste0(colnames(mydata)[-which(colnames(mydata) == "Y")], collapse = " + "))
    myform <- as.formula(formula)
    mod_lm <- glm(myform, data = mydata, family = "binomial")
    coef_est <- mod_lm$coefficients
    
  } else {
    mydata <- data.frame(Y = y, as.matrix(X_temp))
    formula <- paste0("Y~-1 +", paste0(colnames(mydata)[-which(colnames(mydata) == "Y")], collapse = " + "))
    myform <- as.formula(formula)
    
    if(length(which(beta_pred != 0)) >= length(y)){
      mod_ridge <- cv.glmnet(x=as.matrix(X_temp), y=y, alpha=0, intercept=FALSE, family="binomial")
      opt_lambda <- mod_ridge$lambda[which.min(mod_ridge$cvm)]
      coef_est <- as.vector(glmnet(x=as.matrix(X), y=y, alpha=0, intercept=FALSE, family="binomial", lambda = opt_lambda)$beta)
    } else {
      mod_lm <- glm(myform, data = mydata, family = "binomial")
      coef_est <- mod_lm$coefficients
    }
  }
  
  beta_refit <- rep(0, length(beta_pred))
  beta_refit[which(beta_pred != 0)] <- coef_est
  return(beta_refit)
}
