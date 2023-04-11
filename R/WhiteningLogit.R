WhiteningLogit <-
function(X=X, y=y,
                   nlambda=50,
                   maxit=100,
                   gamma=0.9999,
         top_grill=c(1:100)){
  
  p=ncol(X)
  n=nrow(X)
  
  mod_ridge <- cv.glmnet(x=as.matrix(X), y=y, alpha=0.5, intercept=FALSE, family="binomial")
  pr_est <- predict(mod_ridge, as.matrix(X), s = "lambda.min", type="response")
  beta_ini <- predict(mod_ridge, as.matrix(X), s = "lambda.min", type="coefficients")[-1]
  diag_w <- pr_est*(1-pr_est)
  square_root_w <- diag(sqrt(as.vector(diag_w)), nrow=n)
  X_new <- square_root_w%*%X
  
  Cov_est <- cvCovEst(
    dat = X_new,
    estimators = c(
      linearShrinkLWEst, thresholdingEst, sampleCovEst
    ),
    estimator_params = list(
      thresholdingEst = list(gamma = seq(0.1, 0.3, 0.1))
    ),
    center = TRUE,
    scale = TRUE
  )
  
  Sigma_est <- Cov_est$estimate
  
  SVD_new <- fast.svd(Sigma_est)
  U_sigma_new <- SVD_new$u
  D_sigma_new <- SVD_new$d
  inv_transmat <- U_sigma_new%*%diag(sqrt(D_sigma_new))%*%t(U_sigma_new)
  inv_diag_new <- ifelse(D_sigma_new<0.000001, 0, 1/sqrt(D_sigma_new))
  trans_mat <- U_sigma_new%*%diag(inv_diag_new)%*%t(U_sigma_new)
  
  
  if (p <= 50) {
    top_grill <- seq(1, p, 2)
  }else if (p <= 200) {
    top_grill <- c(1:50, seq(52, p, 2))
  }else if (p <= 300) {
    top_grill <- c(1:50, seq(52, 100, 2), seq(105, 200, 5), 
                   seq(210, p, 10))
  }else {
    top_grill <- c(1:50, seq(52, 100, 2), seq(105, 200, 5), 
                   seq(210, 300, 10))
  }
  
  X_tilde <- X%*%trans_mat
  
  beta_tilde_ini <-  inv_transmat%*%beta_ini
  Px <- CalculPx(X_tilde, beta=beta_tilde_ini)
  wt <- CalculWeight(Px)
  # wt <- ifelse(wt0==0, 0.0001, wt0)
  ystar <- WorkingResp(y=y, Px=Px, X=X_tilde, beta=beta_tilde_ini)
  X_tilde_weighted <- sweep(X, MARGIN=1, sqrt(wt), `*`)
  ystar_weighted <- sqrt(wt)*ystar
  
  gen.model0 <- genlasso(y=ystar_weighted, X=X_tilde_weighted, 
                         D=trans_mat, maxsteps = 50)
  parameter_tmp <- beta_tilde_ini
  beta_final <- matrix(NA, length(gen.model0$lambda), p)
  skip_i <- TRUE
  eval_final <- c()
  defaultW <- getOption("warn") 

  for(i in 1:length(gen.model0$lambda)){
    #inner loop
    epsilon=10
    j=0
    if(skip_i){parameter_tmp <- beta_tilde_ini
    } else {parameter_tmp <- parameter_current}
    skip_i <-FALSE
    
    while(epsilon > 0.001){
      j=j+1
      parameter_current <- parameter_tmp
      Px <- CalculPx(X_tilde, beta=parameter_current)
      wt0 <- CalculWeight(Px)
      wt <- ifelse(round(wt0,4)==0, 0.0001, wt)
      ystar <- WorkingResp(y=y, Px=Px, X=X_tilde, beta=parameter_current)
      X_tilde_weighted <- sweep(X, MARGIN=1, sqrt(wt), `*`)
      ystar_weighted <- sqrt(wt)*ystar
      
      gen.model <- genlasso(y=ystar_weighted, X=X_tilde_weighted, D=trans_mat, maxsteps =   maxit)
      
      if(gen.model0$lambda[i] < min(gen.model$lambda)){
        parameter_tmp <- parameter_current
        break
      } else {
        parameter_tmp <- coef(gen.model, lambda=gen.model0$lambda[i],
                              type = "primal")$beta
        beta_current <- parameter_tmp
        if(sum(is.na(parameter_tmp))>0){
          skip_i <-TRUE 
          parameter_tmp <- rep(0,p)
          break}
        epsilon <- max(abs(parameter_current-parameter_tmp))
        if(epsilon >=100){
          skip_i <-TRUE 
          break}
        if (j==maxit){
          skip_i <-TRUE 
          break}
      }
    }
    
    if(skip_i){
      beta_final[i, ] <- rep(NA, p)
      eval_final[i] <- NA
    } else{
      
      correction <- Thresholding(X_tilde, y, coef=parameter_tmp, TOP=top_grill)
      opt_top_tilde <- correction$opt_top
      beta_tilde_opt <- top_thresh(vect=parameter_tmp, thresh = opt_top_tilde)
      beta_final0 <- trans_mat%*%beta_tilde_opt
      
      correction <- Thresholding(X, y, coef=beta_final0, TOP=top_grill)
      opt_top_final <- correction$opt_top
      beta_final[i, ] <- beta_opt_final <- top(vect=beta_final0, thresh = opt_top_final)
      
      beta_refit <- Refit_glm(X=X, beta_pred = beta_opt_final, y=y)
      pr_est <- CalculPx(X, beta_refit)
      ll <- pr_est^y*(1-pr_est)^(1-y)
      #ll <- ifelse(ll<0.000001, 1, ll)
      eval_final[i] <- -log(prod(ll))
      
    }
    beta.min <- beta_final[which.min(eval_final), ]
  }
  options(warn = defaultW)
  return(list(beta=beta_final, lambda=gen.model0$lambda, beta.min=beta.min, 
              log.likelihood=eval_final))
}
