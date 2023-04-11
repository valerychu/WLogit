\name{Refit_glm}
\alias{Refit_glm}

\title{
Refit the logistic regression with chosen variables
}
\description{
Refit the logistic regression with chosen variables.
}
\usage{
Refit_glm(X, beta_pred, y)
}
\arguments{
  \item{X}{
Design matrix of the logistic model considered.
}
  \item{beta_pred}{
Predicted coefficients to be refited.
}
  \item{y}{
Binary response
}
}

\value{
\item{beta_refit }{The new estimated coefficients}
}

\author{
Wencan Zhu
}

\examples{
data(X)
data(y)
data(beta)
Refit_glm(X=X, beta_pred=beta, y=y)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, beta_pred, y) 
{
    X_temp <- X[, which(beta_pred != 0)]
    if (length(which(beta_pred != 0)) == 0) {
        coef_est <- beta_pred
    }
    else if (is.null(ncol(X_temp))) {
        mydata <- data.frame(Y = y, X_temp)
        colnames(mydata) <- c("Y", "X")
        formula <- paste0("Y~-1 +", paste0(colnames(mydata)[-which(colnames(mydata) == 
            "Y")], collapse = " + "))
        myform <- as.formula(formula)
        mod_lm <- glm(myform, data = mydata, family = "binomial")
        coef_est <- mod_lm$coefficients
    }
    else {
        mydata <- data.frame(Y = y, as.matrix(X_temp))
        formula <- paste0("Y~-1 +", paste0(colnames(mydata)[-which(colnames(mydata) == 
            "Y")], collapse = " + "))
        myform <- as.formula(formula)
        if (length(which(beta_pred != 0)) >= length(y)) {
            mod_ridge <- cv.glmnet(x = as.matrix(X_temp), y = y, 
                alpha = 0, intercept = FALSE, family = "binomial")
            opt_lambda <- mod_ridge$lambda[which.min(mod_ridge$cvm)]
            coef_est <- as.vector(glmnet(x = as.matrix(X), y = y, 
                alpha = 0, intercept = FALSE, family = "binomial", 
                lambda = opt_lambda)$beta)
        }
        else {
            mod_lm <- glm(myform, data = mydata, family = "binomial")
            coef_est <- mod_lm$coefficients
        }
    }
    beta_refit <- rep(0, length(beta_pred))
    beta_refit[which(beta_pred != 0)] <- coef_est
    return(beta_refit)
  }
}
