## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE)
library(WLogit)
library(tibble)
library(ggplot2)
set.seed(123456)

## ----generate Sigma-----------------------------------------------------------
p <- 500 # number of variables 
d <- 10 # number of actives
n <- 100 # number of samples
actives <- c(1:d)
nonacts <- c(1:p)[-actives]
Sigma <- matrix(0, p, p)
Sigma[actives, actives] <- 0.3
Sigma[-actives, actives] <- 0.5
Sigma[actives, -actives] <- 0.5
Sigma[-actives, -actives] <- 0.7
diag(Sigma) <- rep(1,p)

## ----X, eval=FALSE------------------------------------------------------------
#  X <- MASS::mvrnorm(n = n, mu=rep(0,p), Sigma, tol = 1e-6, empirical = FALSE)
#  beta <- rep(0,p)
#  beta[actives] <- 1
#  pr <- CalculPx(X,beta=beta)
#  y <- rbinom(n,1,pr)

## ---- echo=FALSE--------------------------------------------------------------
data(X)
data(y)
data(beta)

## ----WLogit model, eval=FALSE-------------------------------------------------
#  mod <- WhiteningLogit(X = X, y = y)

## ---- echo=FALSE--------------------------------------------------------------
data(test)
mod <- test

## ----variable selection,fig.width=4,fig.height=3------------------------------
beta_min <- mod$beta.min
df_beta <- data.frame(beta_est=beta_min, Status = ifelse(beta==0, "non-active", "active"))
df_plot <- df_beta[which(beta_min!=0), ]
df_plot$index <- which(beta_min!=0)
ggplot2::ggplot(data=df_plot, mapping=aes(y=beta_est, x=index, color=Status))+geom_point()+
  theme_bw()+ylab("Estimated coefficients")+xlab("Indices of selected variables")

