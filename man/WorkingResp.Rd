\name{WorkingResp}
\alias{WorkingResp}

\title{
Calculate the working response  
}
\description{
Calculate the working response in the iterative least square regression 
}
\usage{
WorkingResp(y, Px, X, beta, intercept = 0)
}
\arguments{
   \item{X}{
Design matrix of the logistic model considered.
}
  \item{y}{
Binary response of the logistic model considered.
}
  \item{Px}{
The probability of the reponse to be 1}

  \item{beta}{
Vector of coefficients}
  \item{intercept}{
If there is an intercept}
}

\value{
This function returns the vector of working response.
}

\seealso{ Please read \url{https://hastie.su.domains/Papers/glmnet.pdf} for more details
}
\author{
Wencan Zhu
}
